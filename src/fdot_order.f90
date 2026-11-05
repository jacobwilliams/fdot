module fdot_order
!! Ordering the nodes of each rank to reduce crossings, with virtual nodes for edges that span several
!! ranks. This is a port of `lib/dotgen/mincross.c` and `class2.c` of Graphviz, which are under the
!! Eclipse Public License 2.0, without clusters and flat edges.
!!
!! Flat edges come from `rank=same`, which is not supported. An edge between two different nodes has a
!! length of at least 1 in the graph that is ranked, so it can never join nodes of one rank, and
!! the `stat` of [[dot_mincross]] reports it if that ever turns out to be false.
use iso_fortran_env, only: int64
use fdot_structure
implicit none
private

public :: dot_mincross

integer, parameter :: MC_SCALE = 256            !! the scale of medians
integer, parameter :: MIN_QUIT = 8              !! iterations without improvement before stopping
integer, parameter :: MAX_ITER = 24             !! the limit of iterations of a pass
real(wp), parameter :: CONVERGENCE = 0.995_wp   !! what counts as an improvement

contains

!> Order the ranks. `stat` is not 0 if a flat edge turns up.
subroutine dot_mincross(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    integer(int64) :: mc
    integer :: k, r, i, v, total

    call class2(c, stat)
    if (stat /= 0) return
    call decompose(c)
    call allocate_ranks(c)
    do k = 1, size(c%comp)
        if (k > 1) then
            do r = c%minrank, c%maxrank
                c%rk(r)%base = c%rk(r)%base + c%rk(r)%n
                c%rk(r)%n = 0
            end do
        end if
        call mincross(c, k, mc)
    end do

    ! one list of all nodes, and ranks holding all components
    total = sum(c%comp%n)
    allocate (c%nlist(total))
    i = 0
    do k = 1, size(c%comp)
        c%nlist(i + 1:i + c%comp(k)%n) = c%comp(k)%a(1:c%comp(k)%n)
        i = i + c%comp(k)%n
    end do
    do r = c%minrank, c%maxrank
        c%rk(r)%n = c%rk(r)%an
        c%rk(r)%base = 0
        do i = 0, c%rk(r)%n - 1
            v = c%rk(r)%v(i)
            if (v == 0) then
                c%rk(r)%n = i
                exit
            end if
            c%nd(v)%order = i
        end do
    end do
end subroutine dot_mincross

!> Build the graph of ranks: virtual nodes for long edges, and a list of the others.
subroutine class2(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    integer :: n, i, e, prev, t, h, opp, k
    logical :: found

    stat = 0
    do n = 1, c%nn
        do i = 1, c%gout(n)%n
            e = c%gout(n)%a(i)
            associate (hd => c%nd(c%ed(e)%head), tl => c%nd(c%ed(e)%tail))
                if (hd%weight_class <= 2) hd%weight_class = hd%weight_class + 1
                if (tl%weight_class <= 2) tl%weight_class = tl%weight_class + 1
            end associate
        end do
    end do

    do n = 1, c%nn
        prev = 0
        do i = 1, c%gout(n)%n
            e = c%gout(n)%a(i)
            if (c%ed(e)%to_virt /= 0) then
                prev = e
                cycle
            end if
            t = c%ed(e)%tail
            h = c%ed(e)%head

            ! multi-edges
            if (prev /= 0) then
                if (t == c%ed(prev)%tail .and. h == c%ed(prev)%head) then
                    if (c%nd(t)%rank == c%nd(h)%rank) then
                        call merge_oneway(c, e, prev)
                        call other_edge(c, e)
                        cycle
                    end if
                    if (.not. c%labelled(e) .and. .not. c%labelled(prev)) then
                        if (c%concentrate) then
                            c%ed(e)%etype = IGNORED
                        else
                            call merge_chain(c, e, c%ed(prev)%to_virt, .true.)
                            call other_edge(c, e)
                        end if
                        cycle
                    end if
                end if
            end if

            ! self-loops
            if (t == h) then
                call other_edge(c, e)
                prev = e
                cycle
            end if

            ! flat edges
            if (c%nd(t)%rank == c%nd(h)%rank) then
                stat = 3
                return
            end if

            ! forward edges
            if (c%nd(h)%rank > c%nd(t)%rank) then
                call make_chain(c, t, h, e)
                prev = e
                cycle
            end if

            ! backward edges
            found = .false.
            do k = 1, c%gout(h)%n
                opp = c%gout(h)%a(k)
                if (c%ed(opp)%head /= t .or. c%ed(opp)%head == h .or. c%ed(opp)%etype == IGNORED) cycle
                ! shadows a forward edge
                if (c%ed(opp)%to_virt == 0) call make_chain(c, c%ed(opp)%tail, c%ed(opp)%head, opp)
                if (.not. c%labelled(e) .and. .not. c%labelled(opp)) then
                    if (c%concentrate) then
                        c%ed(e)%etype = IGNORED
                        c%ed(opp)%conc_opp = .true.
                    else
                        call other_edge(c, e)
                        call merge_chain(c, e, c%ed(opp)%to_virt, .true.)
                    end if
                    found = .true.
                    exit
                end if
            end do
            if (found) cycle
            call make_chain(c, h, t, e)
            prev = e
        end do
    end do
end subroutine class2

!> Create the chain of virtual nodes for edge `orig`.
subroutine make_chain(c, from, to, orig)
    type(core_t), intent(inout) :: c    !! the state
    integer, value :: from              !! the end with the lower rank
    integer, value :: to                !! the end with the higher rank
    integer, value :: orig              !! the edge
    integer :: r, label_rank, u, v, e, last

    u = from
    label_rank = -1
    if (c%labelled(orig)) label_rank = (c%nd(from)%rank + c%nd(to)%rank)/2
    last = c%nd(to)%rank
    do r = c%nd(from)%rank + 1, last
        if (r < last) then
            call virtual_node(c, v)
            if (r == label_rank) then
                c%nd(v)%label = orig
                c%nd(v)%lw = NODESEP
                if (c%flip) then
                    c%nd(v)%ht = c%label_w(orig)
                    c%nd(v)%rw = c%label_h(orig)
                else
                    c%nd(v)%ht = c%label_h(orig)
                    c%nd(v)%rw = c%label_w(orig)
                end if
            else
                c%nd(v)%lw = c%nd(v)%lw + NODESEP/2
                c%nd(v)%rw = c%nd(v)%rw + NODESEP/2
            end if
            c%nd(v)%rank = r
        else
            v = to
        end if
        call virtual_edge(c, u, v, orig, e)
        call virtual_weight(c, e)
        u = v
    end do
end subroutine make_chain

!> Multiply the weight of `e`, by the kinds of its ends.
subroutine virtual_weight(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    integer, parameter :: TABLE(0:2, 0:2) = &  !! the factor of a weight, by the kinds of its ends
        reshape([1, 1, 1, 1, 2, 2, 1, 2, 4], [3, 3])

    c%ed(e)%weight = c%ed(e)%weight*TABLE(endpoint_class(c%nd(c%ed(e)%head)), endpoint_class(c%nd(c%ed(e)%tail)))
end subroutine virtual_weight

!> The kind of an end of an edge: 0 for ordinary, 1 for a node with one edge, 2 for virtual.
integer function endpoint_class(n)
    type(vnode_t), intent(in) :: n  !! the node

    if (n%kind == VIRTUAL) then
        endpoint_class = 2
    else if (n%weight_class <= 1) then
        endpoint_class = 1
    else
        endpoint_class = 0
    end if
end function endpoint_class

!> Make room in each rank.
subroutine allocate_ranks(c)
    type(core_t), intent(inout) :: c    !! the state
    integer, allocatable :: cn(:)
    integer :: n, i, e, low, high, r

    allocate (cn(0:c%maxrank + 1), source=0)
    do n = 1, c%nn
        cn(c%nd(n)%rank) = cn(c%nd(n)%rank) + 1
        do i = 1, c%gout(n)%n
            e = c%gout(n)%a(i)
            low = min(c%nd(c%ed(e)%tail)%rank, c%nd(c%ed(e)%head)%rank)
            high = max(c%nd(c%ed(e)%tail)%rank, c%nd(c%ed(e)%head)%rank)
            cn(low + 1:high - 1) = cn(low + 1:high - 1) + 1
        end do
    end do
    if (allocated(c%rk)) deallocate (c%rk)
    allocate (c%rk(0:c%maxrank + 1))
    do r = 0, c%maxrank + 1
        c%rk(r)%an = cn(r) + 1
        c%rk(r)%n = cn(r) + 1
        c%rk(r)%base = 0
        allocate (c%rk(r)%v(0:cn(r) + 1), source=0)
    end do
    c%rk(c%maxrank + 1)%n = 0
end subroutine allocate_ranks

!> Node `i` of rank `r`.
integer function at(c, r, i)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: r        !! the rank
    integer, intent(in) :: i        !! the position

    at = c%rk(r)%v(c%rk(r)%base + i)
end function at

!> Install node `n` at the right end of its rank.
subroutine install_in_rank(c, n)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: n            !! the node
    integer :: r

    r = c%nd(n)%rank
    c%rk(r)%v(c%rk(r)%base + c%rk(r)%n) = n
    c%nd(n)%order = c%rk(r)%n
    c%rk(r)%n = c%rk(r)%n + 1
end subroutine install_in_rank

!> Install the nodes of component `k` in their ranks. The search takes the better of searching from
!> the sources and from the sinks, so that trees are drawn with no crossings.
subroutine build_ranks(c, k, pass)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: k            !! the component
    integer, intent(in) :: pass         !! 0 to search from the sources, 1 from the sinks
    logical, allocatable :: mark(:)
    integer, allocatable :: q(:)
    integer :: i, j, n, n0, first, last, r, e, other, count

    allocate (mark(c%nnodes), source=.false.)
    allocate (q(c%comp(k)%n))
    do r = c%minrank, c%maxrank
        c%rk(r)%n = 0
    end do
    do i = 1, c%comp(k)%n
        n = c%comp(k)%a(i)
        if (pass == 0) then
            if (c%nd(n)%in%n > 0) cycle
        else
            if (c%nd(n)%out%n > 0) cycle
        end if
        if (mark(n)) cycle
        mark(n) = .true.
        first = 1
        last = 1
        q(1) = n
        do while (first <= last)
            n0 = q(first)
            first = first + 1
            call install_in_rank(c, n0)
            if (pass == 0) then
                count = c%nd(n0)%out%n
            else
                count = c%nd(n0)%in%n
            end if
            do j = 1, count
                if (pass == 0) then
                    e = c%nd(n0)%out%a(j)
                    other = c%ed(e)%head
                else
                    e = c%nd(n0)%in%a(j)
                    other = c%ed(e)%tail
                end if
                if (.not. mark(other)) then
                    mark(other) = .true.
                    last = last + 1
                    q(last) = other
                end if
            end do
        end do
    end do
    do r = c%minrank, c%maxrank
        c%rk(r)%valid = .false.
        if (c%flip .and. c%rk(r)%n > 0) then
            do j = 0, (c%rk(r)%n - 1)/2
                call exchange(c, at(c, r, j), at(c, r, c%rk(r)%n - 1 - j))
            end do
        end if
    end do
    if (ncross(c) > 0) call transpose(c, .false.)
end subroutine build_ranks

!> Exchange two nodes of a rank.
subroutine exchange(c, v, w)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: v            !! one node
    integer, intent(in) :: w            !! the other node
    integer :: vi, wi, r

    r = c%nd(v)%rank
    vi = c%nd(v)%order
    wi = c%nd(w)%order
    c%nd(v)%order = wi
    c%rk(r)%v(c%rk(r)%base + wi) = v
    c%nd(w)%order = vi
    c%rk(r)%v(c%rk(r)%base + vi) = w
end subroutine exchange

!> The crossings among the in-edges of `v` and `w`, with `v` to the left.
integer(int64) function in_cross(c, v, w) result(cross)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: v        !! the left node
    integer, intent(in) :: w        !! the right node
    integer :: i, j, e1, e2, inv

    cross = 0
    do j = 1, c%nd(w)%in%n
        e2 = c%nd(w)%in%a(j)
        inv = c%nd(c%ed(e2)%tail)%order
        do i = 1, c%nd(v)%in%n
            e1 = c%nd(v)%in%a(i)
            if (c%nd(c%ed(e1)%tail)%order - inv > 0) then
                cross = cross + int(c%ed(e1)%xpenalty, int64)*c%ed(e2)%xpenalty
            end if
        end do
    end do
end function in_cross

!> The crossings among the out-edges of `v` and `w`, with `v` to the left.
integer(int64) function out_cross(c, v, w) result(cross)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: v        !! the left node
    integer, intent(in) :: w        !! the right node
    integer :: i, j, e1, e2, inv

    cross = 0
    do j = 1, c%nd(w)%out%n
        e2 = c%nd(w)%out%a(j)
        inv = c%nd(c%ed(e2)%head)%order
        do i = 1, c%nd(v)%out%n
            e1 = c%nd(v)%out%a(i)
            if (c%nd(c%ed(e1)%head)%order - inv > 0) then
                cross = cross + int(c%ed(e1)%xpenalty, int64)*c%ed(e2)%xpenalty
            end if
        end do
    end do
end function out_cross

!> Try exchanging neighbours in rank `r`, and return how many crossings this removes.
subroutine transpose_step(c, r, reverse, rv)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: r            !! the rank
    logical, intent(in) :: reverse      !! whether to exchange when it makes no difference
    integer(int64), intent(out) :: rv   !! the result
    integer(int64) :: c0, c1
    integer :: i, v, w

    rv = 0
    c%rk(r)%candidate = .false.
    do i = 0, c%rk(r)%n - 2
        v = at(c, r, i)
        w = at(c, r, i + 1)
        c0 = 0
        c1 = 0
        if (r > 0) then
            c0 = c0 + in_cross(c, v, w)
            c1 = c1 + in_cross(c, w, v)
        end if
        if (c%rk(r + 1)%n > 0) then
            c0 = c0 + out_cross(c, v, w)
            c1 = c1 + out_cross(c, w, v)
        end if
        if (c1 < c0 .or. (c0 > 0 .and. reverse .and. c1 == c0)) then
            call exchange(c, v, w)
            rv = rv + c0 - c1
            c%rk(r)%valid = .false.
            c%rk(r)%candidate = .true.
            if (r > c%minrank) then
                c%rk(r - 1)%valid = .false.
                c%rk(r - 1)%candidate = .true.
            end if
            if (r < c%maxrank) then
                c%rk(r + 1)%valid = .false.
                c%rk(r + 1)%candidate = .true.
            end if
        end if
    end do
end subroutine transpose_step

!> Exchange neighbours until that no longer helps.
subroutine transpose(c, reverse)
    type(core_t), intent(inout) :: c    !! the state
    logical, intent(in) :: reverse      !! whether to exchange when it makes no difference
    integer(int64) :: delta, step
    integer :: r

    c%rk(c%minrank:c%maxrank)%candidate = .true.
    do
        delta = 0
        do r = c%minrank, c%maxrank
            if (c%rk(r)%candidate) then
                call transpose_step(c, r, reverse, step)
                delta = delta + step
            end if
        end do
        if (delta < 1) exit
    end do
end subroutine transpose

!> Reduce the crossings of component `k`, and return how many are left.
subroutine mincross(c, k, best_cross)
    type(core_t), intent(inout) :: c            !! the state
    integer, intent(in) :: k                    !! the component
    integer(int64), intent(out) :: best_cross   !! the result
    integer(int64) :: cur_cross
    integer :: pass, iter, trying, maxthispass

    cur_cross = huge(cur_cross)
    best_cross = huge(best_cross)
    do pass = 0, 2
        if (pass <= 1) then
            maxthispass = min(4, MAX_ITER)
            call build_ranks(c, k, pass)
            cur_cross = ncross(c)
            if (cur_cross <= best_cross) then
                call save_best(c)
                best_cross = cur_cross
            end if
        else
            maxthispass = MAX_ITER
            if (cur_cross > best_cross) call restore_best(c)
            cur_cross = best_cross
        end if
        trying = 0
        do iter = 0, maxthispass - 1
            if (trying >= MIN_QUIT) exit
            trying = trying + 1
            if (cur_cross == 0) exit
            call mincross_step(c, iter)
            cur_cross = ncross(c)
            if (cur_cross <= best_cross) then
                call save_best(c)
                if (real(cur_cross, wp) < CONVERGENCE*real(best_cross, wp)) trying = 0
                best_cross = cur_cross
            end if
        end do
        if (cur_cross == 0) exit
    end do
    if (cur_cross > best_cross) call restore_best(c)
    if (best_cross > 0) then
        call transpose(c, .false.)
        best_cross = ncross(c)
    end if
end subroutine mincross

!> Keep the current order.
subroutine save_best(c)
    type(core_t), intent(inout) :: c    !! the state
    integer :: r, i, n

    do r = c%minrank, c%maxrank
        do i = 0, c%rk(r)%n - 1
            n = at(c, r, i)
            c%nd(n)%save_order = c%nd(n)%order
        end do
    end do
end subroutine save_best

!> Go back to the order that was kept.
subroutine restore_best(c)
    type(core_t), intent(inout) :: c    !! the state
    integer, allocatable :: nodes(:)
    integer :: r, i, n

    do r = c%minrank, c%maxrank
        nodes = c%rk(r)%v(c%rk(r)%base:c%rk(r)%base + c%rk(r)%n - 1)
        do i = 1, size(nodes)
            n = nodes(i)
            c%nd(n)%order = c%nd(n)%save_order
            c%rk(r)%v(c%rk(r)%base + c%nd(n)%order) = n
        end do
        c%rk(r)%valid = .false.
    end do
end subroutine restore_best

!> Reorder rank `r` by the medians, with bubble passes.
subroutine reorder(c, r, reverse, hasfixed)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: r            !! the rank
    logical, intent(in) :: reverse      !! whether to exchange equal medians
    logical, intent(in) :: hasfixed     !! whether some node must stay where it is
    integer :: nelt, lp, rp, ep, changed
    real(wp) :: p1, p2

    changed = 0
    ep = c%rk(r)%n
    do nelt = c%rk(r)%n - 1, 0, -1
        lp = 0
        do while (lp < ep)
            ! find the leftmost node that can be compared
            do while (lp < ep)
                if (c%nd(at(c, r, lp))%mval >= 0) exit
                lp = lp + 1
            end do
            if (lp >= ep) exit
            ! find the next node that can be compared
            rp = lp + 1
            do while (rp < ep)
                if (c%nd(at(c, r, rp))%mval >= 0) exit
                rp = rp + 1
            end do
            if (rp >= ep) exit
            p1 = c%nd(at(c, r, lp))%mval
            p2 = c%nd(at(c, r, rp))%mval
            if (p1 > p2 .or. (p1 >= p2 .and. reverse)) then
                call exchange(c, at(c, r, lp), at(c, r, rp))
                changed = changed + 1
            end if
            lp = rp
        end do
        if (.not. hasfixed .and. .not. reverse) ep = ep - 1
    end do
    if (changed > 0) then
        c%rk(r)%valid = .false.
        if (r > 0) c%rk(r - 1)%valid = .false.
    end if
end subroutine reorder

!> One pass over the ranks.
subroutine mincross_step(c, pass)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: pass         !! the number of the pass
    integer :: r, first, last, dir
    logical :: reverse, hasfixed

    reverse = mod(pass, 4) < 2
    if (mod(pass, 2) == 0) then
        first = c%minrank + 1
        last = c%maxrank
        dir = 1
    else
        first = c%maxrank - 1
        last = c%minrank
        dir = -1
    end if
    do r = first, last, dir
        call medians(c, r, r - dir, hasfixed)
        call reorder(c, r, reverse, hasfixed)
    end do
    call transpose(c, .not. reverse)
end subroutine mincross_step

!> The crossings between ranks `r` and `r + 1`.
function rcross(c, r) result(cross)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: r        !! the rank
    integer(int64) :: cross
    integer, allocatable :: count(:)
    integer :: top, i, k, e, maxi, inv, v

    allocate (count(0:c%rk(r + 1)%n + 1), source=0)
    cross = 0
    maxi = 0
    do top = 0, c%rk(r)%n - 1
        v = at(c, r, top)
        if (maxi > 0) then
            do i = 1, c%nd(v)%out%n
                e = c%nd(v)%out%a(i)
                do k = c%nd(c%ed(e)%head)%order + 1, maxi
                    cross = cross + int(count(k), int64)*c%ed(e)%xpenalty
                end do
            end do
        end if
        do i = 1, c%nd(v)%out%n
            e = c%nd(v)%out%a(i)
            inv = c%nd(c%ed(e)%head)%order
            maxi = max(maxi, inv)
            count(inv) = count(inv) + c%ed(e)%xpenalty
        end do
    end do
end function rcross

!> All crossings, keeping the count of each rank until it changes.
function ncross(c) result(count)
    type(core_t), intent(inout) :: c    !! the state
    integer(int64) :: count
    integer :: r

    count = 0
    do r = c%minrank, c%maxrank - 1
        if (.not. c%rk(r)%valid) then
            c%rk(r)%cache_nc = rcross(c, r)
            c%rk(r)%valid = .true.
        end if
        count = count + c%rk(r)%cache_nc
    end do
end function ncross

!> Set the medians of rank `r0`, from the neighbours in rank `r1`.
subroutine medians(c, r0, r1, hasfixed)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: r0           !! the rank
    integer, intent(in) :: r1           !! the other rank
    logical, intent(out) :: hasfixed    !! always false, without flat edges
    integer, allocatable :: list(:)
    integer :: i, j, j0, n, e, rm, lm, rspan, lspan

    hasfixed = .false.
    allocate (list(16))
    do i = 0, c%rk(r0)%n - 1
        n = at(c, r0, i)
        j = 0
        if (size(list) < max(c%nd(n)%out%n, c%nd(n)%in%n)) then
            deallocate (list)
            allocate (list(2*max(c%nd(n)%out%n, c%nd(n)%in%n)))
        end if
        if (r1 > r0) then
            do j0 = 1, c%nd(n)%out%n
                e = c%nd(n)%out%a(j0)
                if (c%ed(e)%xpenalty > 0) then
                    j = j + 1
                    list(j) = MC_SCALE*c%nd(c%ed(e)%head)%order
                end if
            end do
        else
            do j0 = 1, c%nd(n)%in%n
                e = c%nd(n)%in%a(j0)
                if (c%ed(e)%xpenalty > 0) then
                    j = j + 1
                    list(j) = MC_SCALE*c%nd(c%ed(e)%tail)%order
                end if
            end do
        end if
        select case (j)
        case (0)
            c%nd(n)%mval = -1
        case (1)
            c%nd(n)%mval = list(1)
        case (2)
            c%nd(n)%mval = (list(1) + list(2))/2
        case default
            call sort_ints(list, j)
            if (mod(j, 2) == 1) then
                c%nd(n)%mval = list(j/2 + 1)
            else
                rm = j/2 + 1
                lm = rm - 1
                rspan = list(j) - list(rm)
                lspan = list(lm) - list(1)
                if (lspan == rspan) then
                    c%nd(n)%mval = (list(lm) + list(rm))/2
                else
                    c%nd(n)%mval = (list(lm)*real(rspan, wp) + list(rm)*real(lspan, wp))/(lspan + rspan)
                end if
            end if
        end select
    end do
end subroutine medians

end module fdot_order
