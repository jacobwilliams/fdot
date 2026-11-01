module fdot_rank
!! Assigning ranks: cycle breaking, connected components and network simplex, with the ranks of edge
!! labels doubled. This is a port of `lib/dotgen/rank.c`, `acyclic.c` and `class1.c` of Graphviz, which
!! are under the Eclipse Public License 2.0, without clusters and rank sets.
use fdot_structure
use fdot_simplex, only: ns_t, ns_allocate, ns_solve
implicit none
private

public :: dot_rank

contains

!> Rank the nodes of the drawing. `stat` is not 0 if network simplex fails.
subroutine dot_rank(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    integer :: e, v, i

    if (c%has_labels) then
        do e = 1, c%ne
            c%ed(e)%minlen = c%ed(e)%minlen*2
        end do
        c%ranksep = (c%ranksep + 1)/2
    end if
    call class1(c)
    call decompose(c)
    call acyclic(c)
    call rank1(c, stat)
    if (stat /= 0) return

    c%minrank = huge(c%minrank)
    c%maxrank = -1
    do v = 1, c%nn
        c%maxrank = max(c%maxrank, c%nd(v)%rank)
        c%minrank = min(c%minrank, c%nd(v)%rank)
    end do

    ! drop the graph of this stage
    do v = 1, c%nn
        c%nd(v)%in%n = 0
        c%nd(v)%out%n = 0
    end do
    do i = 1, c%ne
        c%ed(i)%to_virt = 0
    end do
    c%nedges = c%ne
end subroutine dot_rank

!> Build the graph for ranking, with a single edge for edges between the same nodes.
subroutine class1(c)
    type(core_t), intent(inout) :: c    !! the state
    integer :: n, i, e, t, h, rep, ve

    do n = 1, c%nn
        do i = 1, c%gout(n)%n
            e = c%gout(n)%a(i)
            if (c%ed(e)%to_virt /= 0) cycle
            t = c%ed(e)%tail
            h = c%ed(e)%head
            if (t == h) cycle
            rep = find_fast_edge(c, t, h)
            if (rep /= 0) then
                call merge_oneway(c, e, rep)
            else
                call virtual_edge(c, t, h, e, ve)
            end if
        end do
    end do
end subroutine class1

!> Reverse edge `e`.
subroutine reverse_edge(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    integer :: f, ve

    call delete_fast_edge(c, e)
    f = find_fast_edge(c, c%ed(e)%head, c%ed(e)%tail)
    if (f /= 0) then
        call merge_oneway(c, e, f)
    else
        call virtual_edge(c, c%ed(e)%head, c%ed(e)%tail, e, ve)
    end if
end subroutine reverse_edge

!> Break cycles, by reversing the edges that go back in a depth-first search of each component. The
!> search uses an explicit stack, as a recursive one could overflow the stack of a thread.
subroutine acyclic(c)
    type(core_t), intent(inout) :: c    !! the state
    logical, allocatable :: mark(:), onstack(:)
    integer, allocatable :: node(:), next(:)
    integer :: k, j, n, sp, e, w

    allocate (mark(c%nnodes), onstack(c%nnodes), source=.false.)
    allocate (node(c%nnodes + 1), next(c%nnodes + 1))
    do k = 1, size(c%comp)
        do j = 1, c%comp(k)%n
            n = c%comp(k)%a(j)
            if (mark(n)) cycle
            mark(n) = .true.
            onstack(n) = .true.
            sp = 1
            node(1) = n
            next(1) = 1
            do while (sp > 0)
                n = node(sp)
                if (next(sp) > c%nd(n)%out%n) then
                    onstack(n) = .false.
                    sp = sp - 1
                    cycle
                end if
                e = c%nd(n)%out%a(next(sp))
                w = c%ed(e)%head
                if (onstack(w)) then
                    call reverse_edge(c, e)
                else
                    next(sp) = next(sp) + 1
                    if (.not. mark(w)) then
                        mark(w) = .true.
                        onstack(w) = .true.
                        sp = sp + 1
                        node(sp) = w
                        next(sp) = 1
                    end if
                end if
            end do
        end do
    end do
end subroutine acyclic

!> Run network simplex on each component.
subroutine rank1(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    type(ns_t) :: s
    integer, allocatable :: local(:), edge_id(:)
    integer :: k, j, i, n, m, v, e, at

    stat = 0
    allocate (local(c%nnodes), source=0)
    allocate (edge_id(c%nedges), source=0)
    do k = 1, size(c%comp)
        associate (nodes => c%comp(k))
            m = 0
            do j = 1, nodes%n
                v = nodes%a(j)
                local(v) = j
                do i = 1, c%nd(v)%out%n
                    m = m + 1
                    edge_id(c%nd(v)%out%a(i)) = m
                end do
            end do
            n = nodes%n
            call ns_allocate(s, n, m)
            at = 1
            do j = 1, n
                v = nodes%a(j)
                s%out_start(j) = at
                do i = 1, c%nd(v)%out%n
                    e = c%nd(v)%out%a(i)
                    s%tail(edge_id(e)) = j
                    s%head(edge_id(e)) = local(c%ed(e)%head)
                    s%minlen(edge_id(e)) = c%ed(e)%minlen
                    s%weight(edge_id(e)) = c%ed(e)%weight
                    s%out_list(at) = edge_id(e)
                    at = at + 1
                end do
            end do
            s%out_start(n + 1) = at
            at = 1
            do j = 1, n
                v = nodes%a(j)
                s%in_start(j) = at
                do i = 1, c%nd(v)%in%n
                    s%in_list(at) = edge_id(c%nd(v)%in%a(i))
                    at = at + 1
                end do
            end do
            s%in_start(n + 1) = at
            call ns_solve(s, 1, huge(1), stat)
            if (stat /= 0) return
            do j = 1, n
                c%nd(nodes%a(j))%rank = s%rank(j)
            end do
        end associate
    end do
end subroutine rank1

end module fdot_rank
