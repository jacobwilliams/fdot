module fdot_position
!! Setting coordinates: y a rank at a time, and x from a second network simplex on an auxiliary graph.
!! This is a port of `lib/dotgen/position.c` of Graphviz, which is under the Eclipse Public License
!! 2.0, without clusters and flat edges.
use fdot_structure
use fdot_concentrate, only: dot_concentrate
use fdot_simplex, only: ns_t, ns_allocate, ns_finish_lists, ns_solve
implicit none
private

public :: dot_position

integer, parameter :: SELF_EDGE_SIZE = 18   !! the width of a self-loop
integer, parameter :: LABEL_SEP = 5         !! points between nodes on the extra ranks of labels

type :: aux_t
    !! The auxiliary graph. Nodes above 0 are nodes of the layout, and below 0 are slack nodes.
    type(ilist_t) :: tail           !! where each edge starts
    type(ilist_t) :: head           !! where each edge ends
    type(ilist_t) :: len            !! the smallest length of each edge
    type(ilist_t) :: weight         !! the weight of each edge
    type(ilist_t) :: slack_rank     !! the first position of each slack node
end type aux_t

contains

!> Position the nodes. All components start at rank 0, and neighbours in a rank are linked, so the
!> auxiliary graph is connected, and the `connectGraph` of the original is not needed.
subroutine dot_position(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    integer, allocatable :: start(:)
    type(aux_t) :: aux

    stat = 0
    call set_ycoords(c)
    if (c%concentrate) call dot_concentrate(c)
    call make_lr_constraints(c, aux, start)
    call make_edge_pairs(c, aux, start)
    call solve(c, aux, start, stat)
    if (stat /= 0) return
    call compute_bb(c)
end subroutine dot_position

!> Add a constraint.
subroutine make_aux_edge(aux, u, v, len, weight)
    type(aux_t), intent(inout) :: aux   !! the graph
    integer, intent(in) :: u            !! where it starts
    integer, intent(in) :: v            !! where it ends
    real(wp), intent(in) :: len         !! its smallest length
    integer, intent(in) :: weight       !! its weight

    call push(aux%tail, u)
    call push(aux%head, v)
    call push(aux%len, round(len))
    call push(aux%weight, weight)
end subroutine make_aux_edge

!> The space a self-loop needs to the right of its node.
function self_right_space(c, e) result(sw)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e        !! the edge
    real(wp) :: sw

    sw = SELF_EDGE_SIZE
    if (c%labelled(e)) then
        if (c%flip) then
            sw = sw + c%label_h(e)
        else
            sw = sw + c%label_w(e)
        end if
    end if
end function self_right_space

!> Add the constraints that keep the nodes of each rank apart and in order, and the first positions.
subroutine make_lr_constraints(c, aux, start)
    type(core_t), intent(inout) :: c            !! the state
    type(aux_t), intent(inout) :: aux           !! the graph
    integer, allocatable, intent(out) :: start(:)   !! the first position of each node
    integer :: i, j, k, u, v, e, sep(0:1), nodesep_here
    real(wp) :: width, last, sw

    allocate (start(c%nnodes), source=0)
    sep = NODESEP
    if (c%has_labels) sep(1) = LABEL_SEP
    do i = c%minrank, c%maxrank
        if (c%rk(i)%v(0) == 0) cycle
        last = 0
        start(c%rk(i)%v(0)) = 0
        nodesep_here = sep(iand(i, 1))
        do j = 0, c%rk(i)%n - 1
            u = c%rk(i)%v(j)
            c%nd(u)%saved_rw = c%nd(u)%rw
            if (c%nd(u)%other%n > 0) then
                sw = 0
                do k = 1, c%nd(u)%other%n
                    e = c%nd(u)%other%a(k)
                    if (c%ed(e)%tail == c%ed(e)%head) sw = sw + self_right_space(c, e)
                end do
                c%nd(u)%rw = c%nd(u)%rw + sw
            end if
            v = c%rk(i)%v(j + 1)
            if (v /= 0) then
                width = c%nd(u)%rw + c%nd(v)%lw + nodesep_here
                call make_aux_edge(aux, u, v, width, 0)
                start(v) = int(last + width)
                last = start(v)
            end if
        end do
    end do
end subroutine make_lr_constraints

!> Add a slack node for each edge, whose cost is the distance between the ends of the edge.
subroutine make_edge_pairs(c, aux, start)
    type(core_t), intent(in) :: c       !! the state
    type(aux_t), intent(inout) :: aux   !! the graph
    integer, intent(in) :: start(:)     !! the first position of each node
    integer :: i, k, n, e, sn

    do k = 1, size(c%nlist)
        n = c%nlist(k)
        do i = 1, c%nd(n)%out%n
            e = c%nd(n)%out%a(i)
            call push(aux%slack_rank, min(start(c%ed(e)%tail) - 1, start(c%ed(e)%head) - 1))
            sn = -aux%slack_rank%n
            call make_aux_edge(aux, sn, c%ed(e)%tail, 1.0_wp, c%ed(e)%weight)
            call make_aux_edge(aux, sn, c%ed(e)%head, 1.0_wp, c%ed(e)%weight)
        end do
    end do
end subroutine make_edge_pairs

!> Solve for x. New slack nodes go to the front of the list of nodes, as they do in the original.
subroutine solve(c, aux, start, stat)
    type(core_t), intent(inout) :: c    !! the state
    type(aux_t), intent(in) :: aux      !! the graph
    integer, intent(in) :: start(:)     !! the first position of each node
    integer, intent(out) :: stat        !! the result
    type(ns_t) :: s
    integer, allocatable :: local(:)
    integer :: ns, nr, i, e

    ns = aux%slack_rank%n
    nr = size(c%nlist)
    allocate (local(c%nnodes), source=0)
    do i = 1, nr
        local(c%nlist(i)) = ns + i
    end do
    call ns_allocate(s, ns + nr, aux%tail%n)
    do i = 1, ns
        s%rank(ns - i + 1) = aux%slack_rank%a(i)
        s%normal(ns - i + 1) = .false.
    end do
    do i = 1, nr
        s%rank(ns + i) = start(c%nlist(i))
    end do
    do e = 1, aux%tail%n
        s%tail(e) = node_of(aux%tail%a(e))
        s%head(e) = node_of(aux%head%a(e))
        s%minlen(e) = aux%len%a(e)
        s%weight(e) = aux%weight%a(e)
    end do
    call ns_finish_lists(s)
    call ns_solve(s, 2, huge(1), stat)
    if (stat /= 0) return
    do i = 1, nr
        c%nd(c%nlist(i))%x = s%rank(ns + i)
    end do

contains

    !> The node of the solver for a node of the auxiliary graph.
    integer function node_of(id)
        integer, intent(in) :: id   !! the node of the auxiliary graph

        if (id > 0) then
            node_of = local(id)
        else
            node_of = ns + id + 1
        end if
    end function node_of
end subroutine solve

!> Set the y coordinates, a rank at a time.
subroutine set_ycoords(c)
    type(core_t), intent(inout) :: c    !! the state
    integer :: r, i, j, n, e
    real(wp) :: ht2, delta, y

    c%ht1 = 0
    c%ht2 = 0
    do r = c%minrank, c%maxrank
        do i = 0, c%rk(r)%n - 1
            n = c%rk(r)%v(i)
            ht2 = c%nd(n)%ht/2
            do j = 1, c%nd(n)%other%n
                e = c%nd(n)%other%a(j)
                if (c%ed(e)%tail == c%ed(e)%head .and. c%labelled(e)) then
                    ht2 = max(ht2, c%label_h(e)/2)
                end if
            end do
            if (c%rk(r)%pht2 < ht2) then
                c%rk(r)%pht2 = ht2
                c%rk(r)%ht2 = ht2
            end if
            if (c%rk(r)%pht1 < ht2) then
                c%rk(r)%pht1 = ht2
                c%rk(r)%ht1 = ht2
            end if
            if (r == c%minrank) c%ht2 = max(c%ht2, ht2)
            if (r == c%maxrank) c%ht1 = max(c%ht1, ht2)
        end do
    end do

    y = c%rk(c%maxrank)%ht1
    call set_rank_y(c, c%maxrank, y)
    do r = c%maxrank - 1, c%minrank, -1
        delta = c%rk(r + 1)%pht2 + c%rk(r)%pht1 + c%ranksep
        y = y + delta
        call set_rank_y(c, r, y)
    end do
end subroutine set_ycoords

!> Set the y coordinate of a rank.
subroutine set_rank_y(c, r, y)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: r            !! the rank
    real(wp), intent(in) :: y           !! the coordinate
    integer :: i

    do i = 0, c%rk(r)%n - 1
        c%nd(c%rk(r)%v(i))%y = y
    end do
end subroutine set_rank_y

!> Find the extent of the nodes.
subroutine compute_bb(c)
    type(core_t), intent(inout) :: c    !! the state
    integer :: r, i, v, first, last

    c%bb%llx = huge(1)
    c%bb%urx = -huge(1)
    do r = c%minrank, c%maxrank
        first = -1
        last = -1
        do i = 0, c%rk(r)%n - 1
            if (c%nd(c%rk(r)%v(i))%kind == NORMAL) then
                if (first < 0) first = i
                last = i
            end if
        end do
        if (first < 0) cycle
        v = c%rk(r)%v(first)
        c%bb%llx = min(c%bb%llx, c%nd(v)%x - c%nd(v)%lw)
        v = c%rk(r)%v(last)
        c%bb%urx = max(c%bb%urx, c%nd(v)%x + c%nd(v)%rw)
    end do
    c%bb%lly = c%nd(c%rk(c%maxrank)%v(0))%y - c%ht1
    do r = c%minrank, c%maxrank
        if (c%rk(r)%n == 0) cycle
        c%bb%ury = c%nd(c%rk(r)%v(0))%y + c%ht2
        exit
    end do
end subroutine compute_bb

end module fdot_position
