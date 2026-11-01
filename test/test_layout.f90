module test_layout
!! Tests of the stages of layout, through the public interface.
use testdrive, only: new_unittest, unittest_type, error_type, check
use fdot, only: graph_t, point_t, wp, RANK_LR, RANK_RL, FDOT_OK, FDOT_INVALID
implicit none
private
public :: collect_layout

contains

!> Collect the tests.
subroutine collect_layout(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)   !! the tests

    testsuite = [new_unittest("empty", test_empty), new_unittest("one node", test_one), &
                 new_unittest("chain", test_chain), new_unittest("direction", test_direction), &
                 new_unittest("cycle", test_cycle), new_unittest("labels double ranks", test_labels), &
                 new_unittest("self-loop", test_self_loop), new_unittest("components", test_components), &
                 new_unittest("invalid edge", test_invalid), new_unittest("opposite edges", test_opposite), &
                 new_unittest("concentrate", test_concentrate), new_unittest("parallel edges", test_parallel), &
                 new_unittest("long edge", test_long_edge), new_unittest("no overlap", test_no_overlap)]
end subroutine collect_layout

!> Add `n` nodes.
subroutine add_nodes(g, n)
    type(graph_t), intent(inout) :: g   !! the graph
    integer, intent(in) :: n            !! how many
    integer :: k, id

    do k = 1, n
        call g%add_node('node', id)
    end do
end subroutine add_nodes

!> The rank of a node.
integer function rank_of(g, id)
    type(graph_t), intent(in) :: g  !! the graph
    integer, intent(in) :: id       !! the node
    real(wp) :: x, y, w, h

    call g%node_position(id, x, y, w, h, rank=rank_of)
end function rank_of

!> A graph with no nodes has an empty drawing.
subroutine test_empty(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat
    character(len=:), allocatable :: svg

    call g%initialize()
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call g%write_svg(svg)
    call check(error, len(svg), 0)
end subroutine test_empty

!> A node alone has the smallest size.
subroutine test_one(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat
    real(wp) :: w, h

    call g%initialize()
    call add_nodes(g, 1)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call g%size_of_drawing(w, h)
    call check(error, abs(w - 54) < 1.0e-9_wp .and. abs(h - 21.6_wp) < 1.0e-9_wp)
end subroutine test_one

!> Ranks follow the edges, and the drawing goes from left to right.
subroutine test_chain(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat, k
    real(wp) :: x1, x2, y, w, h

    call g%initialize()
    call add_nodes(g, 4)
    do k = 1, 3
        call g%add_edge(k, k + 1)
    end do
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    do k = 1, 4
        call check(error, rank_of(g, k), k - 1)
        if (allocated(error)) return
    end do
    call g%node_position(1, x1, y, w, h)
    call g%node_position(4, x2, y, w, h)
    call check(error, x2 > x1)
end subroutine test_chain

!> Right to left mirrors the drawing.
subroutine test_direction(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat
    real(wp) :: x1, x2, y, w, h

    call g%initialize(direction=RANK_RL)
    call add_nodes(g, 2)
    call g%add_edge(1, 2)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call g%node_position(1, x1, y, w, h)
    call g%node_position(2, x2, y, w, h)
    call check(error, x2 < x1)
end subroutine test_direction

!> A cycle is broken, and every edge still has a curve from tail to head.
subroutine test_cycle(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    type(point_t), allocatable :: points(:)
    integer :: stat, k
    real(wp) :: x, y, w, h

    call g%initialize()
    call add_nodes(g, 3)
    call g%add_edge(1, 2)
    call g%add_edge(2, 3)
    call g%add_edge(3, 1)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, rank_of(g, 3), 2)
    if (allocated(error)) return
    call g%edge_curve(3, 1, points)
    call g%node_position(3, x, y, w, h)
    ! the edge back starts at node 3, at the far end
    call check(error, abs(points(1)%x - x) <= w/2 + 1)
end subroutine test_cycle

!> An edge label gets a rank of its own.
subroutine test_labels(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    logical :: has_label
    integer :: stat
    real(wp) :: x, y

    call g%initialize()
    call add_nodes(g, 2)
    call g%add_edge(1, 2, label='has')
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, rank_of(g, 2), 2)
    if (allocated(error)) return
    call g%edge_label_position(1, x, y, has_label)
    call check(error, has_label .and. x > 0)
end subroutine test_labels

!> A self-loop has a curve of two pieces, and does not disturb the ranks.
subroutine test_self_loop(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    type(point_t), allocatable :: points(:)
    integer :: stat

    call g%initialize()
    call add_nodes(g, 2)
    call g%add_edge(1, 1)
    call g%add_edge(1, 2)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, g%edge_curves(1), 1)
    if (allocated(error)) return
    call g%edge_curve(1, 1, points)
    call check(error, size(points), 7)
    if (allocated(error)) return
    call check(error, rank_of(g, 2), 1)
end subroutine test_self_loop

!> Components are drawn together, each from rank 0.
subroutine test_components(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat

    call g%initialize()
    call add_nodes(g, 5)
    call g%add_edge(1, 2)
    call g%add_edge(3, 4)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, rank_of(g, 3), 0)
    if (allocated(error)) return
    call check(error, rank_of(g, 5), 0)
end subroutine test_components

!> An edge to a node that does not exist is reported.
subroutine test_invalid(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat

    call g%initialize()
    call add_nodes(g, 1)
    call g%add_edge(1, 7)
    call g%layout(stat)
    call check(error, stat, FDOT_INVALID)
end subroutine test_invalid

!> With concentration, opposite edges are one curve with arrowheads at both ends.
subroutine test_opposite(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    type(point_t), allocatable :: points(:)
    logical :: has_start, has_end
    integer :: stat

    call g%initialize(concentrate=.true.)
    call add_nodes(g, 2)
    call g%add_edge(1, 2)
    call g%add_edge(2, 1)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, g%edge_curves(1) + g%edge_curves(2), 1)
    if (allocated(error)) return
    call g%edge_curve(1, 1, points, has_start=has_start, has_end=has_end)
    call check(error, has_start .and. has_end)
end subroutine test_opposite

!> With concentration, edges that run side by side share their virtual nodes.
subroutine test_concentrate(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat, k, plain, merged

    do k = 1, 2
        call g%initialize(concentrate=k == 2)
        call add_nodes(g, 8)
        ! two sources feed long edges to one sink, past a chain that sets the ranks
        call g%add_edge(1, 3)
        call g%add_edge(3, 4)
        call g%add_edge(4, 5)
        call g%add_edge(5, 6)
        call g%add_edge(1, 6)
        call g%add_edge(2, 6)
        call g%add_edge(1, 7)
        call g%add_edge(2, 7)
        call g%add_edge(5, 8)
        call g%layout(stat)
        call check(error, stat, FDOT_OK)
        if (allocated(error)) return
        if (k == 1) then
            plain = total_curves(g)
        else
            merged = total_curves(g)
        end if
    end do
    call check(error, plain, g%edge_count())
    if (allocated(error)) return
    call check(error, merged >= plain)
end subroutine test_concentrate

!> The number of curves of a graph.
integer function total_curves(g)
    type(graph_t), intent(in) :: g  !! the graph
    integer :: e

    total_curves = 0
    do e = 1, g%edge_count()
        total_curves = total_curves + g%edge_curves(e)
    end do
end function total_curves

!> Without concentration, parallel edges each have a curve, apart from each other.
subroutine test_parallel(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    type(point_t), allocatable :: a(:), b(:)
    integer :: stat

    call g%initialize()
    call add_nodes(g, 2)
    call g%add_edge(1, 2)
    call g%add_edge(1, 2)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, g%edge_curves(1), 1)
    if (allocated(error)) return
    call check(error, g%edge_curves(2), 1)
    if (allocated(error)) return
    call g%edge_curve(1, 1, a)
    call g%edge_curve(2, 1, b)
    call check(error, abs(a(2)%y - b(2)%y) > 1)
end subroutine test_parallel

!> An edge over many ranks gets one curve, of several pieces or straight.
subroutine test_long_edge(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    integer :: stat, k

    call g%initialize()
    call add_nodes(g, 12)
    do k = 1, 11
        call g%add_edge(k, k + 1)
    end do
    call g%add_edge(1, 12)
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    call check(error, g%edge_curves(12), 1)
end subroutine test_long_edge

!> Nodes of one rank do not overlap.
subroutine test_no_overlap(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    real(wp) :: x1, y1, x2, y2, w, h
    integer :: stat, k, j

    call g%initialize()
    call add_nodes(g, 21)
    do k = 2, 21
        call g%add_edge(1, k)
    end do
    call g%layout(stat)
    call check(error, stat, FDOT_OK)
    if (allocated(error)) return
    do k = 2, 21
        call g%node_position(k, x1, y1, w, h)
        do j = k + 1, 21
            call g%node_position(j, x2, y2, w, h)
            call check(error, abs(y1 - y2) >= h)
            if (allocated(error)) return
        end do
    end do
end subroutine test_no_overlap

end module test_layout
