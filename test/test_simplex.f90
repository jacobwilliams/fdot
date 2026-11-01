module test_simplex
!! Tests of network simplex.
use testdrive, only: new_unittest, unittest_type, error_type, check
use fdot_simplex, only: ns_t, ns_allocate, ns_finish_lists, ns_solve
implicit none
private
public :: collect_simplex

contains

!> Collect the tests.
subroutine collect_simplex(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)   !! the tests

    testsuite = [new_unittest("diamond", test_diamond), new_unittest("disconnected", test_disconnected), &
                 new_unittest("long chain", test_chain)]
end subroutine collect_simplex

!> Fill in an edge.
subroutine set_edge(s, e, tail, head, minlen, weight)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: e        !! the edge
    integer, intent(in) :: tail     !! where it starts
    integer, intent(in) :: head     !! where it ends
    integer, intent(in) :: minlen   !! its smallest length
    integer, intent(in) :: weight   !! its weight

    s%tail(e) = tail
    s%head(e) = head
    s%minlen(e) = minlen
    s%weight(e) = weight
end subroutine set_edge

!> The cost is the smallest there is, with the short side stretched.
subroutine test_diamond(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(ns_t) :: s
    integer :: stat

    call ns_allocate(s, 4, 4)
    call set_edge(s, 1, 1, 2, 1, 1)
    call set_edge(s, 2, 2, 4, 3, 1)
    call set_edge(s, 3, 1, 3, 1, 5)
    call set_edge(s, 4, 3, 4, 1, 1)
    call ns_finish_lists(s)
    call ns_solve(s, 0, huge(1), stat)
    call check(error, stat, 0)
    if (allocated(error)) return
    call check(error, s%rank(1), 0)
    if (allocated(error)) return
    call check(error, s%rank(4) >= 4)
    if (allocated(error)) return
    ! the heavy edge stays short
    call check(error, s%rank(3) - s%rank(1), 1)
end subroutine test_diamond

!> A graph in two parts is reported.
subroutine test_disconnected(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(ns_t) :: s
    integer :: stat

    call ns_allocate(s, 4, 2)
    call set_edge(s, 1, 1, 2, 1, 1)
    call set_edge(s, 2, 3, 4, 1, 1)
    call ns_finish_lists(s)
    call ns_solve(s, 0, huge(1), stat)
    call check(error, stat, 1)
end subroutine test_disconnected

!> A long chain does not need a deep stack.
subroutine test_chain(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    integer, parameter :: N = 100000
    type(ns_t) :: s
    integer :: stat, k

    call ns_allocate(s, N, N - 1)
    do k = 1, N - 1
        call set_edge(s, k, k, k + 1, 1, 1)
    end do
    call ns_finish_lists(s)
    call ns_solve(s, 0, huge(1), stat)
    call check(error, stat, 0)
    if (allocated(error)) return
    call check(error, s%rank(N), N - 1)
end subroutine test_chain

end module test_simplex
