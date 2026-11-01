module test_text
!! Tests of text and node sizes.
use testdrive, only: new_unittest, unittest_type, error_type, check
use fdot_structure, only: wp
use fdot_text, only: text_width, text_height
use fdot_draw, only: node_size
implicit none
private
public :: collect_text

contains

!> Collect the tests.
subroutine collect_text(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)   !! the tests

    testsuite = [new_unittest("width", test_width), new_unittest("sizes", test_sizes)]
end subroutine collect_text

!> Widths come from the table.
subroutine test_width(error)
    type(error_type), allocatable, intent(out) :: error !! the result

    call check(error, abs(text_width('A', 2048.0_wp) - 1366) < 1.0e-9_wp)
    if (allocated(error)) return
    call check(error, abs(text_width('iii', 2048.0_wp) - 3*455) < 1.0e-9_wp)
    if (allocated(error)) return
    call check(error, text_width('', 13.0_wp) == 0)
    if (allocated(error)) return
    ! bytes that are not ASCII count as spaces
    call check(error, abs(text_width(achar(200), 2048.0_wp) - 569) < 1.0e-9_wp)
    if (allocated(error)) return
    call check(error, abs(text_height(10.0_wp) - 12) < 1.0e-9_wp)
end subroutine test_width

!> A box has the smallest width and height unless its text needs more.
subroutine test_sizes(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    real(wp) :: w, h

    call node_size(10.0_wp, 15.6_wp, 7.2_wp, 2.88_wp, 21.6_wp, w, h)
    call check(error, w == 54 .and. h == 21.6_wp)
    if (allocated(error)) return
    call node_size(100.0_wp, 30.0_wp, 7.2_wp, 2.88_wp, 21.6_wp, w, h)
    call check(error, abs(w - 114.4_wp) < 1.0e-9_wp .and. abs(h - 35.76_wp) < 1.0e-9_wp)
end subroutine test_sizes

end module test_text
