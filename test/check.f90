program check
!! Unit tests.
use iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type
use test_text, only: collect_text
use test_simplex, only: collect_simplex
use test_layout, only: collect_layout
use test_output, only: collect_output
implicit none

type(testsuite_type), allocatable :: testsuites(:)
integer :: stat, is

stat = 0
testsuites = [ &
             new_testsuite("text", collect_text), &
             new_testsuite("simplex", collect_simplex), &
             new_testsuite("layout", collect_layout), &
             new_testsuite("output", collect_output) &
             ]
do is = 1, size(testsuites)
    write (error_unit, '(a)') '# Testing: '//testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
end do
if (stat > 0) then
    write (error_unit, '(i0, a)') stat, ' test(s) failed!'
    error stop
end if
end program check
