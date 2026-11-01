module test_output
!! Tests of what the library writes.
use iso_fortran_env, only: int64
use testdrive, only: new_unittest, unittest_type, error_type, check
use fdot, only: graph_t, wp, RANK_RL, node_style_t, edge_style_t
use fdot_svg, only: buffer_t, add_tenths, add_int, add_xml, take
implicit none
private
public :: collect_output

contains

!> Collect the tests.
subroutine collect_output(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)   !! the tests

    testsuite = [new_unittest("numbers", test_numbers), new_unittest("svg", test_svg), &
                 new_unittest("dot", test_dot), new_unittest("repeatable", test_repeatable), &
                 new_unittest("styles", test_styles)]
end subroutine collect_output

!> The text of a number.
subroutine tenths(value, expected, error)
    real(wp), intent(in) :: value                       !! the value
    character(len=*), intent(in) :: expected            !! what it should give
    type(error_type), allocatable, intent(inout) :: error   !! the result
    type(buffer_t) :: b
    character(len=:), allocatable :: text

    call add_tenths(b, value)
    call take(b, text)
    call check(error, text, expected)
end subroutine tenths

!> Numbers are written without internal I/O.
subroutine test_numbers(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(buffer_t) :: b
    character(len=:), allocatable :: text

    call tenths(0.0_wp, '0', error)
    if (allocated(error)) return
    call tenths(-0.04_wp, '0', error)
    if (allocated(error)) return
    call tenths(12.34_wp, '12.3', error)
    if (allocated(error)) return
    call tenths(-1.26_wp, '-1.3', error)
    if (allocated(error)) return
    call tenths(99.96_wp, '100', error)
    if (allocated(error)) return
    call add_int(b, -120_int64)
    call add_xml(b, 'a<b&"c">')
    call take(b, text)
    call check(error, text, '-120a&lt;b&amp;&quot;c&quot;&gt;')
end subroutine test_numbers

!> The drawing has the compact form that fdoc keeps.
subroutine test_svg(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    character(len=:), allocatable :: svg
    integer :: a, b, stat

    call g%initialize(direction=RANK_RL, concentrate=.true.)
    call g%add_node('json_module', a, fill='#3657c9', link='../module/json_module.html')
    call g%add_node('json<kinds', b, fill='#646d7e', outline='#f0a868')
    call g%add_edge(a, b, colour='#112233', label='uses')
    call g%layout(stat)
    call check(error, stat, 0)
    if (allocated(error)) return
    call g%write_svg(svg)
    call check(error, svg(1:12), '<svg width="')
    if (allocated(error)) return
    call contains_text(svg, '<a href="../module/json_module.html" fill="white" stroke="none"><rect fill="#3657c9" rx="7.2"', error)
    if (allocated(error)) return
    call contains_text(svg, 'stroke="#f0a868" stroke-width="2" rx="7.2"', error)
    if (allocated(error)) return
    call contains_text(svg, '>json&lt;kinds</text>', error)
    if (allocated(error)) return
    call contains_text(svg, '<path stroke="#112233" stroke-dasharray="5,2" d="M', error)
    if (allocated(error)) return
    call contains_text(svg, '<polygon fill="#112233"', error)
    if (allocated(error)) return
    call contains_text(svg, 'font-size="11" fill="#8a93a5">uses</text>', error)
    if (allocated(error)) return
    call check(error, svg(len(svg) - 6:len(svg) - 1), '</svg>')
end subroutine test_svg

!> Fail unless `text` holds `part`.
subroutine contains_text(text, part, error)
    character(len=*), intent(in) :: text                !! the text
    character(len=*), intent(in) :: part                !! what it must hold
    type(error_type), allocatable, intent(inout) :: error   !! the result

    call check(error, index(text, part) > 0, 'missing: '//part)
end subroutine contains_text

!> The graph can be written for `dot`.
subroutine test_dot(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    character(len=:), allocatable :: text
    integer :: a, b

    call g%initialize(direction=RANK_RL, concentrate=.true.)
    call g%add_node('say "hi"', a)
    call g%add_node('b', b, link='x.html')
    call g%add_edge(a, b, label='l')
    call g%write_dot(text)
    call contains_text(text, 'rankdir=RL', error)
    if (allocated(error)) return
    call contains_text(text, 'concentrate=true', error)
    if (allocated(error)) return
    call contains_text(text, 'n1 [label="say \"hi\""', error)
    if (allocated(error)) return
    call contains_text(text, 'URL="x.html"', error)
    if (allocated(error)) return
    call contains_text(text, 'n1 -> n2 [', error)
    if (allocated(error)) return
    call contains_text(text, 'label="l"', error)
end subroutine test_dot

!> The same input gives the same output.
subroutine test_repeatable(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    character(len=:), allocatable :: first, second
    integer :: k, id, stat

    call g%initialize(concentrate=.true.)
    do k = 1, 30
        call g%add_node('n', id)
    end do
    do k = 1, 29
        call g%add_edge(k, mod(k*7, 30) + 1)
        call g%add_edge(k, k + 1)
    end do
    call g%layout(stat)
    call check(error, stat, 0)
    if (allocated(error)) return
    call g%write_svg(first)
    call g%layout(stat)
    call g%write_svg(second)
    call check(error, first == second)
end subroutine test_repeatable

!> Styles change what is drawn.
subroutine test_styles(error)
    type(error_type), allocatable, intent(out) :: error !! the result
    type(graph_t) :: g
    type(edge_style_t) :: edges
    real(wp) :: x, y, w, h
    integer :: a, b, stat

    edges%colour = '#000000'
    call g%initialize(node_style=node_style_t(height=40.0_wp), edge_style=edges)
    call g%add_node('a', a)
    call g%add_node('b', b)
    call g%add_edge(a, b)
    call g%layout(stat)
    call check(error, stat, 0)
    if (allocated(error)) return
    call g%node_position(a, x, y, w, h)
    call check(error, abs(h - 40) < 1.0e-9_wp)
end subroutine test_styles

end module test_output
