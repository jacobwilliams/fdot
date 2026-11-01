program usage
!! The example of the README: draw a small graph and print it.
use fdot, only: graph_t, RANK_RL
implicit none

type(graph_t) :: g
character(len=:), allocatable :: svg
integer :: a, b, stat

call g%initialize(direction=RANK_RL, concentrate=.true.)
call g%add_node('json_module', a, fill='#3657c9', link='../module/json_module.html')
call g%add_node('json_kinds', b, fill='#646d7e')
call g%add_edge(a, b, colour='#8a93a5')
call g%layout(stat)
if (stat /= 0) error stop 'layout failed'
call g%write_svg(svg)
write (*, '(a)') svg
end program usage
