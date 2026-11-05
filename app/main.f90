program fdot_layout
!! Lay out graphs given on standard input, for the checks in `tools/`. Reading and writing files is
!! left out of the library, so it is here. Each line has fields separated by tabs:
!!
!!     G  <LR or RL>  <1 to concentrate, or 0>
!!     N  <name>  <fill>  <link>  <outline>
!!     E  <from>  <to>  <label>  <colour>
!!     L
!!
!! `L` lays out the graph so far. With the argument `svg`, the drawing is written, followed by a line
!! with `%%`. Otherwise the layout is written in lines of numbers, with `end` after each graph.
use iso_fortran_env, only: output_unit, error_unit, iostat_end
use fdot, only: graph_t, point_t, wp, RANK_LR, RANK_RL, FDOT_OK
implicit none

character(len=:), allocatable :: line, fields(:), svg
character(len=16) :: mode
type(graph_t) :: g
integer :: status, id

mode = ''
if (command_argument_count() > 0) call get_command_argument(1, mode)
do
    call read_line(line, status)
    if (status /= 0) exit
    if (len(line) == 0) cycle
    call split(line, fields)
    select case (trim(fields(1)))
    case ('G')
        call g%initialize(direction=merge(RANK_RL, RANK_LR, trim(fields(2)) == 'RL'), concentrate=trim(fields(3)) == '1')
    case ('N')
        call g%add_node(trim(fields(2)), id, fill=trim(fields(3)), link=trim(fields(4)), outline=trim(fields(5)))
    case ('E')
        call g%add_edge(to_int(fields(2)), to_int(fields(3)), label=trim(fields(4)), colour=trim(fields(5)))
    case ('L')
        call g%layout(status)
        if (status /= FDOT_OK) then
            write (error_unit, '(a,i0)') 'layout failed: ', status
            write (output_unit, '(a,i0)') 'failed ', status
        else if (trim(mode) == 'svg') then
            call g%write_svg(svg)
            write (output_unit, '(a)') svg
        else
            call report(g)
        end if
        write (output_unit, '(a)') merge('%%% ', 'end ', trim(mode) == 'svg')
    end select
end do

contains

!> Write the layout as numbers.
subroutine report(g)
    type(graph_t), intent(in) :: g  !! the graph
    type(point_t), allocatable :: points(:)
    real(wp) :: x, y, w, h
    integer :: k, j, rank, order
    logical :: has_label

    call g%size_of_drawing(w, h)
    write (output_unit, '(a,2(1x,f0.3))') 'graph', w, h
    do k = 1, g%node_count()
        call g%node_position(k, x, y, w, h, rank, order)
        write (output_unit, '(a,3(1x,i0),4(1x,f0.3))') 'node', k, rank, order, x, y, w, h
    end do
    do k = 1, g%edge_count()
        do j = 1, g%edge_curves(k)
            call g%edge_curve(k, j, points)
            write (output_unit, '(a,2(1x,i0),*(1x,f0.3))') 'edge', k, size(points), &
                (points(rank)%x, points(rank)%y, rank=1, size(points))
        end do
        call g%edge_label_position(k, x, y, has_label)
        if (has_label) write (output_unit, '(a,1x,i0,2(1x,f0.3))') 'label', k, x, y
    end do
end subroutine report

!> Read a line of any length from standard input.
subroutine read_line(line, status)
    character(len=:), allocatable, intent(out) :: line  !! the line
    integer, intent(out) :: status                      !! not 0 at the end of the input
    character(len=1024) :: chunk
    integer :: size_read, ios

    line = ''
    status = 0
    do
        read (*, '(a)', advance='no', size=size_read, iostat=ios) chunk
        if (ios == iostat_end) then
            status = 1
            return
        end if
        line = line//chunk(:size_read)
        if (ios /= 0) return
    end do
end subroutine read_line

!> Split a line at tabs.
subroutine split(line, fields)
    character(len=*), intent(in) :: line                    !! the line
    character(len=:), allocatable, intent(out) :: fields(:) !! the fields
    integer, parameter :: MOST = 6         !! the most fields a line has
    integer :: k, from, at

    allocate (character(len=max(1, len(line))) :: fields(MOST))
    fields = ''
    from = 1
    do k = 1, MOST
        at = index(line(from:), achar(9))
        if (at == 0) then
            fields(k) = line(from:)
            return
        end if
        fields(k) = line(from:from + at - 2)
        from = from + at
    end do
end subroutine split

!> The value of a number.
integer function to_int(text)
    character(len=*), intent(in) :: text    !! the number

    read (text, *) to_int
end function to_int

end program fdot_layout
