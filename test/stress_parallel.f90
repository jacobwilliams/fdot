program stress_parallel
!! Lay out many graphs at once, and check that the result does not depend on the number of threads.
!! This needs OpenMP, which the library itself does not use: `pixi run test-threads`.
use iso_fortran_env, only: int64, error_unit
use fdot, only: graph_t
implicit none

integer, parameter :: GRAPHS = 400
integer(int64) :: reference(GRAPHS), sums(GRAPHS)
integer, parameter :: THREADS(3) = [1, 4, 10]
integer :: t, k

reference = 0
do t = 1, size(THREADS)
    sums = 0
    !$omp parallel do schedule(dynamic) num_threads(THREADS(t))
    do k = 1, GRAPHS
        sums(k) = lay_out(k)
    end do
    !$omp end parallel do
    if (t == 1) then
        reference = sums
    else if (any(sums /= reference)) then
        write (error_unit, '(a,i0,a,i0)') 'the output differs with ', THREADS(t), ' threads, for graph ', findloc(sums /= reference, .true., 1)
        error stop
    end if
end do
if (any(reference == 0)) error stop 'a layout failed'
write (*, '(a,i0,a)') 'the same output for ', GRAPHS, ' graphs with 1, 4 and 10 threads'

contains

!> Lay out random graph `seed`, and return a checksum of the drawing.
integer(int64) function lay_out(seed) result(sum)
    integer, intent(in) :: seed !! the graph
    type(graph_t) :: g
    character(len=:), allocatable :: svg
    integer(int64) :: state
    integer :: n, m, k, id, stat

    state = 12345_int64 + seed
    n = 2 + next(state, 80)
    m = n + next(state, 2*n)
    call g%initialize(concentrate=mod(seed, 3) /= 0)
    do k = 1, n
        call g%add_node(repeat('x', 1 + next(state, 12)), id, fill='#646d7e', link=merge('a.html', '      ', mod(k, 2) == 0))
    end do
    do k = 1, m
        call g%add_edge(1 + next(state, n), 1 + next(state, n), label=merge('l', ' ', next(state, 6) == 0))
    end do
    call g%layout(stat)
    sum = 0
    if (stat /= 0) return
    call g%write_svg(svg)
    do k = 1, len(svg)
        sum = mod(31*sum + iachar(svg(k:k)), 2147483647_int64)
    end do
    sum = sum + 1
end function lay_out

!> The next random number below `limit`.
integer function next(state, limit)
    integer(int64), intent(inout) :: state  !! the state
    integer, intent(in) :: limit            !! the limit

    state = mod(1103515245_int64*state + 12345, 2147483648_int64)
    next = int(mod(state/65536, int(limit, int64)))
end function next

end program stress_parallel
