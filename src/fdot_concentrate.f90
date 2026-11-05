module fdot_concentrate
!! Merging the virtual nodes of edges that run side by side, so that they are drawn as one. This is a
!! port of `lib/dotgen/conc.c` of Graphviz, which is under the Eclipse Public License 2.0.
use fdot_structure
implicit none
private

public :: dot_concentrate

integer, parameter :: UP = 0       !! merging nodes whose edges go to the same node
integer, parameter :: DOWN = 1     !! merging nodes whose edges come from the same node

contains

!> Whether `e` and `f` go in the same direction.
pure logical function samedir(c, e, f)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e        !! one edge
    integer, intent(in) :: f        !! the other edge
    integer :: e0, f0

    samedir = .false.
    e0 = normal_edge(c, e)
    if (e0 == 0) return
    f0 = normal_edge(c, f)
    if (f0 == 0) return
    if (c%ed(e0)%conc_opp .or. c%ed(f0)%conc_opp) return
    samedir = (c%nd(c%ed(f0)%tail)%rank - c%nd(c%ed(f0)%head)%rank) &
              *(c%nd(c%ed(e0)%tail)%rank - c%nd(c%ed(e0)%head)%rank) > 0
end function samedir

!> The edge of the drawing that `e` comes from, or 0.
pure integer function normal_edge(c, e) result(e0)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e        !! the edge

    e0 = e
    do while (e0 /= 0)
        if (c%ed(e0)%etype == NORMAL) exit
        e0 = c%ed(e0)%to_orig
    end do
end function normal_edge

!> Whether `v` can be merged: a plain virtual node, with an edge on each side.
pure logical function candidate(c, v)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: v        !! the node

    candidate = c%nd(v)%kind == VIRTUAL .and. c%nd(v)%in%n == 1 .and. c%nd(v)%out%n == 1 .and. c%nd(v)%label == 0
end function candidate

!> Whether `u` and `v` come from the same node.
pure logical function bothdowncandidates(c, u, v)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: u        !! one node
    integer, intent(in) :: v        !! the other node
    integer :: e, f

    bothdowncandidates = .false.
    if (.not. candidate(c, v)) return
    e = c%nd(u)%in%a(1)
    f = c%nd(v)%in%a(1)
    if (c%ed(e)%tail == c%ed(f)%tail) bothdowncandidates = samedir(c, e, f)
end function bothdowncandidates

!> Whether `u` and `v` go to the same node.
pure logical function bothupcandidates(c, u, v)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: u        !! one node
    integer, intent(in) :: v        !! the other node
    integer :: e, f

    bothupcandidates = .false.
    if (.not. candidate(c, v)) return
    e = c%nd(u)%out%a(1)
    f = c%nd(v)%out%a(1)
    if (c%ed(e)%head == c%ed(f)%head) bothupcandidates = samedir(c, e, f)
end function bothupcandidates

!> Merge the nodes of rank `r` from `lpos` to `rpos` into the first of them.
pure subroutine mergevirtual(c, r, lpos, rpos, dir)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: r            !! the rank
    integer, intent(in) :: lpos         !! the first position
    integer, intent(in) :: rpos         !! the last position
    integer, intent(in) :: dir          !! `UP` or `DOWN`
    integer :: left, right, i, k, e, f, e0, n

    left = c%rk(r)%v(lpos)
    do i = lpos + 1, rpos
        right = c%rk(r)%v(i)
        if (dir == DOWN) then
            do while (c%nd(right)%out%n > 0)
                e = c%nd(right)%out%a(1)
                f = 0
                do k = 1, c%nd(left)%out%n
                    if (c%ed(c%nd(left)%out%a(k))%head == c%ed(e)%head) then
                        f = c%nd(left)%out%a(k)
                        exit
                    end if
                end do
                if (f == 0) call virtual_edge(c, left, c%ed(e)%head, e, f)
                do while (c%nd(right)%in%n > 0)
                    e0 = c%nd(right)%in%a(1)
                    call merge_oneway(c, e0, f)
                    call delete_fast_edge(c, e0)
                end do
                call delete_fast_edge(c, e)
            end do
        else
            do while (c%nd(right)%in%n > 0)
                e = c%nd(right)%in%a(1)
                f = 0
                do k = 1, c%nd(left)%in%n
                    if (c%ed(c%nd(left)%in%a(k))%tail == c%ed(e)%tail) then
                        f = c%nd(left)%in%a(k)
                        exit
                    end if
                end do
                if (f == 0) call virtual_edge(c, c%ed(e)%tail, left, e, f)
                do while (c%nd(right)%out%n > 0)
                    e0 = c%nd(right)%out%a(1)
                    call merge_oneway(c, e0, f)
                    call delete_fast_edge(c, e0)
                end do
                call delete_fast_edge(c, e)
            end do
        end if
        c%nd(right)%alive = .false.
    end do
    k = lpos + 1
    do i = rpos + 1, c%rk(r)%n - 1
        n = c%rk(r)%v(i)
        c%rk(r)%v(k) = n
        c%nd(n)%order = k
        k = k + 1
    end do
    c%rk(r)%n = k
    c%rk(r)%v(k) = 0
end subroutine mergevirtual

!> Merge edges, looking down and then up the ranks.
pure subroutine dot_concentrate(c)
    type(core_t), intent(inout) :: c    !! the state
    integer :: r, leftpos, rightpos, left

    if (c%maxrank - c%minrank <= 1) return
    r = 1
    do while (c%rk(r + 1)%n > 0)
        leftpos = 0
        do while (leftpos < c%rk(r)%n)
            left = c%rk(r)%v(leftpos)
            if (candidate(c, left)) then
                rightpos = leftpos + 1
                do while (rightpos < c%rk(r)%n)
                    if (.not. bothdowncandidates(c, left, c%rk(r)%v(rightpos))) exit
                    rightpos = rightpos + 1
                end do
                if (rightpos - leftpos > 1) call mergevirtual(c, r, leftpos, rightpos - 1, DOWN)
            end if
            leftpos = leftpos + 1
        end do
        r = r + 1
    end do
    do while (r > 0)
        leftpos = 0
        do while (leftpos < c%rk(r)%n)
            left = c%rk(r)%v(leftpos)
            if (candidate(c, left)) then
                rightpos = leftpos + 1
                do while (rightpos < c%rk(r)%n)
                    if (.not. bothupcandidates(c, left, c%rk(r)%v(rightpos))) exit
                    rightpos = rightpos + 1
                end do
                if (rightpos - leftpos > 1) call mergevirtual(c, r, leftpos, rightpos - 1, UP)
            end if
            leftpos = leftpos + 1
        end do
        r = r - 1
    end do
    c%nlist = pack(c%nlist, c%nd(c%nlist)%alive)
end subroutine dot_concentrate

end module fdot_concentrate
