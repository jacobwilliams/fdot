module fdot_path
!! Shortest paths inside a polygon, and fitting a spline to them. This is a port of
!! `lib/pathplan/shortest.c`, `route.c`, `solvers.c`, `util.c` and part of `triang.c` of Graphviz, which
!! are under the Eclipse Public License 2.0. What is static there is passed through the calls here.
use fdot_structure, only: wp, point_t
implicit none
private

public :: plist_t, add_point, shortest_path, route_spline

real(wp), parameter :: EPSILON1 = 1.0e-3_wp   !! how short a path can be, and how near a barrier a crossing can be
real(wp), parameter :: EPSILON2 = 1.0e-6_wp   !! how near the ends of a curve a crossing is ignored
real(wp), parameter :: PI = acos(-1.0_wp)     !! pi
integer, parameter :: ISCCW = 1     !! three points turn counterclockwise
integer, parameter :: ISCW = 2      !! three points turn clockwise
integer, parameter :: ISON = 3      !! three points are on a line

type :: plist_t
    !! A list of points that grows.
    integer :: n = 0                    !! how many there are
    type(point_t), allocatable :: a(:)  !! the points
end type plist_t

type :: triangle_t
    !! A triangle of the polygon. Its corners are numbers of points.
    integer :: mark = 0     !! where it is in the path
    integer :: p0(3) = 0    !! where each edge starts
    integer :: p1(3) = 0    !! where each edge ends
    integer :: right(3) = 0 !! the triangle on the other side of each edge, or 0
end type triangle_t

contains

!> Append a point.
subroutine add_point(list, p)
    type(plist_t), intent(inout) :: list    !! the list
    type(point_t), value :: p               !! the point, a copy as it can be an element of `list`
    type(point_t), allocatable :: bigger(:)

    if (.not. allocated(list%a)) then
        allocate (list%a(16))
    else if (list%n == size(list%a)) then
        allocate (bigger(2*list%n))
        bigger(1:list%n) = list%a(1:list%n)
        call move_alloc(bigger, list%a)
    end if
    list%n = list%n + 1
    list%a(list%n) = p
end subroutine add_point

!> Which way three points turn.
integer function ccw(p1, p2, p3)
    type(point_t), intent(in) :: p1 !! the first point
    type(point_t), intent(in) :: p2 !! the second point
    type(point_t), intent(in) :: p3 !! the third point
    real(wp) :: d

    d = (p1%y - p2%y)*(p3%x - p2%x) - (p3%y - p2%y)*(p1%x - p2%x)
    if (d > 0) then
        ccw = ISCW
    else if (d < 0) then
        ccw = ISCCW
    else
        ccw = ISON
    end if
end function ccw

!> Whether `pb` is between `pa` and `pc`.
logical function between(pa, pb, pc)
    type(point_t), intent(in) :: pa !! one end
    type(point_t), intent(in) :: pb !! the point
    type(point_t), intent(in) :: pc !! the other end
    real(wp) :: bax, bay, cax, cay

    between = .false.
    if (ccw(pa, pb, pc) /= ISON) return
    bax = pb%x - pa%x
    bay = pb%y - pa%y
    cax = pc%x - pa%x
    cay = pc%y - pa%y
    between = cax*bax + cay*bay >= 0 .and. cax*cax + cay*cay <= bax*bax + bay*bay
end function between

!> Whether two segments intersect.
logical function intersects(pa, pb, pc, pd)
    type(point_t), intent(in) :: pa !! the start of one
    type(point_t), intent(in) :: pb !! the end of one
    type(point_t), intent(in) :: pc !! the start of the other
    type(point_t), intent(in) :: pd !! the end of the other
    logical :: c1, c2, c3, c4

    if (ccw(pa, pb, pc) == ISON .or. ccw(pa, pb, pd) == ISON .or. ccw(pc, pd, pa) == ISON .or. ccw(pc, pd, pb) == ISON) then
        intersects = between(pa, pb, pc) .or. between(pa, pb, pd) .or. between(pc, pd, pa) .or. between(pc, pd, pb)
    else
        c1 = ccw(pa, pb, pc) == ISCCW
        c2 = ccw(pa, pb, pd) == ISCCW
        c3 = ccw(pc, pd, pa) == ISCCW
        c4 = ccw(pc, pd, pb) == ISCCW
        intersects = (c1 .neqv. c2) .and. (c3 .neqv. c4)
    end if
end function intersects

!> Whether the segment from point `i` to point `i + 2` of the polygon is a diagonal. Positions are
!> counted from 0, and `points` holds numbers of points of `p`.
logical function is_diagonal(i, p, points, count)
    integer, intent(in) :: i            !! the position
    type(point_t), intent(in) :: p(:)   !! the points
    integer, intent(in) :: points(0:)   !! the polygon
    integer, intent(in) :: count        !! how many points it has
    integer :: ip1, ip2, im1, j, jp1
    logical :: res

    ip1 = mod(i + 1, count)
    ip2 = mod(i + 2, count)
    im1 = mod(i + count - 1, count)
    if (ccw(p(points(im1)), p(points(i)), p(points(ip1))) == ISCCW) then
        res = ccw(p(points(i)), p(points(ip2)), p(points(im1))) == ISCCW
        if (res) res = ccw(p(points(ip2)), p(points(i)), p(points(ip1))) == ISCCW
    else
        res = ccw(p(points(i)), p(points(ip2)), p(points(ip1))) == ISCW
    end if
    is_diagonal = .false.
    if (.not. res) return
    do j = 0, count - 1
        jp1 = mod(j + 1, count)
        if (.not. (j == i .or. jp1 == i .or. j == ip2 .or. jp1 == ip2)) then
            if (intersects(p(points(i)), p(points(ip2)), p(points(j)), p(points(jp1)))) return
        end if
    end do
    is_diagonal = .true.
end function is_diagonal

!> Cut the polygon into triangles.
subroutine triangulate(p, points, count, tris, ntri)
    type(point_t), intent(in) :: p(:)       !! the points
    integer, intent(inout) :: points(0:)    !! the polygon, which this uses up
    integer, intent(in) :: count            !! how many points it has
    type(triangle_t), intent(inout) :: tris(:)  !! the triangles
    integer, intent(out) :: ntri            !! how many there are
    integer :: left, i, ip1, ip2
    logical :: found

    ntri = 0
    left = count
    do while (left > 3)
        found = .false.
        do i = 0, left - 1
            ip1 = mod(i + 1, left)
            ip2 = mod(i + 2, left)
            if (is_diagonal(i, p, points, left)) then
                call load_triangle(tris, ntri, points(i), points(ip1), points(ip2))
                points(ip1:left - 2) = points(ip1 + 1:left - 1)
                found = .true.
                exit
            end if
        end do
        if (.not. found) return
        left = left - 1
    end do
    call load_triangle(tris, ntri, points(0), points(1), points(2))
end subroutine triangulate

!> Add a triangle.
subroutine load_triangle(tris, ntri, a, b, c)
    type(triangle_t), intent(inout) :: tris(:)  !! the triangles
    integer, intent(inout) :: ntri              !! how many there are
    integer, intent(in) :: a                    !! the first corner
    integer, intent(in) :: b                    !! the second corner
    integer, intent(in) :: c                    !! the third corner

    ntri = ntri + 1
    tris(ntri)%p0 = [a, b, c]
    tris(ntri)%p1 = [b, c, a]
end subroutine load_triangle

!> Connect two triangles at their common edge, if they have one.
subroutine connect_tris(tris, t1, t2)
    type(triangle_t), intent(inout) :: tris(:)  !! the triangles
    integer, intent(in) :: t1                   !! one triangle
    integer, intent(in) :: t2                   !! the other triangle
    integer :: ei, ej

    do ei = 1, 3
        do ej = 1, 3
            if ((tris(t1)%p0(ei) == tris(t2)%p0(ej) .and. tris(t1)%p1(ei) == tris(t2)%p1(ej)) .or. &
                (tris(t1)%p0(ei) == tris(t2)%p1(ej) .and. tris(t1)%p1(ei) == tris(t2)%p0(ej))) then
                tris(t1)%right(ei) = t2
                tris(t2)%right(ej) = t1
            end if
        end do
    end do
end subroutine connect_tris

!> Mark the path of triangles from `i` to `j`.
recursive function mark_path(tris, i, j) result(found)
    type(triangle_t), intent(inout) :: tris(:)  !! the triangles
    integer, intent(in) :: i                    !! where to start
    integer, intent(in) :: j                    !! where to end
    logical :: found
    integer :: ei

    found = .false.
    if (tris(i)%mark /= 0) return
    tris(i)%mark = 1
    found = i == j
    if (found) return
    do ei = 1, 3
        if (tris(i)%right(ei) /= 0) then
            found = mark_path(tris, tris(i)%right(ei), j)
            if (found) return
        end if
    end do
    tris(i)%mark = 0
end function mark_path

!> Whether point `q` is in triangle `t`.
logical function point_in_tri(t, p, q)
    type(triangle_t), intent(in) :: t   !! the triangle
    type(point_t), intent(in) :: p(:)   !! the points
    type(point_t), intent(in) :: q      !! the point
    integer :: ei, sum

    sum = 0
    do ei = 1, 3
        if (ccw(p(t%p0(ei)), p(t%p1(ei)), q) /= ISCW) sum = sum + 1
    end do
    point_in_tri = sum == 3 .or. sum == 0
end function point_in_tri

!> Find a shortest path inside the polygon `poly`, from `eps(1)` to `eps(2)`. `stat` is not 0 if an end
!> is not in the polygon.
subroutine shortest_path(poly, eps, output, stat)
    type(point_t), intent(in) :: poly(:)    !! the polygon
    type(point_t), intent(in) :: eps(2)     !! the ends of the path
    type(plist_t), intent(out) :: output    !! the path
    integer, intent(out) :: stat            !! the result
    type(point_t), allocatable :: p(:)
    type(triangle_t), allocatable :: tris(:)
    integer, allocatable :: order(:), link(:), dq(:), reversed(:)
    integer :: pn, count, i, ntri, ti, tj, ftri, ltri, minpi, ep0, ep1, front, back, apex, ei, q, l, r, split, next
    real(wp) :: minx
    type(point_t) :: p1, p2, p3

    stat = 0
    pn = size(poly)
    allocate (p(pn + 2), order(0:pn), link(pn + 2), dq(0:2*pn + 8))
    link = 0
    minx = huge(minx)
    minpi = 1
    do i = 1, pn
        if (minx > poly(i)%x) then
            minx = poly(i)%x
            minpi = i
        end if
    end do
    p2 = poly(minpi)
    p1 = poly(merge(pn, minpi - 1, minpi == 1))
    p3 = poly(mod(minpi, pn) + 1)
    count = 0
    if ((p1%x == p2%x .and. p2%x == p3%x .and. p3%y > p2%y) .or. ccw(p1, p2, p3) /= ISCCW) then
        do i = pn, 1, -1
            if (i < pn) then
                if (poly(i)%x == poly(i + 1)%x .and. poly(i)%y == poly(i + 1)%y) cycle
            end if
            count = count + 1
            p(count) = poly(i)
            order(count - 1) = count
        end do
    else
        do i = 1, pn
            if (i > 1) then
                if (poly(i)%x == poly(i - 1)%x .and. poly(i)%y == poly(i - 1)%y) cycle
            end if
            count = count + 1
            p(count) = poly(i)
            order(count - 1) = count
        end do
    end if
    ep0 = count + 1
    ep1 = count + 2
    p(ep0) = eps(1)
    p(ep1) = eps(2)

    allocate (tris(count))
    call triangulate(p, order, count, tris, ntri)
    do ti = 1, ntri
        do tj = ti + 1, ntri
            call connect_tris(tris, ti, tj)
        end do
    end do
    ftri = 0
    do ti = 1, ntri
        if (point_in_tri(tris(ti), p, eps(1))) then
            ftri = ti
            exit
        end if
    end do
    ltri = 0
    do ti = 1, ntri
        if (point_in_tri(tris(ti), p, eps(2))) then
            ltri = ti
            exit
        end if
    end do
    if (ftri == 0 .or. ltri == 0) then
        stat = 1
        return
    end if
    if (.not. mark_path(tris, ftri, ltri) .or. ftri == ltri) then
        ! a straight line is better than failing
        call add_point(output, eps(1))
        call add_point(output, eps(2))
        return
    end if

    ! build the funnel, and the shortest path as links
    front = pn + 4
    back = front - 1
    call add_front(ep0)
    apex = front
    ti = ftri
    do while (ti /= 0)
        tris(ti)%mark = 2
        ei = exit_edge(tris, ti)
        if (ei == 0) then
            if (ccw(eps(2), p(dq(front)), p(dq(back))) == ISCCW) then
                l = dq(back)
                r = ep1
            else
                l = ep1
                r = dq(back)
            end if
        else
            q = tris(ti)%p1(mod(ei, 3) + 1)
            if (ccw(p(tris(ti)%p0(ei)), p(q), p(tris(ti)%p1(ei))) == ISCCW) then
                l = tris(ti)%p1(ei)
                r = tris(ti)%p0(ei)
            else
                l = tris(ti)%p0(ei)
                r = tris(ti)%p1(ei)
            end if
        end if
        if (ti == ftri) then
            call add_back(l)
            call add_front(r)
        else if (dq(front) /= r .and. dq(back) /= r) then
            split = find_split(r)
            front = split
            call add_front(r)
            if (split > apex) apex = split
        else
            split = find_split(l)
            back = split
            call add_back(l)
            if (split < apex) apex = split
        end if
        next = 0
        if (ei /= 0) next = tris(ti)%right(ei)
        ti = next
    end do

    i = 0
    q = ep1
    do while (q /= 0)
        i = i + 1
        q = link(q)
    end do
    allocate (reversed(i))
    q = ep1
    do while (q /= 0)
        reversed(i) = q
        i = i - 1
        q = link(q)
    end do
    do i = 1, size(reversed)
        call add_point(output, p(reversed(i)))
    end do

contains

    !> Add a point at the front of the funnel.
    subroutine add_front(point)
        integer, intent(in) :: point    !! the point

        if (back >= front) link(point) = dq(front)
        front = front - 1
        dq(front) = point
    end subroutine add_front

    !> Add a point at the back of the funnel.
    subroutine add_back(point)
        integer, intent(in) :: point    !! the point

        if (back >= front) link(point) = dq(back)
        back = back + 1
        dq(back) = point
    end subroutine add_back

    !> Where the funnel must split for a new point.
    integer function find_split(point) result(index)
        integer, intent(in) :: point    !! the point

        do index = front, apex - 1
            if (ccw(p(dq(index + 1)), p(dq(index)), p(point)) == ISCCW) return
        end do
        do index = back, apex + 1, -1
            if (ccw(p(dq(index - 1)), p(dq(index)), p(point)) == ISCW) return
        end do
        index = apex
    end function find_split
end subroutine shortest_path

!> The edge of triangle `ti` that leads to the next triangle of the path, or 0.
integer function exit_edge(tris, ti) result(ei)
    type(triangle_t), intent(in) :: tris(:) !! the triangles
    integer, intent(in) :: ti               !! the triangle

    do ei = 1, 3
        if (tris(ti)%right(ei) /= 0) then
            if (tris(tris(ti)%right(ei))%mark == 1) return
        end if
    end do
    ei = 0
end function exit_edge

!> Fit a spline to `route`, staying inside the barriers.
subroutine route_spline(barriers_a, barriers_b, route, slope0, slope1, output)
    type(point_t), intent(in) :: barriers_a(:)  !! where each barrier starts
    type(point_t), intent(in) :: barriers_b(:)  !! where each barrier ends
    type(point_t), intent(in) :: route(:)       !! the path to follow
    type(point_t), intent(in) :: slope0         !! the direction at the start, or zero
    type(point_t), intent(in) :: slope1         !! the direction at the end, or zero
    type(plist_t), intent(out) :: output        !! the control points

    call add_point(output, route(1))
    call really_route(barriers_a, barriers_b, route, normv(slope0), normv(slope1), output)
end subroutine route_spline

!> Fit a curve to the whole path, or else split the path where it is furthest from the curve.
recursive subroutine really_route(ea, eb, inps, ev0, ev1, output)
    type(point_t), intent(in) :: ea(:)      !! where each barrier starts
    type(point_t), intent(in) :: eb(:)      !! where each barrier ends
    type(point_t), intent(in) :: inps(:)    !! the path
    type(point_t), intent(in) :: ev0        !! the direction at the start
    type(point_t), intent(in) :: ev1        !! the direction at the end
    type(plist_t), intent(inout) :: output  !! the control points
    real(wp), allocatable :: t(:)
    type(point_t), allocatable :: a0(:), a1(:)
    type(point_t) :: p1, p2, v1, v2, cp1, cp2, q, splitv
    real(wp) :: d, maxd
    integer :: i, n, maxi
    logical :: fit

    n = size(inps)
    allocate (t(n), a0(n), a1(n))
    t(1) = 0
    do i = 2, n
        t(i) = t(i - 1) + dist(inps(i), inps(i - 1))
    end do
    if (t(n) > 0) t(2:n) = t(2:n)/t(n)
    do i = 1, n
        a0(i) = times(ev0, b1(t(i)))
        a1(i) = times(ev1, b2(t(i)))
    end do
    call mkspline(inps, t, a0, a1, ev0, ev1, p1, v1, p2, v2)
    call spline_fits(ea, eb, p1, v1, p2, v2, inps, output, fit)
    if (fit) return

    cp1 = add(p1, times(v1, 1/3.0_wp))
    cp2 = sub(p2, times(v2, 1/3.0_wp))
    maxi = -1
    maxd = -1
    do i = 2, n - 1
        q = point_t(b0(t(i))*p1%x + b1(t(i))*cp1%x + b2(t(i))*cp2%x + b3(t(i))*p2%x, &
                    b0(t(i))*p1%y + b1(t(i))*cp1%y + b2(t(i))*cp2%y + b3(t(i))*p2%y)
        d = dist(q, inps(i))
        if (d > maxd) then
            maxd = d
            maxi = i
        end if
    end do
    splitv = normv(add(normv(sub(inps(maxi), inps(maxi - 1))), normv(sub(inps(maxi + 1), inps(maxi)))))
    call really_route(ea, eb, inps(1:maxi), ev0, splitv, output)
    call really_route(ea, eb, inps(maxi:n), splitv, ev1, output)
end subroutine really_route

!> The curve from the first to the last point that best fits the points, by least squares.
subroutine mkspline(inps, t, a0, a1, ev0, ev1, sp0, sv0, sp1, sv1)
    type(point_t), intent(in) :: inps(:)    !! the path
    real(wp), intent(in) :: t(:)            !! where each point is along the path
    type(point_t), intent(in) :: a0(:)      !! the effect of the first scale on each point
    type(point_t), intent(in) :: a1(:)      !! the effect of the second scale on each point
    type(point_t), intent(in) :: ev0        !! the direction at the start
    type(point_t), intent(in) :: ev1        !! the direction at the end
    type(point_t), intent(out) :: sp0       !! the start
    type(point_t), intent(out) :: sv0       !! the tangent at the start
    type(point_t), intent(out) :: sp1       !! the end
    type(point_t), intent(out) :: sv1       !! the tangent at the end
    type(point_t) :: tmp
    real(wp) :: c00, c01, c11, x0, x1, det01, det0x, detx1, d01, scale0, scale3
    integer :: i, n

    n = size(inps)
    scale0 = 0
    scale3 = 0
    c00 = 0
    c01 = 0
    c11 = 0
    x0 = 0
    x1 = 0
    do i = 1, n
        c00 = c00 + dot(a0(i), a0(i))
        c01 = c01 + dot(a0(i), a1(i))
        c11 = c11 + dot(a1(i), a1(i))
        tmp = sub(inps(i), add(times(inps(1), b01(t(i))), times(inps(n), b23(t(i)))))
        x0 = x0 + dot(a0(i), tmp)
        x1 = x1 + dot(a1(i), tmp)
    end do
    det01 = c00*c11 - c01*c01
    det0x = c00*x1 - c01*x0
    detx1 = x0*c11 - x1*c01
    if (abs(det01) >= 1.0e-6_wp) then
        scale0 = detx1/det01
        scale3 = det0x/det01
    end if
    if (abs(det01) < 1.0e-6_wp .or. scale0 <= 0 .or. scale3 <= 0) then
        d01 = dist(inps(1), inps(n))/3
        scale0 = d01
        scale3 = d01
    end if
    sp0 = inps(1)
    sv0 = times(ev0, scale0)
    sp1 = inps(n)
    sv1 = times(ev1, scale3)
end subroutine mkspline

!> The length of a line through points.
function dist_n(p) result(rv)
    type(point_t), intent(in) :: p(:)   !! the points
    real(wp) :: rv
    integer :: i

    rv = 0
    do i = 2, size(p)
        rv = rv + hypot(p(i)%x - p(i - 1)%x, p(i)%y - p(i - 1)%y)
    end do
end function dist_n

!> Look for a curve that stays inside, with tangents of decreasing length, and keep it if there is one.
subroutine spline_fits(ea, eb, pa, va, pb, vb, inps, output, fit)
    type(point_t), intent(in) :: ea(:)      !! where each barrier starts
    type(point_t), intent(in) :: eb(:)      !! where each barrier ends
    type(point_t), intent(in) :: pa         !! the start
    type(point_t), intent(in) :: va         !! the tangent at the start
    type(point_t), intent(in) :: pb         !! the end
    type(point_t), intent(in) :: vb         !! the tangent at the end
    type(point_t), intent(in) :: inps(:)    !! the path
    type(plist_t), intent(inout) :: output  !! the control points
    logical, intent(out) :: fit             !! whether a curve was found
    type(point_t) :: sps(4)
    real(wp) :: a
    logical :: force, first
    integer :: pi

    fit = .false.
    force = size(inps) == 2
    first = .true.
    a = 4
    do
        sps(1) = pa
        sps(2) = point_t(pa%x + a*va%x/3, pa%y + a*va%y/3)
        sps(3) = point_t(pb%x - a*vb%x/3, pb%y - a*vb%y/3)
        sps(4) = pb
        ! a path shorter than the shortest path must leave the polygon
        if (first) then
            if (dist_n(sps) < dist_n(inps) - EPSILON1) return
        end if
        first = .false.
        if (spline_is_inside(ea, eb, sps)) then
            do pi = 2, 4
                call add_point(output, sps(pi))
            end do
            fit = .true.
            return
        end if
        if (a < 0.005_wp) then
            if (force) then
                do pi = 2, 4
                    call add_point(output, sps(pi))
                end do
                fit = .true.
            end if
            return
        end if
        if (a > 0.01_wp) then
            a = a/2
        else
            a = 0
        end if
    end do
end subroutine spline_fits

!> Whether a curve crosses no barrier.
logical function spline_is_inside(ea, eb, sps)
    type(point_t), intent(in) :: ea(:)  !! where each barrier starts
    type(point_t), intent(in) :: eb(:)  !! where each barrier ends
    type(point_t), intent(in) :: sps(4) !! the curve
    real(wp) :: roots(9), t, ta, tb, tc, td
    type(point_t) :: ip
    integer :: ei, rooti, rootn

    spline_is_inside = .false.
    do ei = 1, size(ea)
        call spline_intersects_line(sps, ea(ei), eb(ei), roots, rootn)
        if (rootn == 4) cycle
        do rooti = 1, rootn
            if (roots(rooti) < EPSILON2 .or. roots(rooti) > 1 - EPSILON2) cycle
            t = roots(rooti)
            td = t*t*t
            tc = 3*t*t*(1 - t)
            tb = 3*t*(1 - t)*(1 - t)
            ta = (1 - t)*(1 - t)*(1 - t)
            ip%x = ta*sps(1)%x + tb*sps(2)%x + tc*sps(3)%x + td*sps(4)%x
            ip%y = ta*sps(1)%y + tb*sps(2)%y + tc*sps(3)%y + td*sps(4)%y
            if (distsq(ip, ea(ei)) < EPSILON1 .or. distsq(ip, eb(ei)) < EPSILON1) cycle
            return
        end do
    end do
    spline_is_inside = .true.
end function spline_is_inside

!> Where a curve meets a segment. A count of 4 means that they overlap.
subroutine spline_intersects_line(sps, l0, l1, roots, rootn)
    type(point_t), intent(in) :: sps(4) !! the curve
    type(point_t), intent(in) :: l0     !! the start of the segment
    type(point_t), intent(in) :: l1     !! the end of the segment
    real(wp), intent(out) :: roots(9)   !! the values of the parameter
    integer, intent(out) :: rootn       !! how many there are
    real(wp) :: scoeff(0:3), xc0, xc1, yc0, yc1, xroots(3), yroots(3), tv, sv, rat
    integer :: xrootn, yrootn, i, j

    xc0 = l0%x
    xc1 = l1%x - l0%x
    yc0 = l0%y
    yc1 = l1%y - l0%y
    rootn = 0
    if (xc1 == 0) then
        if (yc1 == 0) then
            call points2coeff(sps%x, scoeff)
            scoeff(0) = scoeff(0) - xc0
            call solve3(scoeff, xroots, xrootn)
            call points2coeff(sps%y, scoeff)
            scoeff(0) = scoeff(0) - yc0
            call solve3(scoeff, yroots, yrootn)
            if (xrootn == 4) then
                if (yrootn == 4) then
                    rootn = 4
                else
                    do j = 1, yrootn
                        call add_root(yroots(j), roots, rootn)
                    end do
                end if
            else if (yrootn == 4) then
                do i = 1, xrootn
                    call add_root(xroots(i), roots, rootn)
                end do
            else
                do i = 1, xrootn
                    do j = 1, yrootn
                        if (xroots(i) == yroots(j)) call add_root(xroots(i), roots, rootn)
                    end do
                end do
            end if
        else
            call points2coeff(sps%x, scoeff)
            scoeff(0) = scoeff(0) - xc0
            call solve3(scoeff, xroots, xrootn)
            if (xrootn == 4) then
                rootn = 4
                return
            end if
            call points2coeff(sps%y, scoeff)
            do i = 1, xrootn
                tv = xroots(i)
                if (tv >= 0 .and. tv <= 1) then
                    sv = scoeff(0) + tv*(scoeff(1) + tv*(scoeff(2) + tv*scoeff(3)))
                    sv = (sv - yc0)/yc1
                    if (0 <= sv .and. sv <= 1) call add_root(tv, roots, rootn)
                end if
            end do
        end if
    else
        rat = yc1/xc1
        call points2coeff(sps%y - rat*sps%x, scoeff)
        scoeff(0) = scoeff(0) + rat*xc0 - yc0
        call solve3(scoeff, xroots, xrootn)
        if (xrootn == 4) then
            rootn = 4
            return
        end if
        call points2coeff(sps%x, scoeff)
        do i = 1, xrootn
            tv = xroots(i)
            if (tv >= 0 .and. tv <= 1) then
                sv = scoeff(0) + tv*(scoeff(1) + tv*(scoeff(2) + tv*scoeff(3)))
                sv = (sv - xc0)/xc1
                if (0 <= sv .and. sv <= 1) call add_root(tv, roots, rootn)
            end if
        end do
    end if
end subroutine spline_intersects_line

!> The coefficients of the polynomial of one coordinate of a curve.
subroutine points2coeff(v, coeff)
    real(wp), intent(in) :: v(4)            !! the coordinates of the control points
    real(wp), intent(out) :: coeff(0:3)     !! the coefficients

    coeff(3) = v(4) + 3*v(2) - (v(1) + 3*v(3))
    coeff(2) = 3*v(1) + 3*v(3) - 6*v(2)
    coeff(1) = 3*(v(2) - v(1))
    coeff(0) = v(1)
end subroutine points2coeff

!> Keep a root if it is on the curve.
subroutine add_root(root, roots, rootn)
    real(wp), intent(in) :: root        !! the root
    real(wp), intent(inout) :: roots(:) !! the roots
    integer, intent(inout) :: rootn     !! how many there are

    if (root >= 0 .and. root <= 1) then
        rootn = rootn + 1
        roots(rootn) = root
    end if
end subroutine add_root

!> The cube root.
elemental function cbrt(x)
    real(wp), intent(in) :: x   !! the value
    real(wp) :: cbrt

    cbrt = sign(abs(x)**(1/3.0_wp), x)
end function cbrt

!> The roots of a cubic. A count of 4 means that all values are roots.
subroutine solve3(coeff, roots, rootn)
    real(wp), intent(in) :: coeff(0:3)  !! the coefficients
    real(wp), intent(out) :: roots(3)   !! the roots
    integer, intent(out) :: rootn       !! how many there are
    real(wp) :: a, b, c, d, p, q, disc, b_over_3a, c_over_a, d_over_a, r, theta, temp, alpha, beta

    a = coeff(3)
    b = coeff(2)
    c = coeff(1)
    d = coeff(0)
    if (aeq0(a)) then
        call solve2(coeff, roots, rootn)
        return
    end if
    b_over_3a = b/(3*a)
    c_over_a = c/a
    d_over_a = d/a
    p = b_over_3a*b_over_3a
    q = 2*b_over_3a*p - b_over_3a*c_over_a + d_over_a
    p = c_over_a/3 - p
    disc = q*q + 4*p*p*p
    if (disc < 0) then
        r = 0.5_wp*sqrt(-disc + q*q)
        theta = atan2(sqrt(-disc), -q)
        temp = 2*cbrt(r)
        roots(1) = temp*cos(theta/3)
        roots(2) = temp*cos((theta + PI + PI)/3)
        roots(3) = temp*cos((theta - PI - PI)/3)
        rootn = 3
    else
        alpha = 0.5_wp*(sqrt(disc) - q)
        beta = -q - alpha
        roots(1) = cbrt(alpha) + cbrt(beta)
        if (disc > 0) then
            rootn = 1
        else
            roots(2) = -0.5_wp*roots(1)
            roots(3) = roots(2)
            rootn = 3
        end if
    end if
    roots(1:rootn) = roots(1:rootn) - b_over_3a
end subroutine solve3

!> The roots of a quadratic.
subroutine solve2(coeff, roots, rootn)
    real(wp), intent(in) :: coeff(0:3)  !! the coefficients
    real(wp), intent(out) :: roots(3)   !! the roots
    integer, intent(out) :: rootn       !! how many there are
    real(wp) :: a, b, c, disc, b_over_2a, c_over_a

    roots = 0
    a = coeff(2)
    b = coeff(1)
    c = coeff(0)
    if (aeq0(a)) then
        if (aeq0(b)) then
            rootn = merge(4, 0, aeq0(c))
        else
            roots(1) = -c/b
            rootn = 1
        end if
        return
    end if
    b_over_2a = b/(2*a)
    c_over_a = c/a
    disc = b_over_2a*b_over_2a - c_over_a
    if (disc < 0) then
        rootn = 0
    else if (disc > 0) then
        roots(1) = -b_over_2a + sqrt(disc)
        roots(2) = -2*b_over_2a - roots(1)
        rootn = 2
    else
        roots(1) = -b_over_2a
        rootn = 1
    end if
end subroutine solve2

!> Whether `x` is zero, nearly.
logical function aeq0(x)
    real(wp), intent(in) :: x   !! the value

    aeq0 = x < 1.0e-7_wp .and. x > -1.0e-7_wp
end function aeq0

!> The vector `v` with length 1, or `v` if it is too short for that.
function normv(v) result(r)
    type(point_t), intent(in) :: v  !! the vector
    type(point_t) :: r
    real(wp) :: d

    r = v
    d = v%x*v%x + v%y*v%y
    if (d > 1.0e-6_wp) then
        d = sqrt(d)
        r%x = r%x/d
        r%y = r%y/d
    end if
end function normv

!> The sum of two points.
function add(a, b) result(r)
    type(point_t), intent(in) :: a  !! one point
    type(point_t), intent(in) :: b  !! the other point
    type(point_t) :: r

    r = point_t(a%x + b%x, a%y + b%y)
end function add

!> The difference of two points.
function sub(a, b) result(r)
    type(point_t), intent(in) :: a  !! one point
    type(point_t), intent(in) :: b  !! the other point
    type(point_t) :: r

    r = point_t(a%x - b%x, a%y - b%y)
end function sub

!> The distance between two points.
function dist(a, b) result(r)
    type(point_t), intent(in) :: a  !! one point
    type(point_t), intent(in) :: b  !! the other point
    real(wp) :: r

    r = hypot(b%x - a%x, b%y - a%y)
end function dist

!> The square of the distance between two points.
function distsq(a, b) result(r)
    type(point_t), intent(in) :: a  !! one point
    type(point_t), intent(in) :: b  !! the other point
    real(wp) :: r

    r = (a%x - b%x)*(a%x - b%x) + (a%y - b%y)*(a%y - b%y)
end function distsq

!> A point times a number.
function times(a, c) result(r)
    type(point_t), intent(in) :: a  !! the point
    real(wp), intent(in) :: c       !! the number
    type(point_t) :: r

    r = point_t(a%x*c, a%y*c)
end function times

!> The dot product of two points.
function dot(a, b) result(r)
    type(point_t), intent(in) :: a  !! one point
    type(point_t), intent(in) :: b  !! the other point
    real(wp) :: r

    r = a%x*b%x + a%y*b%y
end function dot

!> A Bernstein polynomial.
real(wp) function b0(t)
    real(wp), intent(in) :: t   !! the parameter

    b0 = (1 - t)*(1 - t)*(1 - t)
end function b0

!> A Bernstein polynomial.
real(wp) function b1(t)
    real(wp), intent(in) :: t   !! the parameter

    b1 = 3*t*(1 - t)*(1 - t)
end function b1

!> A Bernstein polynomial.
real(wp) function b2(t)
    real(wp), intent(in) :: t   !! the parameter

    b2 = 3*t*t*(1 - t)
end function b2

!> A Bernstein polynomial.
real(wp) function b3(t)
    real(wp), intent(in) :: t   !! the parameter

    b3 = t*t*t
end function b3

!> The sum of the first two Bernstein polynomials.
real(wp) function b01(t)
    real(wp), intent(in) :: t   !! the parameter

    b01 = (1 - t)*(1 - t)*((1 - t) + 3*t)
end function b01

!> The sum of the last two Bernstein polynomials.
real(wp) function b23(t)
    real(wp), intent(in) :: t   !! the parameter

    b23 = t*t*(3*(1 - t) + t)
end function b23

end module fdot_path
