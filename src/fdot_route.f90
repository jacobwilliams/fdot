module fdot_route
!! Routing edges: boxes between the ranks, a shortest path through them, and a spline along it. This is
!! a port of `lib/dotgen/dotsplines.c`, `lib/common/routespl.c` and parts of `splines.c` of Graphviz,
!! which are under the Eclipse Public License 2.0, without clusters, ports and flat edges.
use fdot_structure
use fdot_path, only: plist_t, add_point, shortest_path, route_spline
implicit none
private

public :: dot_splines

integer, parameter :: MINW = 16                 !! the smallest width of a box of a path
integer, parameter :: HALFMINW = 8              !! half of `MINW`
integer, parameter :: FUDGE = 4                 !! room beside a node
integer, parameter :: REGULAREDGE = 1           !! flags of an edge: one between ranks
integer, parameter :: SELFNPEDGE = 8            !! flags of an edge: a self-loop without ports
integer, parameter :: EDGETYPEMASK = 15         !! the flags of the kind of an edge
integer, parameter :: FWDEDGE = 16              !! flags of an edge: it goes down the ranks
integer, parameter :: BWDEDGE = 32              !! flags of an edge: it goes up the ranks
integer, parameter :: MAINGRAPH = 64            !! flags of an edge: it is in the graph of ranks
integer, parameter :: AUXGRAPH = 128            !! flags of an edge: it is routed apart from that graph
integer, parameter :: GRAPHTYPEMASK = 192       !! the flags of the graph of an edge
integer, parameter :: INIT_DELTA = 10           !! samples of a curve for each box, to shrink the boxes
integer, parameter :: LOOP_TRIES = 15           !! how many times to sample more finely
real(wp), parameter :: SMALL_GAP = 0.0001_wp    !! a tolerance in points, for tests of position
real(wp), parameter :: MILLIPOINT = 0.001_wp    !! the distance at which two points are the same
real(wp), parameter :: PI = acos(-1.0_wp)       !! pi
real(wp), parameter :: ARROW_NOMINAL = 10.0_wp  !! the length of an arrowhead of size 1
real(wp), parameter :: EDGE_PEN = 1.0_wp        !! the width of the lines of edges
real(wp), parameter :: HW = 2.0_wp              !! how straight a part of a curve must be to bound it as it is

type :: end_t
    !! An end of a path.
    type(point_t) :: p              !! the point
    real(wp) :: theta = 0           !! the direction
    logical :: constrained = .false.    !! whether `theta` applies
end type end_t

type :: path_t
    !! A path through boxes.
    type(end_t) :: start                !! where it starts
    type(end_t) :: end                  !! where it ends
    integer :: nbox = 0                 !! how many boxes there are
    type(box_t), allocatable :: boxes(:)    !! the boxes
end type path_t

type :: pathend_t
    !! The boxes near an end of a path.
    type(box_t) :: nb                   !! the box around the node
    integer :: boxn = 0                 !! how many there are
    type(box_t) :: boxes(4)             !! the boxes
end type pathend_t

type :: spline_info_t
    !! What is known about the whole drawing.
    real(wp) :: left_bound = 0              !! the left side of the space for paths
    real(wp) :: right_bound = 0             !! the right side of the space for paths
    real(wp) :: splinesep = 0               !! the space beside a virtual node
    real(wp) :: multisep = 0                !! the space between edges drawn together
    type(box_t), allocatable :: rank_box(:) !! the box of the space below each rank, once it is known
end type spline_info_t

type :: blist_t
    !! A list of boxes that grows.
    integer :: n = 0                    !! how many there are
    type(box_t), allocatable :: a(:)    !! the boxes
end type blist_t

contains

!> Whether `n` is a virtual node where edges meet.
logical function spline_merge(c, n)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: n        !! the node

    spline_merge = c%nd(n)%kind == VIRTUAL .and. (c%nd(n)%in%n > 1 .or. c%nd(n)%out%n > 1)
end function spline_merge

!> Whether the ends of a curve for `e` must be swapped to go from tail to head.
logical function swap_ends(c, e0)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e0       !! the edge
    integer :: e

    e = e0
    do while (c%ed(e)%to_orig /= 0)
        e = c%ed(e)%to_orig
    end do
    swap_ends = .false.
    if (c%nd(c%ed(e)%head)%rank > c%nd(c%ed(e)%tail)%rank) return
    if (c%nd(c%ed(e)%head)%rank < c%nd(c%ed(e)%tail)%rank) then
        swap_ends = .true.
        return
    end if
    swap_ends = c%nd(c%ed(e)%head)%order < c%nd(c%ed(e)%tail)%order
end function swap_ends

!> The main routing code.
subroutine dot_splines(c, stat)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: stat        !! the result
    type(spline_info_t) :: sd
    type(path_t) :: p
    type(ilist_t) :: edges
    integer :: i, j, k, n, e, n_nodes, l, ind, cnt, e0, le0, r, t, reals
    real(wp) :: sizey, upy, dwny, saved

    stat = 0
    reals = c%ne
    allocate (c%spl(reals))
    sd%splinesep = 4    ! a quarter of the separation of nodes, in integers
    sd%multisep = NODESEP
    n_nodes = 0
    do i = c%minrank, c%maxrank
        n_nodes = n_nodes + c%rk(i)%n
        n = c%rk(i)%v(0)
        if (n /= 0) sd%left_bound = min(sd%left_bound, c%nd(n)%x - c%nd(n)%lw)
        if (c%rk(i)%n > 0) then
            n = c%rk(i)%v(c%rk(i)%n - 1)
            if (n /= 0) sd%right_bound = max(sd%right_bound, c%nd(n)%x + c%nd(n)%rw)
        end if
        sd%left_bound = sd%left_bound - MINW
        sd%right_bound = sd%right_bound + MINW
        do j = 0, c%rk(i)%n - 1
            n = c%rk(i)%v(j)
            if (c%nd(n)%kind /= NORMAL .and. .not. spline_merge(c, n)) cycle
            do k = 1, c%nd(n)%out%n
                e = c%nd(n)%out%a(k)
                if (c%ed(e)%etype == IGNORED) cycle
                c%ed(e)%flags = ior(ior(REGULAREDGE, FWDEDGE), MAINGRAPH)
                call push(edges, e)
            end do
            if (c%nd(n)%other%n > 0) then
                if (c%nd(n)%kind == NORMAL) then
                    saved = c%nd(n)%rw
                    c%nd(n)%rw = c%nd(n)%saved_rw
                    c%nd(n)%saved_rw = saved
                end if
                do k = 1, c%nd(n)%other%n
                    e = c%nd(n)%other%a(k)
                    call set_other_flags(c, e)
                    call push(edges, e)
                end do
            end if
        end do
    end do

    call sort_edges(c, edges)

    allocate (p%boxes(n_nodes + 200))
    allocate (sd%rank_box(0:c%maxrank + 1))
    l = 1
    do while (l <= edges%n)
        ind = l
        e0 = edges%a(l)
        le0 = main_edge(c, e0)
        l = l + 1
        cnt = 1
        do while (l <= edges%n)
            if (le0 /= main_edge(c, edges%a(l))) exit
            if (iand(c%ed(edges%a(l))%flags, MAINGRAPH) /= 0) exit
            cnt = cnt + 1
            l = l + 1
        end do

        if (c%ed(e0)%tail == c%ed(e0)%head) then
            n = c%ed(e0)%tail
            r = c%nd(n)%rank
            if (r == c%maxrank .or. c%rk(r + 1)%v(0) == 0) then
                t = 0
                if (r > 0) t = c%rk(r - 1)%v(0)
                if (t /= 0) then
                    sizey = c%nd(t)%y - c%nd(n)%y
                else
                    sizey = c%nd(n)%ht
                end if
            else if (r == c%minrank) then
                sizey = c%nd(n)%y - c%nd(c%rk(r + 1)%v(0))%y
            else if (c%rk(r - 1)%v(0) == 0) then
                sizey = c%nd(n)%y - c%nd(c%rk(r + 1)%v(0))%y
            else
                upy = c%nd(c%rk(r - 1)%v(0))%y - c%nd(n)%y
                dwny = c%nd(n)%y - c%nd(c%rk(r + 1)%v(0))%y
                sizey = min(upy, dwny)
            end if
            call self_right(c, edges%a(ind:ind + cnt - 1), sd%multisep, sizey/2)
            do i = ind, ind + cnt - 1
                e = edges%a(i)
                if (c%labelled(e)) call update_bb_label(c, e)
            end do
        else
            call make_regular_edge(c, sd, p, edges%a(ind:ind + cnt - 1), stat)
            if (stat /= 0) return
        end if
    end do

    do i = 1, size(c%nlist)
        n = c%nlist(i)
        if (c%nd(n)%kind == VIRTUAL .and. c%nd(n)%label /= 0) then
            call place_vnlabel(c, n)
        end if
    end do

    do e = 1, reals
        if (swap_ends(c, e)) call swap_spline(c%spl(e))
    end do
end subroutine dot_splines

!> Put the label of a virtual node beside it.
subroutine place_vnlabel(c, n)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: n            !! the node
    integer :: e
    real(wp) :: width

    if (c%nd(n)%in%n == 0) return
    e = c%nd(n)%label
    width = merge(c%label_h(e), c%label_w(e), c%flip)
    c%label_pos(e) = point_t(c%nd(n)%x + width/2, c%nd(n)%y)
    call update_bb_label(c, e)
end subroutine place_vnlabel

!> Grow the extent of the drawing to hold the label of `e`.
subroutine update_bb_label(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    real(wp) :: width, height

    width = merge(c%label_h(e), c%label_w(e), c%flip)
    height = merge(c%label_w(e), c%label_h(e), c%flip)
    c%bb%llx = min(c%bb%llx, c%label_pos(e)%x - width/2)
    c%bb%urx = max(c%bb%urx, c%label_pos(e)%x + width/2)
    c%bb%lly = min(c%bb%lly, c%label_pos(e)%y - height/2)
    c%bb%ury = max(c%bb%ury, c%label_pos(e)%y + height/2)
end subroutine update_bb_label

!> Set how an edge that is not in the graph of ranks is routed.
subroutine set_other_flags(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    integer :: f1, f2

    if (c%ed(e)%tail == c%ed(e)%head) then
        f1 = SELFNPEDGE
        f2 = FWDEDGE
    else
        f1 = REGULAREDGE
        f2 = merge(FWDEDGE, BWDEDGE, c%nd(c%ed(e)%tail)%rank < c%nd(c%ed(e)%head)%rank)
    end if
    c%ed(e)%flags = ior(ior(f1, f2), AUXGRAPH)
end subroutine set_other_flags

!> Whether edge `a` goes before edge `b`, so that equivalent edges are next to each other.
logical function edge_before(c, e0, e1) result(less)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e0       !! one edge
    integer, intent(in) :: e1       !! the other edge
    integer :: et0, et1, le0, le1, v0, v1
    real(wp) :: t0, t1

    ! a tie keeps the order of the list, as the sort is stable
    less = .false.
    et0 = iand(c%ed(e0)%flags, EDGETYPEMASK)
    et1 = iand(c%ed(e1)%flags, EDGETYPEMASK)
    if (et0 /= et1) then
        less = et0 > et1
        return
    end if
    le0 = main_edge(c, e0)
    le1 = main_edge(c, e1)
    v0 = abs(c%nd(c%ed(le0)%tail)%rank - c%nd(c%ed(le0)%head)%rank)
    v1 = abs(c%nd(c%ed(le1)%tail)%rank - c%nd(c%ed(le1)%head)%rank)
    if (v0 /= v1) then
        less = v0 < v1
        return
    end if
    t0 = abs(c%nd(c%ed(le0)%tail)%x - c%nd(c%ed(le0)%head)%x)
    t1 = abs(c%nd(c%ed(le1)%tail)%x - c%nd(c%ed(le1)%head)%x)
    if (t0 /= t1) then
        less = t0 < t1
        return
    end if
    if (c%ed(le0)%seq /= c%ed(le1)%seq) then
        less = c%ed(le0)%seq < c%ed(le1)%seq
        return
    end if
    et0 = iand(c%ed(e0)%flags, GRAPHTYPEMASK)
    et1 = iand(c%ed(e1)%flags, GRAPHTYPEMASK)
    if (et0 /= et1) then
        less = et0 < et1
        return
    end if
    less = c%ed(e0)%seq < c%ed(e1)%seq
end function edge_before

!> Sort the edges with a stable merge sort, as the same input must give the same output everywhere.
subroutine sort_edges(c, edges)
    type(core_t), intent(in) :: c           !! the state
    type(ilist_t), intent(inout) :: edges   !! the edges
    integer, allocatable :: from(:), to(:)
    integer :: n, width, lo, mid, hi, i, j, k

    n = edges%n
    if (n < 2) return
    allocate (from(n), to(n))
    from = edges%a(1:n)
    width = 1
    do while (width < n)
        lo = 1
        do while (lo <= n)
            mid = min(lo + width, n + 1)
            hi = min(lo + 2*width, n + 1)
            i = lo
            j = mid
            do k = lo, hi - 1
                if (i >= mid) then
                    to(k) = from(j)
                    j = j + 1
                else if (j >= hi) then
                    to(k) = from(i)
                    i = i + 1
                else if (edge_before(c, from(j), from(i))) then
                    to(k) = from(j)
                    j = j + 1
                else
                    to(k) = from(i)
                    i = i + 1
                end if
            end do
            lo = hi
        end do
        call move_alloc(to, edges%a)
        allocate (to(n))
        from = edges%a(1:n)
        width = 2*width
    end do
end subroutine sort_edges

!> Add a box that is not empty to `p`.
subroutine add_box(p, b)
    type(path_t), intent(inout) :: p    !! the path
    type(box_t), intent(in) :: b        !! the box
    type(box_t), allocatable :: bigger(:)

    if (b%llx < b%urx .and. b%lly < b%ury) then
        if (p%nbox == size(p%boxes)) then
            allocate (bigger(2*p%nbox))
            bigger(1:p%nbox) = p%boxes
            call move_alloc(bigger, p%boxes)
        end if
        p%nbox = p%nbox + 1
        p%boxes(p%nbox) = b
    end if
end subroutine add_box

!> Add a box to a list.
subroutine add_to_list(list, b)
    type(blist_t), intent(inout) :: list    !! the list
    type(box_t), intent(in) :: b            !! the box
    type(box_t), allocatable :: bigger(:)

    if (.not. allocated(list%a)) then
        allocate (list%a(16))
    else if (list%n == size(list%a)) then
        allocate (bigger(2*list%n))
        bigger(1:list%n) = list%a
        call move_alloc(bigger, list%a)
    end if
    list%n = list%n + 1
    list%a(list%n) = b
end subroutine add_to_list

!> The average direction of the edges that meet at `n`.
function conc_slope(c, n) result(slope)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: n        !! the node
    real(wp) :: slope, s_in, s_out, x1, y1, x2, y2
    integer :: i

    s_in = 0
    s_out = 0
    do i = 1, c%nd(n)%in%n
        s_in = s_in + c%nd(c%ed(c%nd(n)%in%a(i))%tail)%x
    end do
    do i = 1, c%nd(n)%out%n
        s_out = s_out + c%nd(c%ed(c%nd(n)%out%a(i))%head)%x
    end do
    x1 = c%nd(n)%x - s_in/c%nd(n)%in%n
    y1 = c%nd(n)%y - c%nd(c%ed(c%nd(n)%in%a(1))%tail)%y
    x2 = s_out/c%nd(n)%out%n - c%nd(n)%x
    y2 = c%nd(c%ed(c%nd(n)%out%a(1))%head)%y - c%nd(n)%y
    slope = (atan2(y1, x1) + atan2(y2, x2))/2
end function conc_slope

!> Set up the box at the tail of a path.
subroutine begin_path(c, p, n, endp, merged)
    type(core_t), intent(in) :: c       !! the state
    type(path_t), intent(inout) :: p    !! the path
    integer, intent(in) :: n            !! the node at the start
    type(pathend_t), intent(inout) :: endp  !! the end of the path
    logical, intent(in) :: merged       !! whether edges meet at `n`

    p%start%p = point_t(c%nd(n)%x, c%nd(n)%y)
    p%start%constrained = merged
    if (merged) p%start%theta = conc_slope(c, n)
    p%nbox = 0
    endp%boxes(1) = endp%nb
    endp%boxn = 1
    endp%boxes(1)%ury = p%start%p%y
    p%start%p%y = p%start%p%y - 1
end subroutine begin_path

!> Set up the box at the head of a path.
subroutine end_path(c, p, n, endp, merged)
    type(core_t), intent(in) :: c       !! the state
    type(path_t), intent(inout) :: p    !! the path
    integer, intent(in) :: n            !! the node at the end
    type(pathend_t), intent(inout) :: endp  !! the end of the path
    logical, intent(in) :: merged       !! whether edges meet at `n`

    p%end%p = point_t(c%nd(n)%x, c%nd(n)%y)
    p%end%constrained = merged
    if (merged) p%end%theta = conc_slope(c, n) + PI
    endp%boxes(1) = endp%nb
    endp%boxn = 1
    endp%boxes(1)%lly = p%end%p%y
    p%end%p%y = p%end%p%y + 1
end subroutine end_path

!> Add a box to an end, if it is not empty.
subroutine end_box(endp, b)
    type(pathend_t), intent(inout) :: endp  !! the end of the path
    type(box_t), intent(in) :: b            !! the box

    if (b%llx < b%urx .and. b%lly < b%ury) then
        endp%boxn = endp%boxn + 1
        endp%boxes(endp%boxn) = b
    end if
end subroutine end_box

!> Route the edges `edges`, which are drawn the same, the first of which is in the graph of ranks.
subroutine make_regular_edge(c, sp, p, edges, stat)
    type(core_t), intent(inout) :: c    !! the state
    type(spline_info_t), intent(inout) :: sp    !! what is known about the drawing
    type(path_t), intent(inout) :: p    !! the path, which is kept to reuse its room
    integer, intent(in) :: edges(:)     !! the edges
    integer, intent(out) :: stat        !! the result
    type(pathend_t) :: tend, hend
    type(blist_t) :: boxes
    type(plist_t) :: pointfs, ps, copy
    type(box_t) :: b
    integer :: e, fe, tn, hn, sl, segfirst, j, k, count_edges, first_out
    logical :: smode, si
    real(wp) :: dx

    stat = 0
    sl = 0
    e = edges(1)
    fe = e
    segfirst = e
    tn = c%ed(e)%tail
    hn = c%ed(e)%head
    b = maximal_bbox(c, sp, tn, 0, e)
    tend%nb = b
    call begin_path(c, p, tn, tend, spline_merge(c, tn))
    b%ury = tend%boxes(tend%boxn)%ury
    b%lly = tend%boxes(tend%boxn)%lly
    call end_box(tend, regular_end_bottom(b, c%nd(tn)%y - c%rk(c%nd(tn)%rank)%ht1))
    smode = .false.
    si = .false.
    do while (c%nd(hn)%kind == VIRTUAL .and. .not. spline_merge(c, hn))
        call rank_box(c, sp, c%nd(tn)%rank, b)
        call add_to_list(boxes, b)
        if (.not. smode) then
            sl = straight_len(c, hn)
            if (sl >= merge(5, 3, c%has_labels)) then
                smode = .true.
                si = .true.
                sl = sl - 2
            end if
        end if
        if ((.not. smode .or. si) .and. c%nd(hn)%out%n > 0) then
            si = .false.
            first_out = c%nd(hn)%out%a(1)
            call add_to_list(boxes, maximal_bbox(c, sp, hn, e, first_out))
            e = first_out
            tn = c%ed(e)%tail
            hn = c%ed(e)%head
            cycle
        end if
        first_out = c%nd(hn)%out%a(1)
        hend%nb = maximal_bbox(c, sp, hn, e, first_out)
        call end_path(c, p, hn, hend, spline_merge(c, c%ed(e)%head))
        call end_box(hend, regular_end_top(hend%boxes(hend%boxn), c%nd(hn)%y + c%rk(c%nd(hn)%rank)%ht2))
        p%end%theta = PI/2
        p%end%constrained = .true.
        call complete_regular_path(p, tend, hend, boxes)
        call routesplines(p, ps, stat)
        if (stat /= 0) return
        do k = 1, ps%n
            call add_point(pointfs, ps%a(k))
        end do
        e = straight_path(c, first_out, sl, pointfs)
        call recover_slack(c, segfirst, p)
        segfirst = e
        tn = c%ed(e)%tail
        hn = c%ed(e)%head
        boxes%n = 0
        tend%nb = maximal_bbox(c, sp, tn, c%nd(tn)%in%a(1), e)
        call begin_path(c, p, tn, tend, spline_merge(c, tn))
        call end_box(tend, regular_end_bottom(tend%boxes(tend%boxn), c%nd(tn)%y - c%rk(c%nd(tn)%rank)%ht1))
        p%start%theta = -PI/2
        p%start%constrained = .true.
        smode = .false.
    end do
    call rank_box(c, sp, c%nd(tn)%rank, b)
    call add_to_list(boxes, b)
    b = maximal_bbox(c, sp, hn, e, 0)
    hend%nb = b
    call end_path(c, p, hn, hend, spline_merge(c, c%ed(e)%head))
    b%ury = hend%boxes(hend%boxn)%ury
    b%lly = hend%boxes(hend%boxn)%lly
    call end_box(hend, regular_end_top(b, c%nd(hn)%y + c%rk(c%nd(hn)%rank)%ht2))
    call complete_regular_path(p, tend, hend, boxes)
    call routesplines(p, ps, stat)
    if (stat /= 0) return
    do k = 1, ps%n
        call add_point(pointfs, ps%a(k))
    end do
    call recover_slack(c, segfirst, p)
    hn = c%ed(e)%head

    ! one copy of the points for each edge
    count_edges = size(edges)
    if (count_edges == 1) then
        call clip_and_install(c, fe, hn, pointfs%a(1:pointfs%n))
        return
    end if
    dx = sp%multisep*(count_edges - 1)/2
    pointfs%a(2:pointfs%n - 1)%x = pointfs%a(2:pointfs%n - 1)%x - dx
    copy = pointfs
    call clip_and_install(c, fe, hn, copy%a(1:copy%n))
    do j = 2, count_edges
        e = edges(j)
        pointfs%a(2:pointfs%n - 1)%x = pointfs%a(2:pointfs%n - 1)%x + sp%multisep
        copy = pointfs
        if (iand(c%ed(e)%flags, BWDEDGE) /= 0) then
            call new_edge(c, c%ed(e)%head, c%ed(e)%tail, k)
            c%ed(k)%to_orig = e
            call clip_and_install(c, k, c%ed(k)%head, copy%a(1:copy%n))
        else
            call clip_and_install(c, e, c%ed(e)%head, copy%a(1:copy%n))
        end if
    end do
end subroutine make_regular_edge

!> Gather the boxes of a path, in order.
subroutine complete_regular_path(p, tendp, hendp, boxes)
    type(path_t), intent(inout) :: p        !! the path
    type(pathend_t), intent(in) :: tendp    !! the end at the tail
    type(pathend_t), intent(in) :: hendp    !! the end at the head
    type(blist_t), intent(in) :: boxes      !! the boxes between
    integer :: i, fb, lb

    do i = 1, tendp%boxn
        call add_box(p, tendp%boxes(i))
    end do
    ! counted from 0, as in the original
    fb = p%nbox + 1
    lb = fb + boxes%n - 3
    do i = 1, boxes%n
        call add_box(p, boxes%a(i))
    end do
    do i = hendp%boxn, 1, -1
        call add_box(p, hendp%boxes(i))
    end do
    call adjust_regular_path(p, fb, lb)
end subroutine complete_regular_path

!> The box between a node and the bottom of its rank.
function regular_end_bottom(b, y) result(r)
    type(box_t), intent(in) :: b    !! the box of the node
    real(wp), intent(in) :: y       !! the bottom of the rank
    type(box_t) :: r

    r = box_t(b%llx, y, b%urx, b%lly)
end function regular_end_bottom

!> The box between a node and the top of its rank.
function regular_end_top(b, y) result(r)
    type(box_t), intent(in) :: b    !! the box of the node
    real(wp), intent(in) :: y       !! the top of the rank
    type(box_t) :: r

    r = box_t(b%llx, b%ury, b%urx, y)
end function regular_end_top

!> Make sure the path is wide enough. `fb` and `lb` are the first and last boxes between ranks.
subroutine adjust_regular_path(p, fb, lb)
    type(path_t), intent(inout) :: p    !! the path
    integer, intent(in) :: fb           !! the first, counted from 0
    integer, intent(in) :: lb           !! the last, counted from 0
    real(wp) :: x
    integer :: i

    do i = fb - 1, lb
        associate (b => p%boxes(i + 1))
            if (modulo(i - fb, 2) == 0) then
                if (b%llx >= b%urx) then
                    x = (b%llx + b%urx)/2
                    b%llx = x - HALFMINW
                    b%urx = x + HALFMINW
                end if
            else if (b%llx + MINW > b%urx) then
                x = (b%llx + b%urx)/2
                b%llx = x - HALFMINW
                b%urx = x + HALFMINW
            end if
        end associate
    end do
    do i = 0, p%nbox - 2
        associate (bp1 => p%boxes(i + 1), bp2 => p%boxes(i + 2))
            if (i >= fb .and. i <= lb .and. modulo(i - fb, 2) == 0) then
                if (bp1%llx + MINW > bp2%urx) bp2%urx = bp1%llx + MINW
                if (bp1%urx - MINW < bp2%llx) bp2%llx = bp1%urx - MINW
            else if (i + 1 >= fb .and. i < lb .and. modulo(i + 1 - fb, 2) == 0) then
                if (bp1%llx + MINW > bp2%urx) bp1%llx = bp2%urx - MINW
                if (bp1%urx - MINW < bp2%llx) bp1%urx = bp2%llx + MINW
            end if
        end associate
    end do
end subroutine adjust_regular_path

!> The box of the space below rank `r`, kept for the next edge.
subroutine rank_box(c, sp, r, b)
    type(core_t), intent(in) :: c               !! the state
    type(spline_info_t), intent(inout) :: sp    !! what is known about the drawing
    integer, intent(in) :: r                    !! the rank
    type(box_t), intent(out) :: b               !! the box

    b = sp%rank_box(r)
    if (b%llx == b%urx) then
        b%llx = sp%left_bound
        b%lly = c%nd(c%rk(r + 1)%v(0))%y + c%rk(r + 1)%ht2
        b%urx = sp%right_bound
        b%ury = c%nd(c%rk(r)%v(0))%y - c%rk(r)%ht1
        sp%rank_box(r) = b
    end if
end subroutine rank_box

!> How many edges from `n` go straight down.
integer function straight_len(c, n) result(cnt)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: n        !! the node
    integer :: v

    cnt = 0
    v = n
    do while (c%nd(v)%out%n > 0)
        v = c%ed(c%nd(v)%out%a(1))%head
        if (c%nd(v)%kind /= VIRTUAL) exit
        if (c%nd(v)%out%n /= 1 .or. c%nd(v)%in%n /= 1) exit
        if (c%nd(v)%x /= c%nd(n)%x) exit
        cnt = cnt + 1
    end do
end function straight_len

!> Skip a straight part of a chain, and return the edge after it.
integer function straight_path(c, e, cnt, plist) result(f)
    type(core_t), intent(in) :: c           !! the state
    integer, intent(in) :: e                !! the first edge
    integer, intent(in) :: cnt              !! how many to skip
    type(plist_t), intent(inout) :: plist   !! the points so far
    integer :: k

    f = e
    do k = 1, cnt
        f = c%nd(c%ed(f)%head)%out%a(1)
    end do
    call add_point(plist, plist%a(plist%n))
    call add_point(plist, plist%a(plist%n))
end function straight_path

!> Give the room that a path did not use back to the virtual nodes along it.
subroutine recover_slack(c, e, p)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the first edge of the path
    type(path_t), intent(in) :: p       !! the path
    integer :: vn, b

    b = 1
    vn = c%ed(e)%head
    do while (c%nd(vn)%kind == VIRTUAL .and. .not. spline_merge(c, vn))
        do while (b <= p%nbox)
            if (p%boxes(b)%lly <= c%nd(vn)%y) exit
            b = b + 1
        end do
        if (b > p%nbox) exit
        if (p%boxes(b)%ury >= c%nd(vn)%y) then
            if (c%nd(vn)%label /= 0) then
                call resize_vn(c%nd(vn), p%boxes(b)%llx, p%boxes(b)%urx, p%boxes(b)%urx + c%nd(vn)%rw)
            else
                call resize_vn(c%nd(vn), p%boxes(b)%llx, (p%boxes(b)%llx + p%boxes(b)%urx)/2, p%boxes(b)%urx)
            end if
        end if
        if (c%nd(vn)%out%n == 0) exit
        vn = c%ed(c%nd(vn)%out%a(1))%head
    end do
end subroutine recover_slack

!> Set the position and width of a virtual node.
subroutine resize_vn(vn, lx, cx, rx)
    type(vnode_t), intent(inout) :: vn  !! the node
    real(wp), intent(in) :: lx          !! the left side
    real(wp), intent(in) :: cx          !! the centre
    real(wp), intent(in) :: rx          !! the right side

    vn%x = cx
    vn%lw = cx - lx
    vn%rw = rx - cx
end subroutine resize_vn

!> The box for a path to go through at node `vn`: all the space up to the neighbours.
function maximal_bbox(c, sp, vn, ie, oe) result(rv)
    type(core_t), intent(in) :: c           !! the state
    type(spline_info_t), intent(in) :: sp   !! what is known about the drawing
    integer, intent(in) :: vn               !! the node
    integer, intent(in) :: ie               !! the edge in, or 0
    integer, intent(in) :: oe               !! the edge out, or 0
    type(box_t) :: rv
    real(wp) :: b, nb
    integer :: left, right

    b = c%nd(vn)%x - c%nd(vn)%lw - FUDGE
    left = neighbor(c, vn, ie, oe, -1)
    if (left /= 0) then
        nb = c%nd(left)%x + c%nd(left)%saved_rw
        if (c%nd(left)%kind == NORMAL) then
            nb = nb + NODESEP/2.0_wp
        else
            nb = nb + sp%splinesep
        end if
        if (nb < b) b = nb
        rv%llx = c_round(b)
    else
        rv%llx = min(c_round(b), sp%left_bound)
    end if

    if (c%nd(vn)%kind == VIRTUAL .and. c%nd(vn)%label /= 0) then
        b = c%nd(vn)%x + 10
    else
        b = c%nd(vn)%x + c%nd(vn)%rw + FUDGE
    end if
    right = neighbor(c, vn, ie, oe, 1)
    if (right /= 0) then
        nb = c%nd(right)%x - c%nd(right)%lw
        if (c%nd(right)%kind == NORMAL) then
            nb = nb - NODESEP/2.0_wp
        else
            nb = nb - sp%splinesep
        end if
        if (nb > b) b = nb
        rv%urx = c_round(b)
    else
        rv%urx = max(c_round(b), sp%right_bound)
    end if

    if (c%nd(vn)%kind == VIRTUAL .and. c%nd(vn)%label /= 0) then
        rv%urx = rv%urx - c%nd(vn)%rw
        if (rv%urx < rv%llx) rv%urx = c%nd(vn)%x
    end if
    rv%lly = c%nd(vn)%y - c%rk(c%nd(vn)%rank)%ht1
    rv%ury = c%nd(vn)%y + c%rk(c%nd(vn)%rank)%ht2
end function maximal_bbox

!> `round` of the C library: halves go away from zero.
elemental function c_round(x)
    real(wp), intent(in) :: x   !! the value
    real(wp) :: c_round

    c_round = anint(x)
end function c_round

!> The nearest node of the rank that a path at `vn` must stay clear of.
integer function neighbor(c, vn, ie, oe, dir) result(rv)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: vn       !! the node
    integer, intent(in) :: ie       !! the edge in, or 0
    integer, intent(in) :: oe       !! the edge out, or 0
    integer, intent(in) :: dir      !! -1 for left, 1 for right
    integer :: i, n, r

    rv = 0
    r = c%nd(vn)%rank
    i = c%nd(vn)%order + dir
    do while (i >= 0 .and. i < c%rk(r)%n)
        n = c%rk(r)%v(i)
        if (c%nd(n)%kind == VIRTUAL .and. c%nd(n)%label /= 0) then
            rv = n
            return
        end if
        if (c%nd(n)%kind == NORMAL) then
            rv = n
            return
        end if
        if (.not. pathscross(c, n, vn, ie, oe)) then
            rv = n
            return
        end if
        i = i + dir
    end do
end function neighbor

!> Whether the paths through `n0` and `n1` cross.
logical function pathscross(c, n0, n1, ie1, oe1)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: n0       !! one node
    integer, intent(in) :: n1       !! the other node
    integer, intent(in) :: ie1      !! the edge in of `n1`, or 0
    integer, intent(in) :: oe1      !! the edge out of `n1`, or 0
    integer :: e0, e1, na, nb, cnt
    logical :: order

    pathscross = .false.
    order = c%nd(n0)%order > c%nd(n1)%order
    if (c%nd(n0)%out%n /= 1 .and. c%nd(n1)%out%n /= 1) return
    e1 = oe1
    if (c%nd(n0)%out%n == 1 .and. e1 /= 0) then
        e0 = c%nd(n0)%out%a(1)
        do cnt = 1, 2
            na = c%ed(e0)%head
            nb = c%ed(e1)%head
            if (na == nb) exit
            if (order .neqv. (c%nd(na)%order > c%nd(nb)%order)) then
                pathscross = .true.
                return
            end if
            if (c%nd(na)%out%n /= 1 .or. c%nd(na)%kind == NORMAL) exit
            e0 = c%nd(na)%out%a(1)
            if (c%nd(nb)%out%n /= 1 .or. c%nd(nb)%kind == NORMAL) exit
            e1 = c%nd(nb)%out%a(1)
        end do
    end if
    e1 = ie1
    if (c%nd(n0)%in%n == 1 .and. e1 /= 0) then
        e0 = c%nd(n0)%in%a(1)
        do cnt = 1, 2
            na = c%ed(e0)%tail
            nb = c%ed(e1)%tail
            if (na == nb) exit
            if (order .neqv. (c%nd(na)%order > c%nd(nb)%order)) then
                pathscross = .true.
                return
            end if
            if (c%nd(na)%in%n /= 1 .or. c%nd(na)%kind == NORMAL) exit
            e0 = c%nd(na)%in%a(1)
            if (c%nd(nb)%in%n /= 1 .or. c%nd(nb)%kind == NORMAL) exit
            e1 = c%nd(nb)%in%a(1)
        end do
    end if
end function pathscross

!> Route a path through its boxes. This turns the boxes into a polygon, finds a shortest path inside
!> it and fits a spline. The boxes are then shrunk to what the spline uses, for later edges.
subroutine routesplines(pp, ps, stat)
    type(path_t), intent(inout) :: pp   !! the path
    type(plist_t), intent(out) :: ps    !! the control points
    integer, intent(out) :: stat        !! the result
    type(point_t), allocatable :: polypoints(:), ba(:), bb(:)
    type(plist_t) :: pl, spl
    type(point_t) :: eps(2), ev0, ev1
    real(wp) :: v, delta
    logical :: flip, unbounded, is_horizontal, is_vertical, is_bad
    integer :: bi, boxn, pi, prev, next, i, loopcnt

    call checkpath(pp, stat)
    if (stat /= 0) return
    boxn = pp%nbox
    allocate (polypoints(boxn*8))
    flip = .false.
    if (boxn > 1) flip = pp%boxes(1)%lly > pp%boxes(2)%lly
    if (flip) then
        do bi = 1, boxn
            v = pp%boxes(bi)%ury
            pp%boxes(bi)%ury = -pp%boxes(bi)%lly
            pp%boxes(bi)%lly = -v
        end do
    end if

    pi = 0
    do bi = 1, boxn
        next = 0
        prev = 0
        if (bi > 1) prev = merge(-1, 1, pp%boxes(bi)%lly > pp%boxes(bi - 1)%lly)
        if (bi < boxn) next = merge(1, -1, pp%boxes(bi + 1)%lly > pp%boxes(bi)%lly)
        if (prev /= next) then
            if (next == -1 .or. prev == 1) then
                call left_side(pp%boxes(bi))
            else
                call right_side(pp%boxes(bi))
            end if
        else if (prev == 0) then
            call left_side(pp%boxes(bi))
        else if (.not. (prev == -1 .and. next == -1)) then
            stat = 4
            return
        end if
    end do
    do bi = boxn, 1, -1
        next = 0
        prev = 0
        if (bi < boxn) prev = merge(-1, 1, pp%boxes(bi)%lly > pp%boxes(bi + 1)%lly)
        if (bi > 1) next = merge(1, -1, pp%boxes(bi - 1)%lly > pp%boxes(bi)%lly)
        if (prev /= next) then
            if (next == -1 .or. prev == 1) then
                call left_side(pp%boxes(bi))
            else
                call right_side(pp%boxes(bi))
            end if
        else if (prev == 0) then
            call right_side(pp%boxes(bi))
        else
            if (.not. (prev == -1 .and. next == -1)) then
                stat = 4
                return
            end if
            call right_side(pp%boxes(bi))
            call left_side(pp%boxes(bi))
        end if
    end do

    if (flip) then
        do bi = 1, boxn
            v = pp%boxes(bi)%ury
            pp%boxes(bi)%ury = -pp%boxes(bi)%lly
            pp%boxes(bi)%lly = -v
        end do
        polypoints(1:pi)%y = -polypoints(1:pi)%y
    end if

    pp%boxes(1:boxn)%llx = huge(v)
    pp%boxes(1:boxn)%urx = -huge(v)
    eps = [pp%start%p, pp%end%p]
    call shortest_path(polypoints(1:pi), eps, pl, stat)
    if (stat /= 0) return

    allocate (ba(pi), bb(pi))
    do i = 1, pi
        ba(i) = polypoints(i)
        bb(i) = polypoints(mod(i, pi) + 1)
    end do
    ev0 = point_t(0, 0)
    ev1 = point_t(0, 0)
    if (pp%start%constrained) ev0 = point_t(cos(pp%start%theta), sin(pp%start%theta))
    if (pp%end%constrained) ev1 = point_t(-cos(pp%end%theta), -sin(pp%end%theta))
    call route_spline(ba, bb, pl%a(1:pl%n), ev0, ev1, spl)

    ps = spl
    unbounded = .true.
    is_horizontal = all(abs(spl%a(1)%y - spl%a(1:spl%n)%y) <= SMALL_GAP)
    is_vertical = all(abs(spl%a(1)%x - spl%a(1:spl%n)%x) <= SMALL_GAP)
    if (is_horizontal .or. is_vertical) then
        pp%boxes(1:boxn)%llx = spl%a(1)%x
        pp%boxes(1:boxn)%urx = spl%a(1)%x
        unbounded = .false.
    end if
    delta = INIT_DELTA
    do loopcnt = 1, LOOP_TRIES
        if (.not. unbounded) exit
        call limit_boxes(pp%boxes(1:boxn), spl%a(1:spl%n), delta)
        is_bad = any(pp%boxes(1:boxn)%llx == huge(v)) .or. any(pp%boxes(1:boxn)%urx == -huge(v))
        if (is_bad) then
            delta = delta*2
        else
            unbounded = .false.
        end if
    end do
    if (unbounded) call limit_boxes(pp%boxes(1:boxn), polyline_of(pl), real(INIT_DELTA, wp))

contains

    !> Add the left side of a box to the polygon.
    subroutine left_side(b)
        type(box_t), intent(in) :: b    !! the box

        polypoints(pi + 1) = point_t(b%llx, b%ury)
        polypoints(pi + 2) = point_t(b%llx, b%lly)
        pi = pi + 2
    end subroutine left_side

    !> Add the right side of a box to the polygon.
    subroutine right_side(b)
        type(box_t), intent(in) :: b    !! the box

        polypoints(pi + 1) = point_t(b%urx, b%lly)
        polypoints(pi + 2) = point_t(b%urx, b%ury)
        pi = pi + 2
    end subroutine right_side
end subroutine routesplines

!> The control points of a curve that follows a polyline.
function polyline_of(line) result(points)
    type(plist_t), intent(in) :: line   !! the polyline
    type(point_t), allocatable :: points(:)
    integer :: i, n

    allocate (points(3*line%n - 2))
    n = 0
    points(1:2) = line%a(1)
    n = 2
    do i = 2, line%n - 1
        points(n + 1:n + 3) = line%a(i)
        n = n + 3
    end do
    points(n + 1:n + 2) = line%a(line%n)
end function polyline_of

!> Shrink the boxes to what the curve uses.
subroutine limit_boxes(boxes, pps, delta)
    type(box_t), intent(inout) :: boxes(:)  !! the boxes
    type(point_t), intent(in) :: pps(:)     !! the curve
    real(wp), intent(in) :: delta           !! how finely to look
    type(point_t) :: sp(0:3)
    real(wp) :: t, num_div
    integer :: splinepi, si, bi, steps

    num_div = delta*size(boxes)
    steps = int(num_div)
    do splinepi = 1, size(pps) - 3, 3
        do si = 0, steps
            t = si/num_div
            sp = pps(splinepi:splinepi + 3)
            sp(0)%x = sp(0)%x + t*(sp(1)%x - sp(0)%x)
            sp(0)%y = sp(0)%y + t*(sp(1)%y - sp(0)%y)
            sp(1)%x = sp(1)%x + t*(sp(2)%x - sp(1)%x)
            sp(1)%y = sp(1)%y + t*(sp(2)%y - sp(1)%y)
            sp(2)%x = sp(2)%x + t*(sp(3)%x - sp(2)%x)
            sp(2)%y = sp(2)%y + t*(sp(3)%y - sp(2)%y)
            sp(0)%x = sp(0)%x + t*(sp(1)%x - sp(0)%x)
            sp(0)%y = sp(0)%y + t*(sp(1)%y - sp(0)%y)
            sp(1)%x = sp(1)%x + t*(sp(2)%x - sp(1)%x)
            sp(1)%y = sp(1)%y + t*(sp(2)%y - sp(1)%y)
            sp(0)%x = sp(0)%x + t*(sp(1)%x - sp(0)%x)
            sp(0)%y = sp(0)%y + t*(sp(1)%y - sp(0)%y)
            do bi = 1, size(boxes)
                if (sp(0)%y <= boxes(bi)%ury + SMALL_GAP .and. sp(0)%y >= boxes(bi)%lly - SMALL_GAP) then
                    boxes(bi)%llx = min(boxes(bi)%llx, sp(0)%x)
                    boxes(bi)%urx = max(boxes(bi)%urx, sp(0)%x)
                end if
            end do
        end do
    end do
end subroutine limit_boxes

!> The length of the overlap of two intervals.
function overlap(i0, i1, j0, j1) result(r)
    real(wp), intent(in) :: i0  !! the start of one
    real(wp), intent(in) :: i1  !! the end of one
    real(wp), intent(in) :: j0  !! the start of the other
    real(wp), intent(in) :: j1  !! the end of the other
    real(wp) :: r

    if (i1 <= j0 .or. i0 >= j1) then
        r = 0
    else if (i0 <= j0 .and. i1 >= j1) then
        r = i1 - i0
    else if (j0 <= i0 .and. j1 >= i1) then
        r = j1 - j0
    else if (j0 <= i0 .and. i0 <= j1) then
        r = j1 - i0
    else
        r = i1 - j0
    end if
end function overlap

!> Repair minor errors in the boxes, such as boxes that do not join or that overlap slightly. Unlike
!> the original, the count of boxes follows the removal of degenerate ones.
subroutine checkpath(p, stat)
    type(path_t), intent(inout) :: p    !! the path
    integer, intent(out) :: stat        !! the result
    integer :: i, bi, boxn, l, r, d, u, errs, j
    real(wp) :: xy, xoverlap, yoverlap

    stat = 0
    i = 0
    do bi = 1, p%nbox
        if (abs(p%boxes(bi)%lly - p%boxes(bi)%ury) < 0.01_wp) cycle
        if (abs(p%boxes(bi)%llx - p%boxes(bi)%urx) < 0.01_wp) cycle
        i = i + 1
        p%boxes(i) = p%boxes(bi)
    end do
    boxn = i
    p%nbox = boxn
    if (boxn == 0) then
        stat = 5
        return
    end if
    if (p%boxes(1)%llx > p%boxes(1)%urx .or. p%boxes(1)%lly > p%boxes(1)%ury) then
        stat = 5
        return
    end if
    do bi = 1, boxn - 1
        associate (ba => p%boxes(bi), bb => p%boxes(bi + 1))
            if (bb%llx > bb%urx .or. bb%lly > bb%ury) then
                stat = 5
                return
            end if
            l = merge(1, 0, ba%urx < bb%llx)
            r = merge(1, 0, ba%llx > bb%urx)
            d = merge(1, 0, ba%ury < bb%lly)
            u = merge(1, 0, ba%lly > bb%ury)
            errs = l + r + d + u
            if (errs > 0) then
                if (l == 1) then
                    xy = ba%urx
                    ba%urx = bb%llx
                    bb%llx = xy
                    l = 0
                else if (r == 1) then
                    xy = ba%llx
                    ba%llx = bb%urx
                    bb%urx = xy
                    r = 0
                else if (d == 1) then
                    xy = ba%ury
                    ba%ury = bb%lly
                    bb%lly = xy
                    d = 0
                else if (u == 1) then
                    xy = ba%lly
                    ba%lly = bb%ury
                    bb%ury = xy
                    u = 0
                end if
                do j = 1, errs - 1
                    if (l == 1) then
                        xy = (ba%urx + bb%llx)/2 + 0.5_wp
                        ba%urx = xy
                        bb%llx = xy
                        l = 0
                    else if (r == 1) then
                        xy = (ba%llx + bb%urx)/2 + 0.5_wp
                        ba%llx = xy
                        bb%urx = xy
                        r = 0
                    else if (d == 1) then
                        xy = (ba%ury + bb%lly)/2 + 0.5_wp
                        ba%ury = xy
                        bb%lly = xy
                        d = 0
                    else if (u == 1) then
                        xy = (ba%lly + bb%ury)/2 + 0.5_wp
                        ba%lly = xy
                        bb%ury = xy
                        u = 0
                    end if
                end do
            end if
            xoverlap = overlap(ba%llx, ba%urx, bb%llx, bb%urx)
            yoverlap = overlap(ba%lly, ba%ury, bb%lly, bb%ury)
            if (xoverlap > 0 .and. yoverlap > 0) then
                if (xoverlap < yoverlap) then
                    if (ba%urx - ba%llx > bb%urx - bb%llx) then
                        if (ba%urx < bb%urx) then
                            ba%urx = bb%llx
                        else
                            ba%llx = bb%urx
                        end if
                    else if (ba%urx < bb%urx) then
                        bb%llx = ba%urx
                    else
                        bb%urx = ba%llx
                    end if
                else if (ba%ury - ba%lly > bb%ury - bb%lly) then
                    if (ba%ury < bb%ury) then
                        ba%ury = bb%lly
                    else
                        ba%lly = bb%ury
                    end if
                else if (ba%ury < bb%ury) then
                    bb%lly = ba%ury
                else
                    bb%ury = ba%lly
                end if
            end if
        end associate
    end do
    p%start%p%x = min(max(p%start%p%x, p%boxes(1)%llx), p%boxes(1)%urx)
    p%start%p%y = min(max(p%start%p%y, p%boxes(1)%lly), p%boxes(1)%ury)
    p%end%p%x = min(max(p%end%p%x, p%boxes(boxn)%llx), p%boxes(boxn)%urx)
    p%end%p%y = min(max(p%end%p%y, p%boxes(boxn)%lly), p%boxes(boxn)%ury)
end subroutine checkpath

!> Draw self-loops to the right of their node.
subroutine self_right(c, edges, stepx, sizey)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: edges(:)     !! the edges
    real(wp), intent(in) :: stepx       !! the distance between loops
    real(wp), intent(in) :: sizey       !! the room above and below
    type(point_t) :: tp, hp, np, points(7)
    real(wp) :: stepy, dx, dy, hx, tx, width
    integer :: i, e, n

    n = c%ed(edges(1))%tail
    stepy = max(sizey/2/size(edges), 2.0_wp)
    np = point_t(c%nd(n)%x, c%nd(n)%y)
    tp = np
    hp = np
    dx = c%nd(n)%rw
    dy = 0
    tx = min(dx, 3*(np%x + dx - tp%x))
    hx = min(dx, 3*(np%x + dx - hp%x))
    do i = 1, size(edges)
        e = edges(i)
        dx = dx + stepx
        tx = tx + stepx
        hx = hx + stepx
        dy = dy + stepy
        points = [tp, point_t(tp%x + tx/3, tp%y + dy), point_t(np%x + dx, tp%y + dy), &
                  point_t(np%x + dx, (tp%y + hp%y)/2), point_t(np%x + dx, hp%y - dy), &
                  point_t(hp%x + hx/3, hp%y - dy), hp]
        if (c%labelled(e)) then
            width = merge(c%label_h(e), c%label_w(e), c%flip)
            c%label_pos(e) = point_t(c%nd(n)%x + dx + width/2, c%nd(n)%y)
            if (width > stepx) dx = dx + width - stepx
        end if
        call clip_and_install(c, e, c%ed(e)%head, points)
    end do
end subroutine self_right

!> Whether a point relative to node `n` is inside the node, as drawn.
logical function inside_node(c, n, p)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: n        !! the node
    type(point_t), intent(in) :: p  !! the point

    if (c%flip) then
        inside_node = abs(p%x) <= c%out_h(n) .and. abs(p%y) <= c%out_w(n)
    else
        inside_node = abs(p%x) <= c%out_w(n) .and. abs(p%y) <= c%out_h(n)
    end if
end function inside_node

!> A point on a curve, and optionally the two parts on each side of it.
subroutine bezier(v, t, pt, left, right)
    type(point_t), intent(in) :: v(4)                   !! the curve
    real(wp), intent(in) :: t                           !! the parameter
    type(point_t), intent(out) :: pt                    !! the point
    type(point_t), intent(out), optional :: left(4)     !! the part before
    type(point_t), intent(out), optional :: right(4)    !! the part after
    type(point_t) :: w(0:3, 0:3)
    integer :: i, j

    w(0, :) = v
    do i = 1, 3
        do j = 0, 3 - i
            w(i, j)%x = (1 - t)*w(i - 1, j)%x + t*w(i - 1, j + 1)%x
            w(i, j)%y = (1 - t)*w(i - 1, j)%y + t*w(i - 1, j + 1)%y
        end do
    end do
    if (present(left)) then
        do j = 0, 3
            left(j + 1) = w(j, 0)
        end do
    end if
    if (present(right)) then
        do j = 0, 3
            right(j + 1) = w(3 - j, j)
        end do
    end if
    pt = w(3, 0)
end subroutine bezier

!> Clip a curve to a region by a binary search. `node` is the node of the region, or 0 for a circle
!> around `centre`. `left_inside` says that the first point is inside.
subroutine bezier_clip(c, node, centre, radius, sp, left_inside)
    type(core_t), intent(in) :: c       !! the state
    integer, intent(in) :: node         !! the node, in whose coordinates the curve is, or 0
    type(point_t), intent(in) :: centre !! the centre of the circle
    real(wp), intent(in) :: radius      !! its radius
    type(point_t), intent(inout) :: sp(4)   !! the curve
    logical, intent(in) :: left_inside  !! which end is inside
    type(point_t) :: seg(4), best(4), pt, opt
    real(wp) :: low, high, t
    logical :: found, is_in

    pt = merge(sp(1), sp(4), left_inside)
    found = .false.
    low = 0
    high = 1
    do
        opt = pt
        t = (high + low)/2
        if (left_inside) then
            call bezier(sp, t, pt, right=seg)
        else
            call bezier(sp, t, pt, left=seg)
        end if
        if (node /= 0) then
            is_in = inside_node(c, node, pt)
        else
            is_in = hypot(pt%x - centre%x, pt%y - centre%y) <= radius
        end if
        if (is_in) then
            if (left_inside) then
                low = t
            else
                high = t
            end if
            best = seg
            found = .true.
        else if (left_inside) then
            high = t
        else
            low = t
        end if
        if (.not. (abs(opt%x - pt%x) > 0.5_wp .or. abs(opt%y - pt%y) > 0.5_wp)) exit
    end do
    if (found) then
        sp = best
    else
        sp = seg
    end if
end subroutine bezier_clip

!> Clip a curve to the outline of node `n`.
subroutine shape_clip(c, n, curve, left_inside)
    type(core_t), intent(in) :: c           !! the state
    integer, intent(in) :: n                !! the node
    type(point_t), intent(inout) :: curve(4)    !! the curve
    logical, intent(in) :: left_inside      !! which end is inside
    type(point_t) :: q(4)

    q = curve
    q%x = q%x - c%nd(n)%x
    q%y = q%y - c%nd(n)%y
    call bezier_clip(c, n, point_t(0, 0), 0.0_wp, q, left_inside)
    q%x = q%x + c%nd(n)%x
    q%y = q%y + c%nd(n)%y
    curve = q
end subroutine shape_clip

!> How long an arrowhead is: its tip stands back so that the stroke touches the node.
function arrow_length(c) result(length)
    type(core_t), intent(in) :: c   !! the state
    real(wp) :: length, nominal, tip

    if (c%arrow_size == 0) then
        length = 0
        return
    end if
    nominal = c%arrow_size*ARROW_NOMINAL
    tip = EDGE_PEN/2/sin(atan(0.35_wp))
    length = nominal + tip
end function arrow_length

!> Whether two points are the same, nearly.
logical function approx_eq(p, q)
    type(point_t), intent(in) :: p  !! one point
    type(point_t), intent(in) :: q  !! the other point

    approx_eq = (p%x - q%x)**2 + (p%y - q%y)**2 < MILLIPOINT**2
end function approx_eq

!> Clip a raw curve to the nodes at its ends, make room for arrowheads, and attach it to its edge.
subroutine clip_and_install(c, fe, hn_in, ps)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: fe           !! the edge
    integer, intent(in) :: hn_in        !! the node where it ends
    type(point_t), intent(inout) :: ps(0:)  !! the control points
    type(spline_t) :: newspl
    type(point_t) :: q(4)
    integer :: tn, hn, orig, e, start, end, pn, i
    logical :: is_in, sflag, eflag, s, j

    tn = c%ed(fe)%tail
    hn = hn_in
    pn = size(ps)
    orig = fe
    do while (c%ed(orig)%to_orig /= 0 .and. c%ed(orig)%etype /= NORMAL)
        orig = c%ed(orig)%to_orig
    end do

    start = 0
    if (c%nd(tn)%kind == NORMAL) then
        do while (start < pn - 4)
            is_in = inside_node(c, tn, point_t(ps(start + 3)%x - c%nd(tn)%x, ps(start + 3)%y - c%nd(tn)%y))
            if (.not. is_in) exit
            start = start + 3
        end do
        q = ps(start:start + 3)
        call shape_clip(c, tn, q, .true.)
        ps(start:start + 3) = q
    end if
    end = pn - 4
    if (c%nd(hn)%kind == NORMAL) then
        do while (end > 0)
            is_in = inside_node(c, hn, point_t(ps(end)%x - c%nd(hn)%x, ps(end)%y - c%nd(hn)%y))
            if (.not. is_in) exit
            end = end - 3
        end do
        q = ps(end:end + 3)
        call shape_clip(c, hn, q, .false.)
        ps(end:end + 3) = q
    end if
    do while (start < pn - 4)
        if (.not. approx_eq(ps(start), ps(start + 3))) exit
        start = start + 3
    end do
    do while (end > 0)
        if (.not. approx_eq(ps(end), ps(end + 3))) exit
        end = end - 3
    end do

    ! arrowheads
    e = fe
    do while (c%ed(e)%to_orig /= 0)
        e = c%ed(e)%to_orig
    end do
    j = swap_ends(c, e)
    sflag = c%ed(e)%conc_opp
    eflag = .true.
    if (spline_merge(c, hn)) eflag = .false.
    if (spline_merge(c, c%ed(fe)%tail)) sflag = .false.
    if (j) then
        s = sflag
        sflag = eflag
        eflag = s
    end if
    if (sflag) call arrow_start_clip(c, ps, start, end, newspl)
    if (eflag) call arrow_end_clip(c, ps, start, end, newspl)

    newspl%list = ps(start:end + 3)
    i = start
    do while (i + 3 <= end + 3)
        call update_bb_bz(c%bb, ps(i:i + 3))
        i = i + 3
    end do
    call add_spline(c%spl(orig), newspl)
end subroutine clip_and_install

!> Attach a curve to an edge.
subroutine add_spline(path, spline)
    type(edge_path_t), intent(inout) :: path    !! the path
    type(spline_t), intent(in) :: spline        !! the curve
    type(spline_t), allocatable :: bigger(:)
    integer :: k

    allocate (bigger(path%n + 1))
    do k = 1, path%n
        bigger(k) = path%b(k)
    end do
    bigger(path%n + 1) = spline
    call move_alloc(bigger, path%b)
    path%n = path%n + 1
end subroutine add_spline

!> Make room for an arrowhead at the end of a curve.
subroutine arrow_end_clip(c, ps, startp, endp, spl)
    type(core_t), intent(in) :: c       !! the state
    type(point_t), intent(inout) :: ps(0:)  !! the control points
    integer, intent(in) :: startp       !! the first point in use
    integer, intent(inout) :: endp      !! the last curve in use
    type(spline_t), intent(inout) :: spl    !! the result
    type(point_t) :: sp(4)
    real(wp) :: elen

    elen = arrow_length(c)
    spl%eflag = .true.
    spl%ep = ps(endp + 3)
    if (endp > startp) then
        if (hypot(ps(endp)%x - ps(endp + 3)%x, ps(endp)%y - ps(endp + 3)%y) < elen) endp = endp - 3
    end if
    sp = [spl%ep, ps(endp + 2), ps(endp + 1), ps(endp)]
    if (elen > 0) call bezier_clip(c, 0, sp(1), elen, sp, .true.)
    ps(endp) = sp(4)
    ps(endp + 1) = sp(3)
    ps(endp + 2) = sp(2)
    ps(endp + 3) = sp(1)
end subroutine arrow_end_clip

!> Make room for an arrowhead at the start of a curve.
subroutine arrow_start_clip(c, ps, startp, endp, spl)
    type(core_t), intent(in) :: c       !! the state
    type(point_t), intent(inout) :: ps(0:)  !! the control points
    integer, intent(inout) :: startp    !! the first point in use
    integer, intent(in) :: endp         !! the last curve in use
    type(spline_t), intent(inout) :: spl    !! the result
    type(point_t) :: sp(4)
    real(wp) :: slen

    slen = arrow_length(c)
    spl%sflag = .true.
    spl%sp = ps(startp)
    if (endp > startp) then
        if (hypot(ps(startp)%x - ps(startp + 3)%x, ps(startp)%y - ps(startp + 3)%y) < slen) startp = startp + 3
    end if
    sp = [ps(startp + 3), ps(startp + 2), ps(startp + 1), spl%sp]
    if (slen > 0) call bezier_clip(c, 0, sp(4), slen, sp, .false.)
    ps(startp) = sp(4)
    ps(startp + 1) = sp(3)
    ps(startp + 2) = sp(2)
    ps(startp + 3) = sp(1)
end subroutine arrow_start_clip

!> The square of the distance from point `p` to the line through `a` and `b`.
function pt_to_line2(a, b, p) result(r)
    type(point_t), intent(in) :: a  !! a point of the line
    type(point_t), intent(in) :: b  !! another point of the line
    type(point_t), intent(in) :: p  !! the point
    real(wp) :: r, dx, dy, a2

    dx = b%x - a%x
    dy = b%y - a%y
    a2 = (p%y - a%y)*dx - (p%x - a%x)*dy
    a2 = a2*a2
    if (a2 < 1.0e-7_wp) then
        r = 0
    else
        r = a2/(dx*dx + dy*dy)
    end if
end function pt_to_line2

!> Grow a box to hold a curve, refining the curve until it is nearly straight.
recursive subroutine update_bb_bz(bb, cp)
    type(box_t), intent(inout) :: bb    !! the box
    type(point_t), intent(in) :: cp(4)  !! the curve
    type(point_t) :: left(4), right(4), pt

    if (any(cp%x > bb%urx) .or. any(cp%x < bb%llx) .or. any(cp%y > bb%ury) .or. any(cp%y < bb%lly)) then
        if (pt_to_line2(cp(1), cp(4), cp(2)) < HW*HW .and. pt_to_line2(cp(1), cp(4), cp(3)) < HW*HW) then
            bb%urx = max(bb%urx, maxval(cp%x))
            bb%llx = min(bb%llx, minval(cp%x))
            bb%ury = max(bb%ury, maxval(cp%y))
            bb%lly = min(bb%lly, minval(cp%y))
        else
            call bezier(cp, 0.5_wp, pt, left, right)
            call update_bb_bz(bb, left)
            call update_bb_bz(bb, right)
        end if
    end if
end subroutine update_bb_bz

!> Reverse a path, so that it goes from tail to head.
subroutine swap_spline(s)
    type(edge_path_t), intent(inout) :: s   !! the path
    type(spline_t) :: t
    type(point_t) :: q
    logical :: f
    integer :: i

    do i = 1, s%n/2
        t = s%b(i)
        s%b(i) = s%b(s%n + 1 - i)
        s%b(s%n + 1 - i) = t
    end do
    do i = 1, s%n
        s%b(i)%list = s%b(i)%list(size(s%b(i)%list):1:-1)
        f = s%b(i)%sflag
        s%b(i)%sflag = s%b(i)%eflag
        s%b(i)%eflag = f
        q = s%b(i)%sp
        s%b(i)%sp = s%b(i)%ep
        s%b(i)%ep = q
    end do
end subroutine swap_spline

end module fdot_route
