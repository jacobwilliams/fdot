module fdot_structure
!! The state of a layout, kept in integer-indexed arrays and passed through the calls, and the basic
!! operations on the working graph. This is a port of `lib/dotgen/fastgr.c`, `decomp.c` and parts of
!! `class2.c` of Graphviz, which are under the Eclipse Public License 2.0. Graphviz keeps this state in
!! global variables and in records attached to graph objects, which cannot be safe in threads.
use iso_fortran_env, only: real64, int64
implicit none
private

public :: wp, ilist_t, vnode_t, vedge_t, rank_t, core_t, point_t, box_t, spline_t, edge_path_t
public :: push, remove, sort_ints
public :: new_node, new_edge, fast_edge, delete_fast_edge, find_fast_edge
public :: new_virtual_edge, virtual_edge, virtual_node, merge_oneway, merge_chain, other_edge
public :: decompose, round, main_edge

integer, parameter :: wp = real64   !! the kind of reals

integer, parameter, public :: NORMAL = 0        !! a node or edge of the drawing
integer, parameter, public :: VIRTUAL = 1       !! a node or edge added for the layout
integer, parameter, public :: IGNORED = 6       !! an edge merged away by concentration

integer, parameter, public :: NODESEP = 18      !! points between nodes of a rank
integer, parameter, public :: RANKSEP = 36      !! points between ranks

type :: point_t
    !! A point.
    real(wp) :: x = 0   !! x
    real(wp) :: y = 0   !! y
end type point_t

type :: box_t
    !! A rectangle.
    real(wp) :: llx = 0     !! the left side
    real(wp) :: lly = 0     !! the bottom side
    real(wp) :: urx = 0     !! the right side
    real(wp) :: ury = 0     !! the top side
end type box_t

type :: spline_t
    !! A piecewise cubic Bezier curve, with its arrowheads, in the order from tail to head.
    type(point_t), allocatable :: list(:)   !! the control points
    logical :: sflag = .false.              !! whether there is an arrowhead at the start
    logical :: eflag = .false.              !! whether there is an arrowhead at the end
    type(point_t) :: sp                     !! where the arrowhead at the start ends
    type(point_t) :: ep                     !! where the arrowhead at the end ends
end type spline_t

type :: edge_path_t
    !! The curves of an edge. Concentrated edges share parts, so an edge may have several or none.
    integer :: n = 0                    !! how many there are
    type(spline_t), allocatable :: b(:) !! the curves
end type edge_path_t

type :: ilist_t
    !! A list of integers that grows.
    integer :: n = 0                !! how many there are
    integer, allocatable :: a(:)    !! the values
end type ilist_t

type :: vnode_t
    !! A node of the working graph.
    integer :: kind = NORMAL        !! `NORMAL` or `VIRTUAL`
    integer :: rank = 0             !! the rank
    integer :: order = 0            !! the position in the rank
    integer :: save_order = 0       !! the order of the best result so far
    integer :: weight_class = 0     !! for the weights of edges
    integer :: label = 0            !! the edge whose label this node holds, or 0
    real(wp) :: lw = 1              !! the width left of the centre
    real(wp) :: rw = 1              !! the width right of the centre
    real(wp) :: ht = 1              !! the height
    real(wp) :: mval = 0            !! the median, in crossing reduction
    real(wp) :: saved_rw = 0        !! `rw`, without the space for self-loops
    real(wp) :: x = 0               !! x
    real(wp) :: y = 0               !! y
    logical :: alive = .true.       !! false for a node merged away
    type(ilist_t) :: in             !! the edges in
    type(ilist_t) :: out            !! the edges out
    type(ilist_t) :: other          !! the self-loops
end type vnode_t

type :: vedge_t
    !! An edge of the working graph.
    integer :: tail = 0             !! where it starts
    integer :: head = 0             !! where it ends
    integer :: count = 1            !! how many edges it stands for
    integer :: xpenalty = 1         !! the cost of a crossing
    integer :: weight = 1           !! the weight
    integer :: minlen = 1           !! the smallest number of ranks it spans
    integer :: to_virt = 0          !! the virtual edge that stands for this one
    integer :: to_orig = 0          !! the edge this one stands for
    integer :: etype = VIRTUAL      !! `NORMAL`, `VIRTUAL` or `IGNORED`
    integer :: flags = 0            !! how this edge is routed
    integer :: seq = 0              !! the edge of the drawing this one was made from
    logical :: conc_opp = .false.   !! whether an opposite edge is merged into this one
end type vedge_t

type :: rank_t
    !! The nodes of a rank.
    integer :: n = 0                    !! how many there are now
    integer :: an = 0                   !! how many there can be
    integer :: base = 0                 !! where the current component starts in `v`
    integer, allocatable :: v(:)        !! the nodes, ending with 0
    real(wp) :: ht1 = 0                 !! the height below
    real(wp) :: ht2 = 0                 !! the height above
    real(wp) :: pht1 = 0                !! the height below of the nodes themselves
    real(wp) :: pht2 = 0                !! the height above of the nodes themselves
    logical :: valid = .false.          !! whether `cache_nc` is valid
    logical :: candidate = .false.      !! whether a transposition may help
    integer(int64) :: cache_nc = 0      !! the crossings below this rank
end type rank_t

type :: core_t
    !! All the state of a layout.
    logical :: concentrate = .false.    !! whether to merge edges
    logical :: flip = .false.           !! true when ranks go across the page
    logical :: has_labels = .false.     !! whether some edge has a label
    integer :: nn = 0                   !! number of nodes of the drawing, the first in `nd`
    integer :: ne = 0                   !! number of edges of the drawing, the first in `ed`
    integer :: nnodes = 0               !! number of nodes in `nd`
    integer :: nedges = 0               !! number of edges in `ed`
    integer :: ranksep = RANKSEP        !! points between ranks
    integer :: minrank = 0              !! the first rank
    integer :: maxrank = 0              !! the last rank
    type(vnode_t), allocatable :: nd(:) !! the nodes
    type(vedge_t), allocatable :: ed(:) !! the edges
    type(rank_t), allocatable :: rk(:)  !! the ranks, from 0 to `maxrank + 1`
    type(ilist_t), allocatable :: comp(:)   !! the nodes of each connected component
    type(ilist_t), allocatable :: gout(:)   !! the edges out of each node, as `agfstout` gives them
    integer, allocatable :: nlist(:)    !! all the nodes, in the order of the components
    real(wp), allocatable :: label_w(:) !! the width of the label of each edge
    real(wp), allocatable :: label_h(:) !! the height of the label of each edge
    logical, allocatable :: labelled(:) !! whether each edge has a label
    type(point_t), allocatable :: label_pos(:)  !! where the labels of edges go
    type(edge_path_t), allocatable :: spl(:)    !! the path of each edge of the drawing
    real(wp), allocatable :: out_w(:)   !! half the width of each node's outline, as drawn
    real(wp), allocatable :: out_h(:)   !! half the height of each node's outline, as drawn
    type(box_t) :: bb                   !! the extent of the drawing
    real(wp) :: ht1 = 0                 !! the height of the last rank below its centre
    real(wp) :: ht2 = 0                 !! the height of the first rank above its centre
    real(wp) :: arrow_size = 0.7_wp     !! the size of arrowheads, as a part of the standard size
end type core_t

contains

!> Round like the C macro `ROUND`, to an integer.
elemental integer function round(x)
    real(wp), intent(in) :: x   !! the value

    if (x >= 0) then
        round = int(x + 0.5_wp)
    else
        round = int(x - 0.5_wp)
    end if
end function round

!> Append `value` to `list`.
pure subroutine push(list, value)
    type(ilist_t), intent(inout) :: list    !! the list
    integer, value :: value                 !! the value, a copy as it can be an element of `list`
    integer, allocatable :: bigger(:)

    if (.not. allocated(list%a)) then
        allocate (list%a(4))
    else if (list%n == size(list%a)) then
        allocate (bigger(2*list%n))
        bigger(1:list%n) = list%a(1:list%n)
        call move_alloc(bigger, list%a)
    end if
    list%n = list%n + 1
    list%a(list%n) = value
end subroutine push

!> Remove `value` from `list`, filling the hole with the last member.
pure subroutine remove(list, value)
    type(ilist_t), intent(inout) :: list    !! the list
    integer, intent(in) :: value            !! the value
    integer :: i

    do i = 1, list%n
        if (list%a(i) == value) then
            list%a(i) = list%a(list%n)
            list%n = list%n - 1
            exit
        end if
    end do
end subroutine remove

!> Sort the first `n` integers of `a` in increasing order.
pure subroutine sort_ints(a, n)
    integer, intent(inout) :: a(:)  !! the values
    integer, intent(in) :: n        !! how many there are
    integer :: gap, i, j, t

    gap = n/2
    do while (gap > 0)
        do i = gap + 1, n
            t = a(i)
            j = i
            do while (j > gap)
                if (a(j - gap) <= t) exit
                a(j) = a(j - gap)
                j = j - gap
            end do
            a(j) = t
        end do
        gap = gap/2
    end do
end subroutine sort_ints

!> Add a node, and return its number.
pure subroutine new_node(c, v)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: v           !! the node
    type(vnode_t), allocatable :: bigger(:)
    integer :: k

    if (.not. allocated(c%nd)) then
        allocate (c%nd(16))
    else if (c%nnodes == size(c%nd)) then
        allocate (bigger(2*c%nnodes))
        do k = 1, c%nnodes
            bigger(k)%kind = c%nd(k)%kind
            bigger(k)%rank = c%nd(k)%rank
            bigger(k)%order = c%nd(k)%order
            bigger(k)%save_order = c%nd(k)%save_order
            bigger(k)%weight_class = c%nd(k)%weight_class
            bigger(k)%label = c%nd(k)%label
            bigger(k)%lw = c%nd(k)%lw
            bigger(k)%rw = c%nd(k)%rw
            bigger(k)%ht = c%nd(k)%ht
            bigger(k)%mval = c%nd(k)%mval
            bigger(k)%saved_rw = c%nd(k)%saved_rw
            bigger(k)%x = c%nd(k)%x
            bigger(k)%y = c%nd(k)%y
            bigger(k)%alive = c%nd(k)%alive
            call take_list(c%nd(k)%in, bigger(k)%in)
            call take_list(c%nd(k)%out, bigger(k)%out)
            call take_list(c%nd(k)%other, bigger(k)%other)
        end do
        call move_alloc(bigger, c%nd)
    end if
    c%nnodes = c%nnodes + 1
    v = c%nnodes
end subroutine new_node

!> Move list `from` to `to`, without copying it.
pure subroutine take_list(from, to)
    type(ilist_t), intent(inout) :: from    !! the list
    type(ilist_t), intent(inout) :: to      !! where it goes

    to%n = from%n
    if (allocated(from%a)) call move_alloc(from%a, to%a)
end subroutine take_list

!> Add an edge that is in no list, and return its number. The ends are copies, as this can move the
!> edges, and callers often pass ends of other edges.
pure subroutine new_edge(c, u, v, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, value :: u                 !! where it starts
    integer, value :: v                 !! where it ends
    integer, intent(out) :: e           !! the edge
    type(vedge_t), allocatable :: bigger(:)

    if (.not. allocated(c%ed)) then
        allocate (c%ed(16))
    else if (c%nedges == size(c%ed)) then
        allocate (bigger(2*c%nedges))
        bigger(1:c%nedges) = c%ed(1:c%nedges)
        call move_alloc(bigger, c%ed)
    end if
    c%nedges = c%nedges + 1
    e = c%nedges
    c%ed(e) = vedge_t(tail=u, head=v)
end subroutine new_edge

!> Connect edge `e` to its nodes.
pure subroutine fast_edge(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge

    call push(c%nd(c%ed(e)%tail)%out, e)
    call push(c%nd(c%ed(e)%head)%in, e)
end subroutine fast_edge

!> Disconnect edge `e` from its nodes.
pure subroutine delete_fast_edge(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge

    call remove(c%nd(c%ed(e)%tail)%out, e)
    call remove(c%nd(c%ed(e)%head)%in, e)
end subroutine delete_fast_edge

!> The edge from `u` to `v`, or 0.
pure integer function find_fast_edge(c, u, v) result(e)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: u        !! where it starts
    integer, intent(in) :: v        !! where it ends
    integer :: i

    e = 0
    if (c%nd(u)%out%n > 0 .and. c%nd(v)%in%n > 0) then
        if (c%nd(u)%out%n < c%nd(v)%in%n) then
            do i = 1, c%nd(u)%out%n
                if (c%ed(c%nd(u)%out%a(i))%head == v) then
                    e = c%nd(u)%out%a(i)
                    return
                end if
            end do
        else
            do i = 1, c%nd(v)%in%n
                if (c%ed(c%nd(v)%in%a(i))%tail == u) then
                    e = c%nd(v)%in%a(i)
                    return
                end if
            end do
        end if
    end if
end function find_fast_edge

!> Create a virtual edge attached to `orig`, which may be 0, and return it.
pure subroutine new_virtual_edge(c, u, v, orig, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, value :: u                 !! where it starts
    integer, value :: v                 !! where it ends
    integer, value :: orig              !! the edge it stands for, or 0
    integer, intent(out) :: e           !! the result

    call new_edge(c, u, v, e)
    if (orig /= 0) then
        c%ed(e)%count = c%ed(orig)%count
        c%ed(e)%xpenalty = c%ed(orig)%xpenalty
        c%ed(e)%weight = c%ed(orig)%weight
        c%ed(e)%minlen = c%ed(orig)%minlen
        if (c%ed(orig)%to_virt == 0) c%ed(orig)%to_virt = e
        c%ed(e)%to_orig = orig
        c%ed(e)%seq = c%ed(orig)%seq
    end if
end subroutine new_virtual_edge

!> Create a virtual edge attached to `orig`, and connect it.
pure subroutine virtual_edge(c, u, v, orig, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, value :: u                 !! where it starts
    integer, value :: v                 !! where it ends
    integer, value :: orig              !! the edge it stands for, or 0
    integer, intent(out) :: e           !! the result

    call new_virtual_edge(c, u, v, orig, e)
    call fast_edge(c, e)
end subroutine virtual_edge

!> Create a virtual node.
pure subroutine virtual_node(c, v)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(out) :: v           !! the result

    call new_node(c, v)
    c%nd(v)%kind = VIRTUAL
end subroutine virtual_node

!> Make `rep` stand for `e` as well.
pure subroutine merge_oneway(c, e, rep)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    integer, intent(in) :: rep          !! the edge that stands for it
    integer :: r

    if (rep == c%ed(e)%to_virt .or. e == c%ed(rep)%to_virt) return
    c%ed(e)%to_virt = rep
    if (c%ed(rep)%minlen < c%ed(e)%minlen) c%ed(rep)%minlen = c%ed(e)%minlen
    r = rep
    do while (r /= 0)
        c%ed(r)%count = c%ed(r)%count + c%ed(e)%count
        c%ed(r)%xpenalty = c%ed(r)%xpenalty + c%ed(e)%xpenalty
        c%ed(r)%weight = c%ed(r)%weight + c%ed(e)%weight
        r = c%ed(r)%to_virt
    end do
end subroutine merge_oneway

!> Make the chain of virtual edges from `f` stand for `e` as well.
pure subroutine merge_chain(c, e, f, update_count)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge
    integer, intent(in) :: f            !! the first edge of the chain
    logical, intent(in) :: update_count !! whether to add to the counts of edges
    integer :: rep, lastrank, width

    lastrank = max(c%nd(c%ed(e)%tail)%rank, c%nd(c%ed(e)%head)%rank)
    c%ed(e)%to_virt = f
    rep = f
    width = NODESEP/2
    do
        if (update_count) c%ed(rep)%count = c%ed(rep)%count + c%ed(e)%count
        c%ed(rep)%xpenalty = c%ed(rep)%xpenalty + c%ed(e)%xpenalty
        c%ed(rep)%weight = c%ed(rep)%weight + c%ed(e)%weight
        if (c%nd(c%ed(rep)%head)%rank == lastrank) exit
        associate (v => c%nd(c%ed(rep)%head))
            v%lw = v%lw + width
            v%rw = v%rw + width
            if (v%out%n == 0) exit
            rep = v%out%a(1)
        end associate
    end do
end subroutine merge_chain

!> Note that `e` is routed apart from the graph of ranks.
pure subroutine other_edge(c, e)
    type(core_t), intent(inout) :: c    !! the state
    integer, intent(in) :: e            !! the edge

    call push(c%nd(c%ed(e)%tail)%other, e)
end subroutine other_edge

!> The edge of the drawing that `e` stands for.
pure integer function main_edge(c, e) result(le)
    type(core_t), intent(in) :: c   !! the state
    integer, intent(in) :: e        !! the edge

    le = e
    do while (c%ed(le)%to_virt /= 0)
        le = c%ed(le)%to_virt
    end do
    do while (c%ed(le)%to_orig /= 0)
        le = c%ed(le)%to_orig
    end do
end function main_edge

!> Find the connected components of the working graph, starting from the nodes of the drawing.
pure subroutine decompose(c)
    type(core_t), intent(inout) :: c    !! the state
    logical, allocatable :: seen(:)
    type(ilist_t) :: stack
    type(ilist_t), allocatable :: comps(:), bigger(:)
    integer :: start, n, i, e, other, ncomp

    allocate (seen(c%nnodes), source=.false.)
    allocate (comps(max(4, c%nn)))
    ncomp = 0
    do start = 1, c%nn
        if (seen(start)) cycle
        ncomp = ncomp + 1
        if (ncomp > size(comps)) then
            allocate (bigger(2*size(comps)))
            do i = 1, ncomp - 1
                call take_list(comps(i), bigger(i))
            end do
            call move_alloc(bigger, comps)
        end if
        call push(stack, start)
        do while (stack%n > 0)
            n = stack%a(stack%n)
            stack%n = stack%n - 1
            if (seen(n)) cycle
            seen(n) = .true.
            call push(comps(ncomp), n)
            do i = c%nd(n)%in%n, 1, -1
                e = c%nd(n)%in%a(i)
                other = c%ed(e)%head
                if (other == n) other = c%ed(e)%tail
                if (.not. seen(other)) call push(stack, other)
            end do
            do i = c%nd(n)%out%n, 1, -1
                e = c%nd(n)%out%a(i)
                other = c%ed(e)%head
                if (other == n) other = c%ed(e)%tail
                if (.not. seen(other)) call push(stack, other)
            end do
        end do
    end do
    if (allocated(c%comp)) deallocate (c%comp)
    allocate (c%comp(ncomp))
    do i = 1, ncomp
        call take_list(comps(i), c%comp(i))
    end do
end subroutine decompose

end module fdot_structure
