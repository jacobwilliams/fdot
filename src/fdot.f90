module fdot
!! The graph layout of Graphviz `dot`, in modern Fortran: a port of the part that draws directed graphs
!! of named boxes, for fdoc. Give it nodes and edges, and get an SVG drawing back as a string. There is
!! no DOT language, no plugins, no files and no external process.
!!
!! Graphviz is under the Eclipse Public License 2.0, and so is this port. It is safe to use different
!! graphs in different threads at once.
use iso_fortran_env, only: int64
use fdot_structure
use fdot_text, only: text_width, text_height
use fdot_rank, only: dot_rank
use fdot_order, only: dot_mincross
use fdot_position, only: dot_position
use fdot_route, only: dot_splines
use fdot_draw, only: node_size, arrow_polygon, translate_drawing
use fdot_svg
implicit none
private

public :: graph_t, node_style_t, edge_style_t, point_t, wp

integer, parameter, public :: RANK_LR = 1   !! ranks go from left to right
integer, parameter, public :: RANK_RL = 2   !! ranks go from right to left

integer, parameter, public :: FDOT_OK = 0               !! success
integer, parameter, public :: FDOT_INVALID = 1          !! a node of an edge does not exist
integer, parameter, public :: FDOT_NO_RANKS = 2         !! network simplex failed
integer, parameter, public :: FDOT_FLAT_EDGE = 3        !! an edge joins nodes of one rank
integer, parameter, public :: FDOT_NO_ROUTE = 4         !! no route was found for an edge

real(wp), parameter :: PAD = 4.0_wp                     !! points around a drawing
real(wp), parameter :: BASELINE_SHARE = 0.36_wp               !! where the middle of text is above its baseline, in font sizes

type :: node_style_t
    !! How nodes look, by default as fdoc draws them.
    real(wp) :: font_size = 13.0_wp     !! the size of names, in points
    real(wp) :: height = 21.6_wp        !! the smallest height, in points
    real(wp) :: margin_x = 7.2_wp       !! the margin beside names, in points
    real(wp) :: margin_y = 2.88_wp      !! the margin above and below names, in points
    real(wp) :: pen_width = 2.0_wp      !! the width of an outline, in points
end type node_style_t

type :: edge_style_t
    !! How edges look, by default as fdoc draws them.
    real(wp) :: font_size = 11.0_wp                 !! the size of labels, in points
    real(wp) :: arrow_size = 0.7_wp                 !! the size of arrowheads
    character(len=:), allocatable :: colour         !! the colour of edges that give none
    character(len=:), allocatable :: label_colour   !! the colour of labels
end type edge_style_t

type :: node_data_t
    !! What was given for a node.
    character(len=:), allocatable :: name, fill, link, outline
    character(len=:), allocatable :: owner, owner_fill, owner_link  !! the left part of a node of two parts
end type node_data_t

type :: edge_data_t
    !! What was given for an edge.
    integer :: from = 0, to = 0
    logical :: dashed = .false.
    character(len=:), allocatable :: label, colour
end type edge_data_t

type :: graph_t
    !! A graph, and its layout.
    private
    integer :: direction = RANK_LR
    logical :: concentrate = .false.
    logical :: done = .false.
    type(node_style_t) :: node_style
    type(edge_style_t) :: edge_style
    integer :: nn = 0, ne = 0
    type(node_data_t), allocatable :: nodes(:)
    type(edge_data_t), allocatable :: edges(:)
    type(core_t) :: c
contains
    private
    procedure, public :: initialize         !! start a new graph
    procedure, public :: add_node           !! add a node
    procedure, public :: add_edge           !! add an edge
    procedure, public :: layout             !! compute the layout
    procedure, public :: write_svg          !! draw the layout
    procedure, public :: write_dot          !! write the graph in the DOT language, for debugging
    procedure, public :: node_count         !! how many nodes there are
    procedure, public :: edge_count         !! how many edges there are
    procedure, public :: size_of_drawing    !! the size of the layout
    procedure, public :: node_position      !! where a node is
    procedure, public :: edge_curves        !! how many curves an edge has
    procedure, public :: edge_curve         !! the control points of a curve of an edge
    procedure, public :: edge_label_position    !! where the label of an edge is
end type graph_t

contains

!> Start a new graph. By default ranks go from left to right, and edges are not concentrated.
pure subroutine initialize(this, direction, concentrate, node_style, edge_style)
    class(graph_t), intent(inout) :: this                   !! the graph
    integer, intent(in), optional :: direction              !! `RANK_LR` or `RANK_RL`
    logical, intent(in), optional :: concentrate            !! whether to merge edges that run together
    type(node_style_t), intent(in), optional :: node_style  !! how nodes look
    type(edge_style_t), intent(in), optional :: edge_style  !! how edges look

    this%direction = RANK_LR
    if (present(direction)) this%direction = direction
    this%concentrate = .false.
    if (present(concentrate)) this%concentrate = concentrate
    this%node_style = node_style_t()
    if (present(node_style)) this%node_style = node_style
    this%edge_style = edge_style_t()
    if (present(edge_style)) this%edge_style = edge_style
    if (.not. allocated(this%edge_style%colour)) this%edge_style%colour = '#8a93a5'
    if (.not. allocated(this%edge_style%label_colour)) this%edge_style%label_colour = '#8a93a5'
    this%nn = 0
    this%ne = 0
    this%done = .false.
    if (allocated(this%nodes)) deallocate (this%nodes)
    if (allocated(this%edges)) deallocate (this%edges)
    allocate (this%nodes(16), this%edges(16))
end subroutine initialize

!> Add a node, and return its number. The colours and the links are passed through as text. With
!> `owner`, the node has two parts joined side by side: the owner on the left, with a colour and a link
!> of its own, and the name on the right. This is for something that belongs to something else.
pure subroutine add_node(this, name, id, fill, link, outline, owner, owner_fill, owner_link)
    class(graph_t), intent(inout) :: this           !! the graph
    character(len=*), intent(in) :: name            !! the text of the node
    integer, intent(out) :: id                      !! its number
    character(len=*), intent(in), optional :: fill  !! the fill colour
    character(len=*), intent(in), optional :: link  !! where it links to
    character(len=*), intent(in), optional :: outline   !! the colour of an outline
    character(len=*), intent(in), optional :: owner     !! the text of the left part
    character(len=*), intent(in), optional :: owner_fill    !! the fill colour of the left part
    character(len=*), intent(in), optional :: owner_link    !! where the left part links to
    type(node_data_t), allocatable :: bigger(:)
    integer :: k

    if (.not. allocated(this%nodes)) call this%initialize()
    if (this%nn == size(this%nodes)) then
        allocate (bigger(2*this%nn))
        do k = 1, this%nn
            call move_alloc(this%nodes(k)%name, bigger(k)%name)
            call move_alloc(this%nodes(k)%fill, bigger(k)%fill)
            call move_alloc(this%nodes(k)%link, bigger(k)%link)
            call move_alloc(this%nodes(k)%outline, bigger(k)%outline)
            call move_alloc(this%nodes(k)%owner, bigger(k)%owner)
            call move_alloc(this%nodes(k)%owner_fill, bigger(k)%owner_fill)
            call move_alloc(this%nodes(k)%owner_link, bigger(k)%owner_link)
        end do
        call move_alloc(bigger, this%nodes)
    end if
    this%nn = this%nn + 1
    id = this%nn
    this%done = .false.
    associate (n => this%nodes(id))
        n%name = name
        n%fill = '#646d7e'
        if (present(fill)) n%fill = fill
        n%link = ''
        if (present(link)) n%link = link
        n%outline = ''
        if (present(outline)) n%outline = outline
        n%owner = ''
        if (present(owner)) n%owner = owner
        n%owner_fill = '#646d7e'
        if (present(owner_fill)) n%owner_fill = owner_fill
        n%owner_link = ''
        if (present(owner_link)) n%owner_link = owner_link
    end associate
end subroutine add_node

!> Add an edge between two nodes. An edge with a label is dashed unless told otherwise.
pure subroutine add_edge(this, from, to, colour, label, dashed)
    class(graph_t), intent(inout) :: this               !! the graph
    integer, intent(in) :: from                         !! where it starts
    integer, intent(in) :: to                           !! where it ends
    character(len=*), intent(in), optional :: colour    !! the colour
    character(len=*), intent(in), optional :: label     !! the label
    logical, intent(in), optional :: dashed             !! whether it is dashed
    type(edge_data_t), allocatable :: bigger(:)
    integer :: k

    if (.not. allocated(this%edges)) call this%initialize()
    if (this%ne == size(this%edges)) then
        allocate (bigger(2*this%ne))
        do k = 1, this%ne
            bigger(k)%from = this%edges(k)%from
            bigger(k)%to = this%edges(k)%to
            bigger(k)%dashed = this%edges(k)%dashed
            call move_alloc(this%edges(k)%label, bigger(k)%label)
            call move_alloc(this%edges(k)%colour, bigger(k)%colour)
        end do
        call move_alloc(bigger, this%edges)
    end if
    this%ne = this%ne + 1
    this%done = .false.
    associate (e => this%edges(this%ne))
        e%from = from
        e%to = to
        e%colour = this%edge_style%colour
        if (present(colour)) e%colour = colour
        e%label = ''
        if (present(label)) e%label = label
        e%dashed = len(e%label) > 0
        if (present(dashed)) e%dashed = dashed
    end associate
end subroutine add_edge

!> Compute the layout. `stat` is `FDOT_OK`, or says what went wrong.
subroutine layout(this, stat)
    class(graph_t), intent(inout) :: this   !! the graph
    integer, intent(out) :: stat            !! the result
    type(core_t) :: fresh
    integer, allocatable :: count(:)
    integer :: k, v, e, h, r

    stat = FDOT_OK
    this%done = .false.
    this%c = fresh
    if (this%nn == 0) then
        this%done = .true.
        return
    end if
    do e = 1, this%ne
        if (min(this%edges(e)%from, this%edges(e)%to) < 1 .or. max(this%edges(e)%from, this%edges(e)%to) > this%nn) then
            stat = FDOT_INVALID
            return
        end if
    end do
    associate (c => this%c)
        c%concentrate = this%concentrate
        c%flip = .true.
        c%arrow_size = this%edge_style%arrow_size
        call set_up_nodes(this)
        call set_up_edges(this)

        ! the edges of each node in the order `agfstout` gives them: by the node they go to
        allocate (c%gout(this%nn), count(this%nn))
        count = 0
        do e = 1, this%ne
            count(this%edges(e)%to) = count(this%edges(e)%to) + 1
        end do
        k = 1
        do v = 1, this%nn
            h = count(v)
            count(v) = k
            k = k + h
        end do
        block
            integer, allocatable :: sorted(:)

            allocate (sorted(this%ne))
            do e = 1, this%ne
                sorted(count(this%edges(e)%to)) = e
                count(this%edges(e)%to) = count(this%edges(e)%to) + 1
            end do
            do k = 1, this%ne
                call push(c%gout(this%edges(sorted(k))%from), sorted(k))
            end do
        end block

        call dot_rank(c, r)
        if (r /= 0) then
            stat = FDOT_NO_RANKS
            return
        end if
        call dot_mincross(c, r)
        if (r /= 0) then
            stat = FDOT_FLAT_EDGE
            return
        end if
        call dot_position(c, r)
        if (r /= 0) then
            stat = FDOT_NO_RANKS
            return
        end if
        call dot_splines(c, r)
        if (r /= 0) then
            stat = FDOT_NO_ROUTE
            return
        end if
        call translate_drawing(c, this%direction == RANK_RL)
    end associate
    this%done = .true.
end subroutine layout

!> Make the nodes of the working graph, with their sizes.
pure subroutine set_up_nodes(this)
    class(graph_t), intent(inout) :: this   !! the graph
    real(wp) :: w, h, pen, text_w
    integer :: k, v

    allocate (this%c%out_w(this%nn), this%c%out_h(this%nn))
    do k = 1, this%nn
        call new_node(this%c, v)
        associate (n => this%c%nd(v), s => this%node_style)
            text_w = text_width(this%nodes(k)%name, s%font_size)
            if (len(this%nodes(k)%owner) > 0) text_w = text_w + owner_width(this, k)
            call node_size(text_w, text_height(s%font_size), s%margin_x, s%margin_y, s%height, w, h)
            ! ranks go across the page, so the layout works with the box on its side
            n%lw = h/2
            n%rw = h/2
            n%ht = w
            pen = merge(s%pen_width, 1.0_wp, len(this%nodes(k)%outline) > 0)
            this%c%out_w(v) = (w + pen)/2
            this%c%out_h(v) = (h + pen)/2
        end associate
    end do
    this%c%nn = this%nn
end subroutine set_up_nodes

!> The width of the left part of node `k`.
pure real(wp) function owner_width(this, k)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: k            !! the node

    owner_width = text_width(this%nodes(k)%owner, this%node_style%font_size) + 2*this%node_style%margin_x
end function owner_width

!> Make the edges of the working graph, with the sizes of their labels.
pure subroutine set_up_edges(this)
    class(graph_t), intent(inout) :: this   !! the graph
    integer :: k, e

    allocate (this%c%label_w(this%ne), this%c%label_h(this%ne), this%c%labelled(this%ne), this%c%label_pos(this%ne))
    do k = 1, this%ne
        call new_edge(this%c, this%edges(k)%from, this%edges(k)%to, e)
        this%c%ed(e)%etype = NORMAL
        this%c%ed(e)%seq = e
        this%c%labelled(e) = len(this%edges(k)%label) > 0
        this%c%label_w(e) = 0
        this%c%label_h(e) = 0
        if (this%c%labelled(e)) then
            this%c%label_w(e) = text_width(this%edges(k)%label, this%edge_style%font_size)
            this%c%label_h(e) = text_height(this%edge_style%font_size)
            this%c%has_labels = .true.
        end if
    end do
    this%c%ne = this%ne
end subroutine set_up_edges

!> Draw the layout as SVG. The drawing is empty if there is no layout.
pure subroutine write_svg(this, svg)
    class(graph_t), intent(in) :: this                  !! the graph
    character(len=:), allocatable, intent(out) :: svg   !! the drawing
    type(buffer_t) :: b
    logical, allocatable :: drawn(:)
    integer :: v, k, e

    svg = ''
    if (.not. this%done .or. this%nn == 0) return
    associate (c => this%c)
        call svg_begin(b, c%bb%urx, c%bb%ury, PAD, this%edge_style%colour, this%node_style%font_size)
        allocate (drawn(this%nn), source=.false.)
        do v = 1, this%nn
            if (.not. drawn(v)) call draw_node(this, v, b)
            drawn(v) = .true.
            do k = 1, c%gout(v)%n
                e = c%gout(v)%a(k)
                if (.not. drawn(c%ed(e)%head)) call draw_node(this, c%ed(e)%head, b)
                drawn(c%ed(e)%head) = .true.
                call draw_edge(this, e, b)
            end do
        end do
        call svg_end(b)
    end associate
    call take(b, svg)
end subroutine write_svg

!> Draw a node.
pure subroutine draw_node(this, v, b)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: v            !! the node
    type(buffer_t), intent(inout) :: b  !! the drawing
    real(wp) :: w, h, x, y

    call node_position(this, v, x, y, w, h)
    associate (n => this%nodes(v))
        if (len(n%owner) > 0) then
            call draw_two_parts(this, v, x, y, w, h, b)
        else if (len(n%link) > 0) then
            call add(b, '<a href="')
            call add_xml(b, n%link)
            call add(b, '" fill="white" stroke="none">')
            call svg_box_in_link(b, x - w/2, -y - h/2, w, h, n%fill, n%outline, this%node_style%pen_width)
            call svg_text(b, x, -y + BASELINE_SHARE*this%node_style%font_size, n%name, 0.0_wp, '')
            call add(b, '</a>')
        else
            call svg_box(b, x - w/2, -y - h/2, w, h, n%fill, n%outline, this%node_style%pen_width)
            call svg_text(b, x, -y + BASELINE_SHARE*this%node_style%font_size, n%name, 0.0_wp, 'white')
        end if
    end associate
end subroutine draw_node

!> Draw a node of two parts: the whole box in the colour of the right part, the left part over it, and
!> last the outline, so that the left part does not cover it.
pure subroutine draw_two_parts(this, v, x, y, w, h, b)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: v            !! the node
    real(wp), intent(in) :: x           !! the centre
    real(wp), intent(in) :: y           !! the centre
    real(wp), intent(in) :: w           !! the width
    real(wp), intent(in) :: h           !! the height
    type(buffer_t), intent(inout) :: b  !! the drawing
    real(wp) :: left, top, split, baseline

    left = x - w/2
    top = -y - h/2
    split = owner_width(this, v)
    baseline = -y + BASELINE_SHARE*this%node_style%font_size
    associate (n => this%nodes(v))
        call open_link(n%link, b)
        call svg_box(b, left, top, w, h, n%fill, '', this%node_style%pen_width)
        call svg_text(b, left + (split + w)/2, baseline, n%name, 0.0_wp, 'white')
        call close_link(n%link, b)
        call open_link(n%owner_link, b)
        call svg_left_box(b, left, top, split, h, n%owner_fill)
        call svg_text(b, left + split/2, baseline, n%owner, 0.0_wp, 'white')
        call close_link(n%owner_link, b)
        if (len(n%outline) > 0) call svg_box(b, left, top, w, h, 'none', n%outline, this%node_style%pen_width)
    end associate
end subroutine draw_two_parts

!> Start a link, if there is one.
pure subroutine open_link(link, b)
    character(len=*), intent(in) :: link    !! where it goes
    type(buffer_t), intent(inout) :: b      !! the drawing

    if (len(link) == 0) return
    call add(b, '<a href="')
    call add_xml(b, link)
    call add(b, '" fill="white" stroke="none">')
end subroutine open_link

!> End a link, if there is one.
pure subroutine close_link(link, b)
    character(len=*), intent(in) :: link    !! where it goes
    type(buffer_t), intent(inout) :: b      !! the drawing

    if (len(link) > 0) call add(b, '</a>')
end subroutine close_link

!> Draw the box of a node inside a link, which has `stroke="none"` already.
pure subroutine svg_box_in_link(b, x, y, width, height, fill, outline, pen)
    type(buffer_t), intent(inout) :: b      !! the drawing
    real(wp), intent(in) :: x               !! the left side
    real(wp), intent(in) :: y               !! the top side, downwards
    real(wp), intent(in) :: width           !! the width
    real(wp), intent(in) :: height          !! the height
    character(len=*), intent(in) :: fill    !! the fill
    character(len=*), intent(in) :: outline !! the colour of the outline
    real(wp), intent(in) :: pen             !! the width of the outline
    type(buffer_t) :: box
    character(len=:), allocatable :: text
    integer :: at

    call svg_box(box, x, y, width, height, fill, outline, pen)
    call take(box, text)
    at = index(text, ' stroke="none"')
    if (at > 0) text = text(1:at - 1)//text(at + 14:)
    call add(b, text)
end subroutine svg_box_in_link

!> Draw an edge.
pure subroutine draw_edge(this, e, b)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: e            !! the edge
    type(buffer_t), intent(inout) :: b  !! the drawing
    type(point_t) :: corners(3)
    character(len=:), allocatable :: stroke
    integer :: j, n

    stroke = ''
    if (this%edges(e)%colour /= this%edge_style%colour) stroke = this%edges(e)%colour
    do j = 1, this%c%spl(e)%n
        associate (s => this%c%spl(e)%b(j))
            n = size(s%list)
            call svg_curve(b, s%list, stroke, this%edges(e)%dashed)
            if (s%sflag) then
                call arrow_polygon(s%sp, s%list(1), this%edge_style%arrow_size, corners)
                call svg_arrow(b, corners, this%edges(e)%colour, stroke)
            end if
            if (s%eflag) then
                call arrow_polygon(s%ep, s%list(n), this%edge_style%arrow_size, corners)
                call svg_arrow(b, corners, this%edges(e)%colour, stroke)
            end if
        end associate
    end do
    if (this%c%labelled(e)) then
        call svg_text(b, this%c%label_pos(e)%x, -this%c%label_pos(e)%y + BASELINE_SHARE*this%edge_style%font_size, &
                      this%edges(e)%label, this%edge_style%font_size, this%edge_style%label_colour)
    end if
end subroutine draw_edge

!> Write the graph in the DOT language, for debugging and for comparing with `dot`.
pure subroutine write_dot(this, text)
    class(graph_t), intent(in) :: this                  !! the graph
    character(len=:), allocatable, intent(out) :: text  !! the text
    character(len=*), parameter :: NL = new_line('a')
    type(buffer_t) :: b
    integer :: k

    call add(b, 'digraph g {'//NL//'rankdir='//merge('RL', 'LR', this%direction == RANK_RL)// &
             '; bgcolor="transparent"; concentrate='//merge('true ', 'false', this%concentrate)//';'//NL)
    call add(b, 'node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=')
    call add_tenths(b, this%node_style%font_size)
    call add(b, ', fontcolor=white, color=transparent, height=0.3, margin="0.1,0.04"];'//NL)
    call add(b, 'edge [color="'//this%edge_style%colour//'", arrowsize=')
    call add_tenths(b, this%edge_style%arrow_size)
    call add(b, '];'//NL)
    do k = 1, this%nn
        call add(b, 'n')
        call add_int(b, int(k, int64))
        call add(b, ' [label="')
        if (len(this%nodes(k)%owner) > 0) call add_dot_string(b, this%nodes(k)%owner//'%')
        call add_dot_string(b, this%nodes(k)%name)
        call add(b, '", fillcolor="'//this%nodes(k)%fill//'"')
        if (len(this%nodes(k)%outline) > 0) then
            call add(b, ', penwidth=')
            call add_tenths(b, this%node_style%pen_width)
            call add(b, ', color="'//this%nodes(k)%outline//'"')
        end if
        if (len(this%nodes(k)%link) > 0) then
            call add(b, ', URL="')
            call add_dot_string(b, this%nodes(k)%link)
            call add(b, '"')
        end if
        call add(b, '];'//NL)
    end do
    do k = 1, this%ne
        call add(b, 'n')
        call add_int(b, int(this%edges(k)%from, int64))
        call add(b, ' -> n')
        call add_int(b, int(this%edges(k)%to, int64))
        call add(b, ' [color="'//this%edges(k)%colour//'"')
        if (this%edges(k)%dashed) call add(b, ', style=dashed')
        if (len(this%edges(k)%label) > 0) then
            call add(b, ', fontname="Helvetica", fontsize=')
            call add_tenths(b, this%edge_style%font_size)
            call add(b, ', fontcolor="'//this%edge_style%label_colour//'", label="')
            call add_dot_string(b, this%edges(k)%label)
            call add(b, '"')
        end if
        call add(b, '];'//NL)
    end do
    call add(b, '}'//NL)
    call take(b, text)
end subroutine write_dot

!> Append text for a quoted string of DOT.
pure subroutine add_dot_string(b, text)
    type(buffer_t), intent(inout) :: b      !! the buffer
    character(len=*), intent(in) :: text    !! the text
    integer :: k

    do k = 1, len(text)
        if (text(k:k) == '"' .or. text(k:k) == '\') call add(b, '\')
        call add(b, text(k:k))
    end do
end subroutine add_dot_string

!> How many nodes there are.
pure integer function node_count(this)
    class(graph_t), intent(in) :: this  !! the graph

    node_count = this%nn
end function node_count

!> How many edges there are.
pure integer function edge_count(this)
    class(graph_t), intent(in) :: this  !! the graph

    edge_count = this%ne
end function edge_count

!> The size of the layout in points, without padding.
pure subroutine size_of_drawing(this, width, height)
    class(graph_t), intent(in) :: this  !! the graph
    real(wp), intent(out) :: width      !! the width
    real(wp), intent(out) :: height     !! the height

    width = this%c%bb%urx
    height = this%c%bb%ury
end subroutine size_of_drawing

!> Where a node is, in points, with `y` going up, and optionally its rank and place in the rank.
pure subroutine node_position(this, id, x, y, width, height, rank, order)
    class(graph_t), intent(in) :: this      !! the graph
    integer, intent(in) :: id               !! the node
    real(wp), intent(out) :: x              !! the centre
    real(wp), intent(out) :: y              !! the centre
    real(wp), intent(out) :: width          !! the width
    real(wp), intent(out) :: height         !! the height
    integer, intent(out), optional :: rank  !! the rank
    integer, intent(out), optional :: order !! the place in the rank

    x = this%c%nd(id)%x
    y = this%c%nd(id)%y
    ! the layout had the box on its side
    width = this%c%nd(id)%ht
    height = this%c%nd(id)%lw + this%c%nd(id)%rw
    if (present(rank)) rank = this%c%nd(id)%rank
    if (present(order)) order = this%c%nd(id)%order
end subroutine node_position

!> How many curves edge `id` has. A concentrated edge can have several, or none.
pure integer function edge_curves(this, id)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: id           !! the edge

    edge_curves = this%c%spl(id)%n
end function edge_curves

!> The control points of curve `k` of edge `id`, from tail to head, and the end points of arrowheads.
pure subroutine edge_curve(this, id, k, points, has_start, start, has_end, end)
    class(graph_t), intent(in) :: this      !! the graph
    integer, intent(in) :: id               !! the edge
    integer, intent(in) :: k                !! the curve
    type(point_t), allocatable, intent(out) :: points(:)    !! the control points
    logical, intent(out), optional :: has_start     !! whether there is an arrowhead at the start
    type(point_t), intent(out), optional :: start   !! where it ends
    logical, intent(out), optional :: has_end       !! whether there is an arrowhead at the end
    type(point_t), intent(out), optional :: end     !! where it ends

    points = this%c%spl(id)%b(k)%list
    if (present(has_start)) has_start = this%c%spl(id)%b(k)%sflag
    if (present(start)) start = this%c%spl(id)%b(k)%sp
    if (present(has_end)) has_end = this%c%spl(id)%b(k)%eflag
    if (present(end)) end = this%c%spl(id)%b(k)%ep
end subroutine edge_curve

!> Where the label of edge `id` is, if it has one.
pure subroutine edge_label_position(this, id, x, y, has_label)
    class(graph_t), intent(in) :: this  !! the graph
    integer, intent(in) :: id           !! the edge
    real(wp), intent(out) :: x          !! the centre
    real(wp), intent(out) :: y          !! the centre
    logical, intent(out) :: has_label   !! whether it has a label

    has_label = this%c%labelled(id)
    x = 0
    y = 0
    if (has_label) then
        x = this%c%label_pos(id)%x
        y = this%c%label_pos(id)%y
    end if
end subroutine edge_label_position

end module fdot
