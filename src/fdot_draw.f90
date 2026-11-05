module fdot_draw
!! Node sizes, arrowheads and the final placement of a drawing. This is a port of parts of
!! `lib/common/shapes.c`, `arrows.c` and `postproc.c` of Graphviz, which are under the Eclipse Public
!! License 2.0, for rounded boxes and `normal` arrows.
use fdot_structure
implicit none
private

public :: node_size, arrow_polygon, translate_drawing

real(wp), parameter :: MIN_WIDTH = 54.0_wp      !! the default smallest width of a node: 0.75 inches
real(wp), parameter :: ARROW_LENGTH = 10.0_wp   !! the length of an arrowhead of size 1
real(wp), parameter :: EPSILON = 0.0001_wp      !! keeps directions stable at small lengths
real(wp), parameter :: EDGE_PEN = 1.0_wp        !! the width of the lines of edges

contains

!> The size of a box around text of width `text_w` and height `text_h`, with margins given in points.
pure subroutine node_size(text_w, text_h, margin_x, margin_y, height, width, node_height)
    real(wp), intent(in) :: text_w          !! the width of the text
    real(wp), intent(in) :: text_h          !! the height of the text
    real(wp), intent(in) :: margin_x        !! the margin on each side
    real(wp), intent(in) :: margin_y        !! the margin above and below
    real(wp), intent(in) :: height          !! the smallest height
    real(wp), intent(out) :: width          !! the width
    real(wp), intent(out) :: node_height    !! the height

    width = MIN_WIDTH
    if (text_w > 0 .or. text_h > 0) then
        width = max(width, text_w + 2*margin_x)
        node_height = max(height, text_h + 2*margin_y)
    else
        node_height = height
    end if
end subroutine node_size

!> The corners of the arrowhead that ends at `p`, on the line from `towards`.
pure subroutine arrow_polygon(p, towards, arrow_size, corners)
    type(point_t), intent(in) :: p          !! where the arrowhead ends
    type(point_t), intent(in) :: towards    !! the point it points away from
    real(wp), intent(in) :: arrow_size      !! its size
    type(point_t), intent(out) :: corners(3)    !! the corners
    real(wp), parameter :: ARROW_WIDTH = 0.35_wp   !! half the width of an arrowhead, as a part of its length
    real(wp) :: ux, uy, s, vx, vy, len, tip, dx, dy

    ux = towards%x - p%x
    uy = towards%y - p%y
    s = ARROW_LENGTH/(hypot(ux, uy) + EPSILON)
    ux = (ux + sign(EPSILON, ux))*s*arrow_size
    uy = (uy + sign(EPSILON, uy))*s*arrow_size
    vx = -uy*ARROW_WIDTH
    vy = ux*ARROW_WIDTH
    len = hypot(ux, uy)
    ! the tip stands back, so that the corner of the stroke touches the node
    tip = EDGE_PEN/2/sin(atan(ARROW_WIDTH))
    dx = 0
    dy = 0
    if (len > 0) then
        dx = tip*ux/len
        dy = tip*uy/len
    end if
    corners(1) = point_t(p%x + dx + ux - vx, p%y + dy + uy - vy)
    corners(2) = point_t(p%x + dx, p%y + dy)
    corners(3) = point_t(p%x + dx + ux + vx, p%y + dy + uy + vy)
end subroutine arrow_polygon

!> Turn the drawing for its direction, and move it so that the corner of its extent is at the origin.
pure subroutine translate_drawing(c, right_to_left)
    type(core_t), intent(inout) :: c        !! the state
    logical, intent(in) :: right_to_left    !! true if ranks go from right to left
    type(point_t) :: a, b, offset
    integer :: n, e, j, k

    a = map(point_t(c%bb%llx, c%bb%lly))
    b = map(point_t(c%bb%urx, c%bb%ury))
    offset = point_t(min(a%x, b%x), min(a%y, b%y))
    do n = 1, c%nn
        a = shift(point_t(c%nd(n)%x, c%nd(n)%y))
        c%nd(n)%x = a%x
        c%nd(n)%y = a%y
    end do
    do e = 1, c%ne
        do j = 1, c%spl(e)%n
            associate (s => c%spl(e)%b(j))
                do k = 1, size(s%list)
                    s%list(k) = shift(s%list(k))
                end do
                if (s%sflag) s%sp = shift(s%sp)
                if (s%eflag) s%ep = shift(s%ep)
            end associate
        end do
        if (c%labelled(e)) c%label_pos(e) = shift(c%label_pos(e))
    end do
    a = map(point_t(c%bb%llx, c%bb%lly))
    b = map(point_t(c%bb%urx, c%bb%ury))
    c%bb = box_t(0, 0, abs(a%x - b%x), abs(a%y - b%y))

contains

    !> Turn a point.
    pure function map(p) result(r)
        type(point_t), intent(in) :: p  !! the point
        type(point_t) :: r

        if (right_to_left) then
            r = point_t(p%y, p%x)
        else
            r = point_t(-p%y, p%x)
        end if
    end function map

    !> Turn and move a point.
    pure function shift(p) result(r)
        type(point_t), intent(in) :: p  !! the point
        type(point_t) :: r

        r = map(p)
        r%x = r%x - offset%x
        r%y = r%y - offset%y
    end function shift
end subroutine translate_drawing

end module fdot_draw
