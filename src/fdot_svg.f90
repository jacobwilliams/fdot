module fdot_svg
!! Writing text without internal `read` or `write` statements, which take a lock, and writing a drawing
!! as SVG. The output is the compact form that fdoc used to get by rewriting the output of `dot`, for
!! putting inline in a page: boxes are rectangles, what most elements share is on the group around them,
!! and coordinates keep one decimal, which also hides small differences in floating-point results
!! between platforms.
use iso_fortran_env, only: int64
use fdot_structure, only: wp, point_t
implicit none
private

public :: buffer_t, add, add_int, add_tenths, add_xml, take
public :: svg_begin, svg_end, svg_box, svg_text, svg_curve, svg_arrow

type :: buffer_t
    !! Text that grows.
    character(len=:), allocatable :: s  !! the room for the text
    integer :: n = 0                    !! how much is in use
end type buffer_t

contains

!> Append text.
pure subroutine add(b, text)
    type(buffer_t), intent(inout) :: b      !! the buffer
    character(len=*), intent(in) :: text    !! the text
    character(len=:), allocatable :: bigger

    if (.not. allocated(b%s)) then
        allocate (character(len=max(4096, 2*len(text))) :: b%s)
    else if (b%n + len(text) > len(b%s)) then
        allocate (character(len=2*(b%n + len(text))) :: bigger)
        bigger(1:b%n) = b%s(1:b%n)
        call move_alloc(bigger, b%s)
    end if
    b%s(b%n + 1:b%n + len(text)) = text
    b%n = b%n + len(text)
end subroutine add

!> Move the text to `text`, and empty the buffer.
pure subroutine take(b, text)
    type(buffer_t), intent(inout) :: b                  !! the buffer
    character(len=:), allocatable, intent(out) :: text  !! the text

    if (allocated(b%s)) then
        text = b%s(1:b%n)
    else
        text = ''
    end if
    b%n = 0
end subroutine take

!> Append an integer.
pure subroutine add_int(b, value)
    type(buffer_t), intent(inout) :: b  !! the buffer
    integer(int64), intent(in) :: value !! the value
    character(len=20) :: digits
    integer(int64) :: rest
    integer :: at

    rest = abs(value)
    at = len(digits)
    do
        digits(at:at) = achar(iachar('0') + int(mod(rest, 10_int64)))
        rest = rest/10
        if (rest == 0) exit
        at = at - 1
    end do
    if (value < 0) call add(b, '-')
    call add(b, digits(at:))
end subroutine add_int

!> Append a number with at most one decimal.
pure subroutine add_tenths(b, value)
    type(buffer_t), intent(inout) :: b  !! the buffer
    real(wp), intent(in) :: value       !! the value
    integer(int64) :: tenths

    tenths = nint(10*value, int64)
    if (tenths < 0) call add(b, '-')
    call add_int(b, abs(tenths)/10)
    if (mod(abs(tenths), 10_int64) /= 0) then
        call add(b, '.')
        call add_int(b, mod(abs(tenths), 10_int64))
    end if
end subroutine add_tenths

!> Append text for XML.
pure subroutine add_xml(b, text)
    type(buffer_t), intent(inout) :: b      !! the buffer
    character(len=*), intent(in) :: text    !! the text
    integer :: k

    do k = 1, len(text)
        select case (text(k:k))
        case ('&'); call add(b, '&amp;')
        case ('<'); call add(b, '&lt;')
        case ('>'); call add(b, '&gt;')
        case ('"'); call add(b, '&quot;')
        case default; call add(b, text(k:k))
        end select
    end do
end subroutine add_xml

!> Start a drawing of the given size, which is padded by `pad` on every side.
pure subroutine svg_begin(b, width, height, pad, edge_colour, font_size)
    type(buffer_t), intent(inout) :: b          !! the buffer
    real(wp), intent(in) :: width               !! the width
    real(wp), intent(in) :: height              !! the height
    real(wp), intent(in) :: pad                 !! the padding
    character(len=*), intent(in) :: edge_colour !! the colour of most edges
    real(wp), intent(in) :: font_size           !! the size of most text
    integer(int64) :: w, h

    w = nint(width + 2*pad, int64)
    h = nint(height + 2*pad, int64)
    call add(b, '<svg width="')
    call add_int(b, w)
    call add(b, 'pt" height="')
    call add_int(b, h)
    call add(b, 'pt"'//new_line('a')//' viewBox="0 0 ')
    call add_int(b, w)
    call add(b, ' ')
    call add_int(b, h)
    call add(b, '" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">')
    call add(b, '<g transform="scale(1 1) rotate(0) translate(')
    call add_tenths(b, pad)
    call add(b, ' ')
    call add_tenths(b, height + pad)
    call add(b, ')" fill="none" stroke="'//edge_colour//'" font-family="Helvetica,sans-Serif" font-size="')
    call add_tenths(b, font_size)
    call add(b, '" text-anchor="middle">')
end subroutine svg_begin

!> End a drawing.
pure subroutine svg_end(b)
    type(buffer_t), intent(inout) :: b  !! the buffer

    call add(b, '</g>'//new_line('a')//'</svg>'//new_line('a'))
end subroutine svg_end

!> Add the rounded box of a node. `outline` may be empty.
pure subroutine svg_box(b, x, y, width, height, fill, outline, pen)
    type(buffer_t), intent(inout) :: b      !! the buffer
    real(wp), intent(in) :: x               !! the left side
    real(wp), intent(in) :: y               !! the top side, downwards
    real(wp), intent(in) :: width           !! the width
    real(wp), intent(in) :: height          !! the height
    character(len=*), intent(in) :: fill    !! the fill
    character(len=*), intent(in) :: outline !! the colour of the outline
    real(wp), intent(in) :: pen             !! the width of the outline

    call add(b, '<rect fill="')
    call add_xml(b, fill)
    call add(b, '" stroke="')
    if (len(outline) > 0) then
        call add_xml(b, outline)
        call add(b, '" stroke-width="')
        call add_tenths(b, pen)
    else
        call add(b, 'none')
    end if
    call add(b, '" rx="7.2" x="')
    call add_tenths(b, x)
    call add(b, '" y="')
    call add_tenths(b, y)
    call add(b, '" width="')
    call add_tenths(b, width)
    call add(b, '" height="')
    call add_tenths(b, height)
    call add(b, '"/>')
end subroutine svg_box

!> Add text. An empty `fill` means the colour of the group around it.
pure subroutine svg_text(b, x, y, text, font_size, fill)
    type(buffer_t), intent(inout) :: b      !! the buffer
    real(wp), intent(in) :: x               !! the centre
    real(wp), intent(in) :: y               !! the baseline, downwards
    character(len=*), intent(in) :: text    !! the text
    real(wp), intent(in) :: font_size       !! the size, or 0 for that of the group
    character(len=*), intent(in) :: fill    !! the colour

    call add(b, '<text x="')
    call add_tenths(b, x)
    call add(b, '" y="')
    call add_tenths(b, y)
    call add(b, '"')
    if (font_size > 0) then
        call add(b, ' font-size="')
        call add_tenths(b, font_size)
        call add(b, '"')
    end if
    if (len(fill) > 0) then
        call add(b, ' fill="')
        call add_xml(b, fill)
        call add(b, '"')
    end if
    call add(b, '>')
    call add_xml(b, text)
    call add(b, '</text>')
end subroutine svg_text

!> Add a curve, with `y` going up in `points`.
pure subroutine svg_curve(b, points, colour, dashed)
    type(buffer_t), intent(inout) :: b      !! the buffer
    type(point_t), intent(in) :: points(:)  !! the control points
    character(len=*), intent(in) :: colour  !! the colour, or empty for that of the group
    logical, intent(in) :: dashed           !! whether the curve is dashed
    integer :: k

    call add(b, '<path')
    if (len(colour) > 0) then
        call add(b, ' stroke="')
        call add_xml(b, colour)
        call add(b, '"')
    end if
    if (dashed) call add(b, ' stroke-dasharray="5,2"')
    call add(b, ' d="M')
    call add_tenths(b, points(1)%x)
    call add(b, ',')
    call add_tenths(b, -points(1)%y)
    call add(b, 'C')
    do k = 2, size(points)
        if (k > 2) call add(b, ' ')
        call add_tenths(b, points(k)%x)
        call add(b, ',')
        call add_tenths(b, -points(k)%y)
    end do
    call add(b, '"/>')
end subroutine svg_curve

!> Add an arrowhead, with `y` going up in `points`.
pure subroutine svg_arrow(b, points, fill, colour)
    type(buffer_t), intent(inout) :: b      !! the buffer
    type(point_t), intent(in) :: points(3)  !! the corners
    character(len=*), intent(in) :: fill    !! the fill
    character(len=*), intent(in) :: colour  !! the colour, or empty for that of the group
    integer :: k

    call add(b, '<polygon fill="')
    call add_xml(b, fill)
    call add(b, '"')
    if (len(colour) > 0) then
        call add(b, ' stroke="')
        call add_xml(b, colour)
        call add(b, '"')
    end if
    call add(b, ' points="')
    do k = 1, 3
        if (k > 1) call add(b, ' ')
        call add_tenths(b, points(k)%x)
        call add(b, ',')
        call add_tenths(b, -points(k)%y)
    end do
    call add(b, '"/>')
end subroutine svg_arrow

end module fdot_svg
