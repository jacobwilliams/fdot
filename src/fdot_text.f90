module fdot_text
!! The width of text, from the table of Helvetica widths that Graphviz keeps for builds without Pango,
!! in `lib/common/textspan_lut.c`, which is under the Eclipse Public License 2.0.
use iso_fortran_env, only: real64
implicit none
private

public :: text_width, text_height

integer, parameter :: wp = real64                   !! the kind of reals
integer, parameter :: UNITS_PER_EM = 2048           !! the units of the table
real(wp), parameter :: LINESPACING = 1.20_wp        !! the height of a line, as a part of the font size

integer, parameter :: WIDTHS(0:127) = [ &          !! the width of each ASCII character, or -1 for none
          -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, &
          -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, &
          -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  569,  569,  727, 1139, &
        1139, 1821, 1366,  391,  682,  682,  797, 1196,  569,  682,  569,  569, &
        1139, 1139, 1139, 1139, 1139, 1139, 1139, 1139, 1139, 1139,  569,  569, &
        1196, 1196, 1196, 1139, 2079, 1366, 1366, 1479, 1479, 1366, 1251, 1593, &
        1479,  569, 1024, 1366, 1139, 1706, 1479, 1593, 1366, 1593, 1479, 1366, &
        1251, 1479, 1366, 1933, 1366, 1366, 1251,  569,  569,  569,  961, 1139, &
         682, 1139, 1139, 1024, 1139, 1139,  569, 1139, 1139,  455,  455, 1024, &
         455, 1706, 1139, 1139, 1139, 1139,  682, 1024,  569, 1139, 1024, 1479, &
        1024, 1024, 1024,  684,  532,  684, 1196,   -1]

contains

!> The width of `text` in points. A byte that is not ASCII counts as a space, as it does in Graphviz.
pure function text_width(text, font_size) result(width)
    character(len=*), intent(in) :: text    !! the text
    real(wp), intent(in) :: font_size       !! the size of the font, in points
    real(wp) :: width
    integer :: k, code, total

    total = 0
    do k = 1, len(text)
        code = iachar(text(k:k))
        if (code < 0 .or. code > 127) code = iachar(' ')
        if (WIDTHS(code) > 0) total = total + WIDTHS(code)
    end do
    width = font_size*(real(total, wp)/UNITS_PER_EM)
end function text_width

!> The height of a line of text in points.
pure function text_height(font_size) result(height)
    real(wp), intent(in) :: font_size   !! the size of the font, in points
    real(wp) :: height

    height = font_size*LINESPACING
end function text_height

end module fdot_text
