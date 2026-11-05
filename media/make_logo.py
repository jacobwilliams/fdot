#!/usr/bin/env python
"""Make the logo of fdoc: a small graph of the kind fdoc draws, whose root, at the top left, holds the F of the
Fortran logo.

    python media/make_logo.py [OUTPUT]

The F is read from `fortran-logo.svg`, next to this script, and the logo is written to OUTPUT, by default
`logo.svg` in the same directory. Change the numbers below to change the design.
"""

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SIZE = 128          # the logo is a square of this size
CORNER = 24         # the radius of its corners
BACK = "#ffffff"    # the background
BORDER = "#e3dfe9"  # a faint border, so that the square shows on a white page
EDGE = "#8a93a5"    # the colour of edges, as in fdoc
WIDTH = 3.5         # the width of edges
ARROW = 8           # the length of arrowheads

# the root of the graph, at the top left, with fdoc's outline for a root, holds the F in the purple of
# the Fortran logo
ROOT = (10, 10, 38, 44)             # x, y, width, height
OUTLINE = "#f0a868"
F_HEIGHT = 28                       # the height of the F, which is centred in the root

# the other nodes: name, x, y, width, height, fill
NODES = [
    ("a", 66, 24, 30, 16, "#3657c9"),
    ("b", 13, 78, 32, 16, "#1d7a55"),
    ("c", 82, 94, 32, 16, "#0b7285"),
]

# the edges, between sides of nodes or points with a direction, with an arrowhead where they end if the
# third value is true. The edges into JOIN meet there, as concentrated edges do, and go on as one.
JOIN = (98, 80, "down")
EDGES = [
    (("root", "right"), ("a", "left"), True),
    (("root", "bottom"), ("b", "top"), True),
    (("a", "bottom"), JOIN, False),
    (("b", "right"), JOIN, False),
    (JOIN, ("c", "top"), True),
]

DIRECTIONS = {"left": (-1, 0), "right": (1, 0), "top": (0, -1), "bottom": (0, 1), "up": (0, -1), "down": (0, 1)}


def the_f(path):
    """The path of the F in the Fortran logo, the transform it is drawn with, and the purple of the logo."""
    text = path.read_text()
    element = re.search(r"<path\b[^>]*>", text, re.S)[0]
    d = re.search(r'\sd="([^"]*)"', element)[1]
    group = re.search(r'<g\b[^>]*transform="translate\(([-\d.]+),([-\d.]+)\)"', text)
    shift = (float(group[1]), float(group[2])) if group else (0.0, 0.0)
    purple = re.search(r"<rect\b[^>]*fill:(#[0-9a-fA-F]{6})", text, re.S)[1]
    return d, shift, purple


def extent(d):
    """The smallest box holding all the points of path data, control points included."""
    tokens = re.findall(r"[A-Za-z]|-?(?:\d+\.?\d*|\.\d+)(?:e-?\d+)?", d)
    xs, ys = [], []
    x = y = start_x = start_y = 0.0
    command = ""
    i = 0
    counts = {"m": 2, "l": 2, "h": 1, "v": 1, "c": 6, "z": 0}
    while i < len(tokens):
        if tokens[i].isalpha():
            command = tokens[i]
            i += 1
            if command in "zZ":
                x, y = start_x, start_y
                continue
        lower = command.lower()
        relative = command.islower()
        values = [float(v) for v in tokens[i:i + counts[lower]]]
        i += counts[lower]
        if lower in "ml":
            x, y = (x + values[0], y + values[1]) if relative else values
            if lower == "m":
                start_x, start_y = x, y
                command = "l" if relative else "L"  # pairs after a move are lines
        elif lower == "h":
            x = x + values[0] if relative else values[0]
        elif lower == "v":
            y = y + values[0] if relative else values[0]
        elif lower == "c":
            for k in range(0, 6, 2):
                px, py = (x + values[k], y + values[k + 1]) if relative else (values[k], values[k + 1])
                xs.append(px)
                ys.append(py)
            x, y = xs[-1], ys[-1]
        else:
            raise ValueError(f"path command {command} is not handled")
        xs.append(x)
        ys.append(y)
    return min(xs), min(ys), max(xs), max(ys)


def point(end):
    """A point and the direction an edge goes there: a side of a node, or a point with a direction."""
    if isinstance(end[0], (int, float)):
        x, y, direction = end
        return (x, y), DIRECTIONS[direction]
    name, side = end
    if name == "root":
        x, y, w, h = ROOT
    else:
        _, x, y, w, h, _ = next(n for n in NODES if n[0] == name)
    middle = {"left": (x, y + h / 2), "right": (x + w, y + h / 2), "top": (x + w / 2, y), "bottom": (x + w / 2, y + h)}
    return middle[side], DIRECTIONS[side]


def edge(start, end, arrow):
    """The path of an edge, which leaves and arrives along the directions of its ends, and its arrowhead
    or `None`."""
    (x0, y0), (dx0, dy0) = point(start)
    (x1, y1), (dx1, dy1) = point(end)
    if isinstance(end[0], str):
        dx1, dy1 = -dx1, -dy1   # an edge into a side of a node goes into the node
    tip = (x1, y1)
    if arrow:
        x1, y1 = x1 - dx1 * ARROW, y1 - dy1 * ARROW     # the curve stops where the arrowhead starts
    bend = max(abs(x1 - x0), abs(y1 - y0)) / 2
    curve = (f"M{x0:g} {y0:g}C{x0 + dx0 * bend:g} {y0 + dy0 * bend:g} "
             f"{x1 - dx1 * bend:g} {y1 - dy1 * bend:g} {x1:g} {y1:g}")
    head = None
    if arrow:
        half = ARROW * 0.45
        bx, by = x1 - dx1 * 0.5, y1 - dy1 * 0.5     # a little back, to cover the end of the curve
        head = (f"M{bx - dy1 * half:g} {by + dx1 * half:g}L{tip[0]:g} {tip[1]:g}"
                f"L{bx + dy1 * half:g} {by - dx1 * half:g}z")
    return curve, head


def logo(d, shift, purple):
    """The text of the logo."""
    left, top, right, bottom = extent(d)
    scale = F_HEIGHT / (bottom - top)
    # the path is drawn shifted, as in its own file, then scaled and centred in the root
    rx, ry, rw, rh = ROOT
    tx = rx + (rw - (right - left) * scale) / 2 - (left + shift[0]) * scale
    ty = ry + (rh - F_HEIGHT) / 2 - (top + shift[1]) * scale
    curves, heads = zip(*(edge(*e) for e in EDGES))
    nodes = "\n".join(f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" fill="{fill}"/>'
                      for _, x, y, w, h, fill in NODES)
    edges = "\n".join(f'    <path d="{c}"/>' for c in curves)
    arrows = "\n".join(f'    <path d="{h}"/>' for h in heads if h)
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE} {SIZE}" role="img" aria-label="fdoc">
  <rect x="1" y="1" width="{SIZE - 2}" height="{SIZE - 2}" rx="{CORNER - 1}" fill="{BACK}" stroke="{BORDER}" stroke-width="2"/>
  <!-- a graph, as fdoc draws them, whose root holds the F of the Fortran logo -->
  <g fill="none" stroke="{EDGE}" stroke-width="{WIDTH}" stroke-linecap="round">
{edges}
  </g>
  <g fill="{EDGE}">
{arrows}
  </g>
  <rect x="{rx}" y="{ry}" width="{rw}" height="{rh}" rx="7" fill="{purple}" stroke="{OUTLINE}" stroke-width="3"/>
  <path transform="translate({tx:.3f} {ty:.3f}) scale({scale:.4f}) translate({shift[0]} {shift[1]})" fill="#fff"
        d="{d}"/>
{nodes}
</svg>
"""


def main():
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "logo.svg"
    d, shift, purple = the_f(HERE / "fortran-logo.svg")
    output.write_text(logo(d, shift, purple))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
