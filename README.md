# fdot

The graph layout of [Graphviz](https://graphviz.org) `dot`, in modern Fortran.

`fdot` is a port of the part of `dot` that [fdoc](https://github.com/jacobwilliams/fdoc) needs. You give
it nodes and edges, and get an SVG drawing back as a string. There is no DOT language, no plugins, no
files and no external process.

- Standard Fortran 2018, with no compiler extensions and no C library routines.
- Safe to call from several threads at once, on different graphs.
- The same input gives the same output.

It draws directed graphs of named rounded boxes, from left to right or right to left, with cycles,
self-loops, edge labels and concentrated edges. It is not a general replacement for Graphviz: clusters,
ports, `rank=same`, record and HTML labels and other shapes and arrows are left out.

## Use

```toml
[dependencies]
fdot = { git = "https://github.com/jacobwilliams/fdot" }
```

```fortran
use fdot, only: graph_t, RANK_RL

type(graph_t) :: g
character(len=:), allocatable :: svg
integer :: a, b, stat

call g%initialize(direction=RANK_RL, concentrate=.true.)
call g%add_node('json_module', a, fill='#3657c9', link='../module/json_module.html')
call g%add_node('json_kinds', b, fill='#646d7e')
call g%add_edge(a, b, colour='#8a93a5')
call g%layout(stat)
call g%write_svg(svg)
```

- Nodes are integers returned by `add_node`.
- `layout` and `write_svg` are separate, so that a caller can read the result with `node_position`,
  `edge_curve` and `edge_label_position`, and draw it another way.
- Colours and links are passed through as text.
- `stat` is `FDOT_OK` or says what went wrong.
- `write_dot` writes the graph in the DOT language, for debugging.
- `node_style_t` and `edge_style_t` change sizes. The defaults are what fdoc draws.

## Development

All tools come from the pixi environment, and fpm builds the code.

| Task | What it does |
|---|---|
| `pixi run build` | Build the library |
| `pixi run test` | Run the unit tests |
| `pixi run test-threads` | Lay out graphs with 1, 4 and 10 threads, and compare |
| `pixi run robust` | Lay out random graphs, with runtime checks |
| `pixi run check-standard` | Check that the code is standard Fortran |
| `pixi run check-threads` | Check for what is not safe in threads, and for I/O and C |
| `pixi run register-dot` | Register the plugins of `dot`, once, before the next two |
| `pixi run corpus FDOC NAME=FILE ...` | Export the graphs of fdoc as a corpus |
| `pixi run compare` | Compare with `dot` on the corpus |

`fdot-layout`, in `app/`, lays out graphs described on standard input, for these tools.

## How it relates to Graphviz

The port keeps the algorithms of Graphviz and their order, so that results stay comparable. The state
moves from global variables and records on graph objects to a derived type passed through the calls,
in integer-indexed arrays, and recursive searches on graphs use explicit stacks.

| Module | Port of |
|---|---|
| `fdot_simplex` | `lib/common/ns.c` |
| `fdot_rank` | `lib/dotgen/rank.c`, `acyclic.c`, `class1.c` |
| `fdot_structure` | `lib/dotgen/fastgr.c`, `decomp.c` |
| `fdot_order` | `lib/dotgen/mincross.c`, `class2.c` |
| `fdot_concentrate` | `lib/dotgen/conc.c` |
| `fdot_position` | `lib/dotgen/position.c` |
| `fdot_route` | `lib/dotgen/dotsplines.c`, `lib/common/routespl.c`, `splines.c` |
| `fdot_path` | `lib/pathplan/shortest.c`, `route.c`, `solvers.c`, `util.c`, `triang.c` |
| `fdot_text` | `lib/common/textspan_lut.c`, `textspan.c` |
| `fdot_draw` | Parts of `lib/common/shapes.c`, `arrows.c`, `postproc.c` |
| `fdot_svg` | Parts of `plugin/core/gvrender_core_svg.c`, `lib/common/emit.c` |

The reference is Graphviz at commit `c83e61e`.

Text is measured with Graphviz's table of Helvetica widths and not with Pango, so node widths differ a
little from those of `dot`, and so do the coordinates that follow from them.

## How it is checked

`pixi run compare` runs on the 6,786 graphs that fdoc generates for fpm, json-fortran, stdlib and
qt-fortran, with `dot` 14.1.2, and measures fdot against it:

| Check | Result |
|---|---|
| Nodes do not overlap, edges reach their nodes, ranks follow the edges | All graphs |
| Ranks | Equal in 6,659 graphs, and equal in total edge length in all |
| Order of nodes in each rank | Equal in all 6,659 graphs with equal ranks |
| Edge length, fdot over `dot` | Median 0.995, worst 1.18 |
| Area, fdot over `dot` | Median 0.979, worst 1.34 |

The graphs that are furthest from `dot` come from ties. `dot` leaves some ties to `qsort`, whose order
differs between C libraries, while fdot sorts stably and gives the same result everywhere. Text widths
matter too, so `compare --table` runs `dot` without Pango, and `tools/stable_qsort.c` gives it a stable
`qsort`. With both, ranks and orders are equal in all 6,786 graphs, and the worst ratio is 1.000 for both
edge length and area, so the differences above do not come from the port. The curves then match to the
precision that `dot` prints, which is five significant digits.

Layout takes about half the time of `dot`, on one thread.

## License

`fdot` is a derivative work of Graphviz, and is under the [Eclipse Public License 2.0](LICENSE), like
Graphviz. Graphviz is Copyright (c) 2011 AT&T Intellectual Property, and its contributors are listed at
https://graphviz.org.
