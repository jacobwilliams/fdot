# fdot: the graph layout fdoc needs, in modern Fortran

`fdot` is a new fpm library, under the Eclipse Public License 2.0: a port of the part of Graphviz's `dot` layout that fdoc
uses. fdoc calls it directly with nodes and edges, and gets an SVG drawing back as a string. There is
no DOT language, no plugins, no files and no external process.

Reference: Graphviz at https://gitlab.com/graphviz/graphviz/, commit `c83e61e` (2026-11-01). Sizes in
this plan come from that commit. fdoc's pixi env has `dot` 14.1.2, which is what results are compared to.

## 1. Requirements

- An fpm package, in its own repository, under EPL-2.0, usable by any Fortran project.
- Standard Fortran, with no compiler extensions, checked with `-std=f2018`. It must be able to build
  with other compilers such as ifx, though only gfortran is tested for now.
- Runs on Windows, macOS and Linux, tested in CI on each.
- No C library routines. Section 4 shows that none are needed.
- SVG is the only output for now.
- No file I/O and no shell commands in the library. Text goes in and out through strings.
- Safe to call from several threads at once, on different graphs: no module variables, no `save`, no
  internal `read` or `write`, and no function returning a deferred-length string.
- The same input gives the same output, on every platform and with any number of threads.
- All tests run in the project's pixi env, which is the configuration that counts.

Not goals: being a general replacement for Graphviz, reading DOT files, or byte-identical output.

## 2. Why

What fdoc does now, from `fdoc_graphs`, `fdoc_dot` and `fdoc_svg`:

1. Writes a `.gv` file for each graph.
2. Starts `dot` on batches of 25 files, and polls for marker files, which is why it needs `usleep`.
3. Reads each SVG back, and parses it to make it smaller.

The port removes all of that. In return fdoc loses a required external program, its last use of
`usleep`, and several of its remaining Windows risks: `dot` batches are expected to run one at a time
there, and shell commands and quoting are needed to run them. Layout becomes an ordinary parallel loop
over graphs. `fdoc_dot` and `fdoc_svg` can be deleted, and the SVG is written in its compact form to
begin with.

Layout is about 1.1 s of stdlib's 3.8 s today, on 10 cores, so speed is a side benefit and not the
reason. Being a single binary with no dependencies is the reason.

## 3. What fdoc actually needs

Measured on the graphs fdoc generates for four projects:

| Project | Graphs | Nodes, median / 95% / max | Edges, median / 95% / max | Cyclic graphs | Self-loops | Edge labels |
|---|---:|---:|---:|---:|---:|---:|
| fpm | 761 | 10 / 56 / 360 | 16 / 235 / 804 | 127 | 0 | 449 |
| json-fortran | 535 | 59 / 225 / 388 | 70 / 317 / 712 | 335 | 3 | 9 |
| qt-fortran | 2721 | 3 / 9 / 967 | 2 / 8 / 975 | 0 | 0 | 790 |
| stdlib | 2769 | 2 / 15 / 250 | 1 / 25 / 476 | 0 | 0 | 15 |

No graph has parallel edges between the same pair of nodes.

**Graph attributes:** `rankdir` of `LR` or `RL` only, `concentrate=true` always, and a transparent
background.

**Node attributes:** a filled box with rounded corners, `height=0.3`, `margin="0.1,0.04"`, 13 pt
Helvetica in white, a fill colour, an optional outline with `penwidth=2` for the root of the graph, and
a link.

**Edge attributes:** a colour, `arrowsize=0.7`, and sometimes dashed with an 11 pt label.

**Never used:** clusters, ports, `rank=same`, record and HTML labels, other node shapes, other arrow
types, `splines` settings, `ordering`, images, or any output but SVG.

What this means for the algorithms:

- **Cycles are common,** so cycle breaking and reversed edges are required.
- **Self-loops are rare but real,** so they are needed.
- **Edge labels are needed.** `dot` lays them out as nodes on doubled ranks.
- **Concentration is needed,** and it is what makes the densest graphs readable.
- **Flat edges are not,** as they arise from `rank=same`, which fdoc never uses. Parts of `flat.c` and
  flat-edge routing in `dotsplines.c` can be left out, once M3 has confirmed this.

## 4. What is ported

About 12,500 lines of C, out of about 58,000 in the Graphviz libraries that `dot` needs, and less once
the parts for clusters and flat edges are left out.

| Graphviz source | Lines | Port | Becomes |
|---|---:|---|---|
| `lib/common/ns.c` | 1479 | All: network simplex | `fdot_simplex` |
| `lib/dotgen/rank.c` | 1113 | Without clusters | `fdot_rank` |
| `lib/dotgen/acyclic.c`, `decomp.c` | 200 | All | `fdot_rank` |
| `lib/dotgen/class1.c`, `class2.c`, `fastgr.c` | 650 | Without clusters | `fdot_structure` |
| `lib/dotgen/conc.c` | 248 | All | `fdot_concentrate` |
| `lib/dotgen/mincross.c` | 1795 | Without clusters and flat edges | `fdot_order` |
| `lib/dotgen/position.c` | 1145 | Without clusters | `fdot_position` |
| `lib/dotgen/dotsplines.c` | 2330 | Regular, reversed and self-loop edges | `fdot_route` |
| `lib/common/routespl.c`, parts of `splines.c` | 1500 | What `dotsplines.c` calls | `fdot_route` |
| `lib/pathplan/shortest.c`, `route.c`, `solvers.c`, `util.c` | 1000 | All | `fdot_path` |
| `lib/common/textspan_lut.c`, `textspan.c` | 850 | Helvetica width table only | `fdot_text` |
| Parts of `shapes.c`, `arrows.c`, `postproc.c` | 800 | Rounded box, `normal` arrow, translation | `fdot_draw` |
| `plugin/core/gvrender_core_svg.c`, parts of `emit.c` | 500 | The elements fdoc keeps | `fdot_svg` |

Left out entirely: `cgraph` and `cdt` (replaced by plain arrays), `gvc` and the plugin system, the DOT
parser, `cluster.c`, `compound.c`, `sameport.c`, `aspect.c`, `htmltable.c`, the other layout engines,
`ortho`, and all other renderers.

**C library routines.** Standard Fortran covers all of this, as the port is arithmetic on arrays and
text. The one thing Graphviz gets from outside is text measurement, through Pango here, and the port
uses the table of Helvetica widths that Graphviz keeps for builds without Pango. So the library has no
`bind(C)`.

**The design is not a line-by-line copy.** Graphviz keeps its state in global variables and in records
attached to graph objects, and that cannot be thread-safe. The algorithms and their order are kept, so
that results stay comparable, while the state moves to a derived type passed through the calls, stored
in integer-indexed arrays.

## 5. Interface

```fortran
use fdot, only: graph_t, node_style_t, edge_style_t, RANK_LR, RANK_RL

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

- Nodes are integers returned by `add_node`, so fdoc does not need to build names for them.
- `layout` and `write_svg` are separate, so a caller can read the positions and draw them another way.
- Colours and links are passed through as text, and the library does not interpret them.
- `write_dot` writes the graph in the DOT language, for debugging and for the comparison in section 6.
  There is no DOT parser.

## 6. How it is checked

Byte-for-byte comparison with `dot` is impossible, and not wanted: `dot` here measures text with Pango
and the port with a table, so node widths differ, and everything after follows from those. So results
are compared by what makes a drawing good. `dot` stays in the pixi env as a development tool for this.

- **Corpus:** the graphs fdoc generates for the four projects in section 3, about 6,800 in all.
- **Validity, which must always hold:** nodes do not overlap, each edge starts and ends at its nodes,
  and ranks follow the edges other than the reversed ones.
- **Quality, compared with `dot` on the same graph:** edge crossings, total edge length and drawing
  area. The median over the corpus should be within 5% of `dot`, and no graph more than 25% worse.
- **Ranking,** which involves no text sizes: for graphs without labels, node ranks should match `dot`
  exactly. Network simplex can have several optimal answers, so where they differ, the total edge length
  must match instead.
- **Determinism:** each graph is laid out twice, and with 1 and with 10 threads, and the output must be
  identical. CI compares the output of all three platforms.
- **Robustness:** random graphs, with deep cycles, long chains and very dense parts, must end, with no
  bounds errors under the debug build's runtime checks.
- **fdoc:** `test-projects` passes, and sample pages are looked at in a browser in both themes.

## 7. Milestones

**M0. Foundation.** New repository, fpm package, pixi env, CI on three platforms, and a `check-standard`
step. `graph_t` with `add_node`, `add_edge` and `write_dot`. A tool in fdoc to export its graphs as the
corpus, and the checker for section 6 in place before there is anything to check.

**M1. Text and sizes.** The width table, and node and label sizes from `margin` and `height`.

**M2. Ranking.** Network simplex, cycle breaking and connected components, with edge-label rank
doubling. This is the first real comparison with `dot`, using node ranks.

**M3. Ordering.** Virtual nodes for long edges, and crossing reduction. Confirm here that flat edges are
never needed on the corpus. Compare crossings with `dot`.

**M4. Positioning.** The auxiliary graph and second network simplex in `position.c`. Compare drawing area
and edge length.

**M5. Edge routing.** Shortest paths and spline fitting from `pathplan`, boxes from `dotsplines.c`,
reversed edges, self-loops, and labels. This is the largest and riskiest milestone.

**M6. Concentration.** Merging edges, in `conc.c`, which needs routing to be working. This comes after
M5 so that its effect can be measured alone.

**M7. Drawing.** SVG for rounded boxes, arrows, links, dashed edges and labels, for the compact form that
`fdoc_svg` produces today, so that fdoc's CSS and zoom keep working.

**M8. Quality and speed.** All of section 6 passes. Stress tests with threads. Compare speed with `dot`.
Release a first version.

**M9. Use in fdoc.** Add it as an fpm dependency, behind a setting, with `dot` still the default. Compare
sites over `test-projects`. Then make it the default, and remove `fdoc_dot`, `fdoc_svg` and the `usleep`
binding. `graph_dir` then has no files to keep, so it is accepted and ignored, with a warning. Keep
`graphviz` in the pixi env only as a development tool.

M1 to M3 give a usable check of whether the port is on track, before the expensive parts. Each
milestone is compared with `dot` as it lands, so differences are found near their cause.

## 8. Risks

| Risk | Handling |
|---|---|
| Licenses of the two projects differ | fdot is EPL-2.0 and fdoc is MIT, so they stay in separate repositories, and fdoc's releases say where fdot's source is (section 9) |
| Edge routing is intricate, with many special cases in `dotsplines.c` | Most cases are for clusters, ports and flat edges, which are left out. Port it last among the layout stages, and stop for a review if the corpus shows gaps |
| Width table differs from what browsers draw | `dot` has this problem too. fdoc already asks for Helvetica and a fallback. Add a small margin if text overflows in testing |
| Recursion on large graphs overflows a thread's stack | Use explicit stacks in the depth-first searches of `ns.c` and in cycle breaking, which are recursive in C |
| Floating-point results differ between platforms | Ranks and order use integers. Keep positions in integers as `dot` does in network simplex, and round coordinates to one decimal when written |
| Crossing reduction gets slow on dense graphs | `dot` has limits for this, `mclimit` and `nslimit`, so port them with the same defaults |
| Corpus misses something a project will hit | `write_dot` makes any failing graph easy to reduce and compare. Random graphs cover the rest |
| fdoc's graphs change, and need a feature left out | Add that feature then. The design keeps the left-out parts easy to add, as with clusters |

## 9. Decisions

1. **License: EPL-2.0.** Graphviz is under the Eclipse Public License 2.0, and a port of its code is a
   derivative work, so fdot uses the same license and keeps Graphviz's notices. fdoc stays MIT and
   links to fdot, which is why they are separate repositories. Binary releases of fdoc must say where
   fdot's source is, and that is worth confirming against the license before the first release.

2. **Name: `fdot`.**

3. **Output: SVG only, for now.** Others could come later from the positions that `layout` leaves.
