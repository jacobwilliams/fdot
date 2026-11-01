#!/usr/bin/env python
"""Compare fdot with `dot` on the corpus, as section 6 of the plan describes.

    pixi run compare [--list N] [NAME ...]

Validity is checked for every graph: nodes do not overlap, edges reach their nodes, and ranks follow the
edges. Then ranks and orders, which do not depend on text sizes, are compared exactly, and crossings, edge
length and area are compared as ratios of fdot to `dot`, given as median, 95th percentile and worst. The
results of `dot` are kept next to the corpus, so a later run needs only fdot. Exits with 1 if a required
check fails.

With `--table`, `dot` runs without its Pango plugin, so that it measures text with the same table as
fdot. Nothing then separates the two but the port, so this also reports how many graphs have the same
coordinates. `dot` leaves some ties to `qsort`, so to rule those out as well, see `tools/stable_qsort.c`.
"""

import argparse
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "corpus"
NODE = re.compile(r'^(n\d+) \[(.*)\];$')
EDGE = re.compile(r'^(n\d+) -> (n\d+) \[(.*)\];$')
ATTR = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,\s\]]+)')
BATCH = 50
EXACT = 0.05    # points: what counts as the same coordinate
FULL = ".full"


def unquote(text):
    if text.startswith('"'):
        text = text[1:-1]
    return re.sub(r"\\(.)", r"\1", text)


def parse(path):
    """The direction, nodes and edges of a graph that fdoc wrote."""
    direction, nodes, edges = "LR", {}, []
    for line in path.read_text(errors="replace").splitlines():
        match = EDGE.match(line)
        if match:
            attrs = {k: unquote(v) for k, v in ATTR.findall(match[3])}
            edges.append((match[1], match[2], attrs))
            continue
        match = NODE.match(line)
        if match:
            nodes[match[1]] = {k: unquote(v) for k, v in ATTR.findall(match[2])}
        elif line.startswith("rankdir="):
            direction = line[8:10]
    return direction, nodes, edges


def clean(text):
    return text.replace("\t", " ").replace("\n", " ")


def script(graph):
    """The input of fdot-layout for a graph."""
    direction, nodes, edges = graph
    ids = {name: k + 1 for k, name in enumerate(nodes)}
    lines = [f"G\t{direction}\t1"]
    for name, attrs in nodes.items():
        lines.append("\t".join(["N", clean(attrs.get("label", name)), attrs.get("fillcolor", ""), clean(attrs.get("URL", "")),
                                attrs.get("color", "") if "penwidth" in attrs else ""]))
    for tail, head, attrs in edges:
        lines.append("\t".join(["E", str(ids[tail]), str(ids[head]), clean(attrs.get("label", "")), attrs.get("color", "#8a93a5")]))
    lines.append("L")
    return "\n".join(lines) + "\n"


def plugins():
    """A directory with the plugins of `dot`, but not Pango, and the environment to use it."""
    source = Path(shutil.which("dot")).resolve().parent.parent / "lib" / "graphviz"
    target = CORPUS / ".plugins"
    target.mkdir(exist_ok=True)
    for pattern in ("libgvplugin_core*", "libgvplugin_dot_layout*"):
        for path in source.glob(pattern):
            if not (target / path.name).exists():
                shutil.copy2(path, target / path.name, follow_symlinks=False)
    env = {**os.environ, "GVBINDIR": str(target)}
    subprocess.run(["dot", "-c"], env=env, check=True)
    return env


def run_dot(args):
    """Run `dot` on a batch of graphs, for what is not yet there."""
    paths, full, env = args
    for suffix, extra in ((full, []), (".rank", ["-Gphase=1"]), (".order", ["-Gphase=2"])):
        for path in paths:
            target = path.with_suffix(path.suffix + suffix)
            if not target.exists():
                subprocess.run(["dot", "-Tdot", *extra, "-o", str(target), str(path)], check=True, capture_output=True, env=env)


def statements(path):
    """The statements of a graph that `dot` wrote."""
    text = path.read_text(errors="replace").replace("\\\n", "")
    return re.split(r";\n", text)


def read_dot(path, full_suffix):
    """Positions, sizes, ranks, orders and curves from `dot`, indexed as in the graph."""
    full = {"nodes": {}, "curves": [], "bb": (0.0, 0.0)}
    rank = {}
    for name, target in (("full", full_suffix), ("rank", ".rank"), ("order", ".order")):
        for part in statements(path.with_suffix(path.suffix + target)):
            part = part.strip()
            match = re.match(r'^(n\d+)\s+\[(.*)$', part, re.S)
            if match:
                attrs = dict(ATTR.findall(match[2]))
                if name == "full":
                    x, y = map(float, unquote(attrs["pos"]).split(","))
                    full["nodes"][match[1]] = (x, y, 72 * float(unquote(attrs["width"])), 72 * float(unquote(attrs.get("height", "0.3"))))
                elif name == "rank":
                    rank[match[1]] = int(unquote(attrs.get("rank", "0")))
                else:
                    full.setdefault("order", {})[match[1]] = int(unquote(attrs.get("order", "0")))
                continue
            match = re.match(r'^(n\d+) -> (n\d+)\s+\[(.*)$', part, re.S)
            if match and name == "full":
                attrs = dict(ATTR.findall(match[3]))
                if "pos" in attrs:
                    full["curves"].append((match[1], match[2], parse_pos(unquote(attrs["pos"]))))
            box = re.search(r'bb="[\d.-]+,[\d.-]+,([\d.]+),([\d.]+)"', part)
            if box and name == "full":
                full["bb"] = (float(box[1]), float(box[2]))
    full["rank"] = rank
    return full


def parse_pos(text):
    """The curves of a `pos` attribute."""
    curves = []
    for piece in text.split(";"):
        points = [tuple(map(float, item.split(","))) for item in piece.split() if not item.startswith(("e,", "s,"))]
        curves.append(points)
    return curves


def sample(points):
    """A curve as a polyline."""
    out = []
    for k in range(0, len(points) - 3, 3):
        p = points[k:k + 4]
        for i in range(0, 5):
            t = i / 4
            s = 1 - t
            out.append((s*s*s*p[0][0] + 3*s*s*t*p[1][0] + 3*s*t*t*p[2][0] + t*t*t*p[3][0],
                        s*s*s*p[0][1] + 3*s*s*t*p[1][1] + 3*s*t*t*p[2][1] + t*t*t*p[3][1]))
    return out


def metrics(curves):
    """Crossings and length of a list of (tail, head, polyline)."""
    seg, owner, length = [], [], 0.0
    ends = []
    for index, (tail, head, line) in enumerate(curves):
        ends.append((tail, head))
        for a, b in zip(line, line[1:]):
            seg.append((*a, *b))
            owner.append(index)
            length += math.hypot(b[0] - a[0], b[1] - a[1])
    if not seg:
        return 0, length
    s = np.array(seg)
    o = np.array(owner)
    keys = {}
    tails = np.zeros(len(curves), dtype=np.int64)
    heads = np.zeros(len(curves), dtype=np.int64)
    for k, (t, h) in enumerate(ends):
        tails[k] = keys.setdefault(t, len(keys))
        heads[k] = keys.setdefault(h, len(keys))
    crossed = set()
    dx = s[:, 2] - s[:, 0]
    dy = s[:, 3] - s[:, 1]
    for i in range(len(s) - 1):
        j = np.arange(i + 1, len(s))
        j = j[o[j] != o[i]]
        if j.size == 0:
            continue
        rx = s[j, 0] - s[i, 0]
        ry = s[j, 1] - s[i, 1]
        d = dx[i] * dy[j] - dy[i] * dx[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (rx * dy[j] - ry * dx[j]) / d
            u = (rx * dy[i] - ry * dx[i]) / d
        hit = j[(d != 0) & (t > 0.001) & (t < 0.999) & (u > 0.001) & (u < 0.999)]
        for k in np.unique(o[hit]):
            a, b = int(min(o[i], k)), int(max(o[i], k))
            share = {tails[a], heads[a]} & {tails[b], heads[b]}
            if not share:
                crossed.add((a, b, int(i)))
    return len(crossed), length


def read_fdot(lines):
    """The layout from fdot-layout."""
    out = {"nodes": {}, "curves": {}, "labels": {}, "bb": (0.0, 0.0), "failed": False}
    for line in lines:
        f = line.split()
        if f[0] == "graph":
            out["bb"] = (float(f[1]), float(f[2]))
        elif f[0] == "node":
            out["nodes"][int(f[1])] = (int(f[2]), int(f[3]), *map(float, f[4:8]))
        elif f[0] == "edge":
            v = list(map(float, f[3:]))
            out["curves"].setdefault(int(f[1]), []).append(list(zip(v[0::2], v[1::2])))
        elif f[0] == "label":
            out["labels"][int(f[1])] = (float(f[2]), float(f[3]))
        elif f[0] == "failed":
            out["failed"] = True
    return out


def check_graph(args):
    """All the checks for one graph."""
    path, lines, full_suffix = args
    graph = parse(path)
    _, nodes, edges = graph
    names = list(nodes)
    ids = {name: k + 1 for k, name in enumerate(names)}
    mine = read_fdot(lines)
    theirs = read_dot(path, full_suffix)
    result = {"path": str(path), "problems": [], "nodes": len(names)}
    if mine["failed"]:
        result["problems"].append("layout failed")
        return result

    # validity
    pos = mine["nodes"]
    items = sorted(pos.values(), key=lambda v: v[2])
    for i, a in enumerate(items):
        for b in items[i + 1:]:
            if b[2] - b[4] / 2 >= a[2] + a[4] / 2:
                break
            if abs(a[2] - b[2]) < (a[4] + b[4]) / 2 - 0.01 and abs(a[3] - b[3]) < (a[5] + b[5]) / 2 - 0.01:
                result["problems"].append("nodes overlap")
    labelled = any("label" in attrs for _, _, attrs in edges)
    touched = set()
    for k, (tail, head, attrs) in enumerate(edges, 1):
        t, h = pos[ids[tail]], pos[ids[head]]
        least = 2 if labelled else 1
        if tail != head and abs(t[0] - h[0]) < least:
            result["problems"].append(f"edge {k} is {abs(t[0] - h[0])} ranks long")
        for curve in mine["curves"].get(k, []):
            for point, node in ((curve[0], t), (curve[-1], h)):
                inside = abs(point[0] - node[2]) <= node[4] / 2 + 12 and abs(point[1] - node[3]) <= node[5] / 2 + 12
                if inside:
                    touched.add(ids[tail] if node is t else ids[head])
        curves = mine["curves"].get(k, [])
        if len(curves) == 1:
            c = curves[0]
            near_t = abs(c[0][0] - t[2]) <= t[4] / 2 + 12 and abs(c[0][1] - t[3]) <= t[5] / 2 + 12
            near_h = abs(c[-1][0] - h[2]) <= h[4] / 2 + 12 and abs(c[-1][1] - h[3]) <= h[5] / 2 + 12
            merged = any(len(mine["curves"].get(j, [])) != 1 for j in range(1, len(edges) + 1))
            if not merged and not (near_t and near_h):
                result["problems"].append(f"edge {k} does not reach its nodes")
    connected = {ids[e[0]] for e in edges} | {ids[e[1]] for e in edges}
    if connected - touched:
        result["problems"].append(f"{len(connected - touched)} nodes without an edge")

    # ranks and orders, which do not depend on text
    rank_d = {n: theirs["rank"].get(n) for n in names}
    rank_f = {n: pos[ids[n]][0] for n in names}
    result["ranks_equal"] = all(rank_d[n] == rank_f[n] for n in names)
    length = lambda r: sum(abs(r[t] - r[h]) for t, h, _ in edges)
    result["rank_length_equal"] = length(rank_d) == length(rank_f)
    order_d = theirs.get("order", {})
    # relative, as `dot` stops this phase before concentration removes virtual nodes
    if result["ranks_equal"]:
        by_rank = {}
        for n in names:
            by_rank.setdefault(rank_f[n], []).append(n)
        result["orders_equal"] = all(sorted(group, key=lambda n: order_d.get(n, 0)) == sorted(group, key=lambda n: pos[ids[n]][1])
                                     for group in by_rank.values())
    else:
        result["orders_equal"] = None

    # the same coordinates
    same = abs(mine["bb"][0] - theirs["bb"][0]) <= 1.01 and abs(mine["bb"][1] - theirs["bb"][1]) <= 1.01
    for n in names:
        p, q = pos[ids[n]], theirs["nodes"][n]
        same = same and abs(p[2] - q[0]) <= EXACT and abs(p[3] - q[1]) <= EXACT
    result["same_nodes"] = same
    wanted = {}
    for tail, head, parts in theirs["curves"]:
        wanted.setdefault((tail, head), []).append(parts)
    same_curves = True
    for k, (tail, head, _) in enumerate(edges, 1):
        have = mine["curves"].get(k, [])
        options_d = wanted.get((tail, head), [])
        matched = not have and not options_d
        for parts in options_d:
            # the order of the parts of an edge is a tie in `edgecmp`, which `qsort` leaves open
            first = lambda c: sorted(c, key=lambda q: q[0])
            if len(parts) == len(have) and all(len(a) == len(b) and all(abs(u[0] - v[0]) <= EXACT and abs(u[1] - v[1]) <= EXACT
                                                                       for u, v in zip(a, b))
                                               for a, b in zip(first(parts), first(have))):
                matched = True
        same_curves = same_curves and matched
    result["same_curves"] = same_curves

    # quality
    curves_f = []
    for k, (tail, head, _) in enumerate(edges, 1):
        for curve in mine["curves"].get(k, []):
            curves_f.append((tail, head, sample(curve)))
    curves_d = []
    for tail, head, parts in theirs["curves"]:
        for curve in parts:
            curves_d.append((tail, head, sample(curve)))
    cross_f, len_f = metrics(curves_f)
    cross_d, len_d = metrics(curves_d)
    result.update(cross_f=cross_f, cross_d=cross_d, len_f=len_f, len_d=len_d,
                  area_f=mine["bb"][0] * mine["bb"][1], area_d=theirs["bb"][0] * theirs["bb"][1], edges=len(edges))
    return result


def layout_all(paths, binary):
    """The output of fdot-layout for each graph, and how long it took."""
    text = "".join(script(parse(path)) for path in paths)
    start = time.time()
    result = subprocess.run([binary], input=text, capture_output=True, text=True, errors="replace")
    seconds = time.time() - start
    if result.returncode != 0:
        print(result.stderr[-2000:])
        sys.exit(f"fdot-layout failed with {result.returncode}")
    groups, group = [], []
    for line in result.stdout.splitlines():
        if line.strip() == "end":
            groups.append(group)
            group = []
        else:
            group.append(line)
    if len(groups) != len(paths):
        sys.exit(f"{len(groups)} results for {len(paths)} graphs")
    return groups, seconds


def binary_path():
    subprocess.run(["fpm", "build", "--profile", "release"], cwd=ROOT, check=True, capture_output=True)
    out = subprocess.run(["fpm", "run", "fdot-layout", "--profile", "release", "--runner", "echo"],
                         cwd=ROOT, check=True, capture_output=True, text=True).stdout
    return str(ROOT / out.split()[-1])


def spread(values):
    values = sorted(values)
    if not values:
        return "no data"
    return f"median {statistics.median(values):.3f}, 95% {values[int(0.95 * (len(values) - 1))]:.3f}, worst {values[-1]:.3f}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="*")
    parser.add_argument("--list", type=int, default=0, help="show the N worst graphs for each measure")
    parser.add_argument("--table", action="store_true", help="run dot with the same table of text widths")
    parser.add_argument("--different", type=int, default=0, help="show N graphs with other coordinates")
    options = parser.parse_args()
    names = options.names or sorted(p.name for p in CORPUS.iterdir() if p.is_dir())
    if not shutil.which("dot"):
        sys.exit("dot is not on the path")
    binary = binary_path()
    env = plugins() if options.table else None
    full = ".table" if options.table else FULL
    all_results = []
    with ProcessPoolExecutor() as pool:
        for name in names:
            paths = sorted((CORPUS / name).rglob("*.gv"))
            todo = [paths[k:k + BATCH] for k in range(0, len(paths), BATCH)]
            list(pool.map(run_dot, [(batch, full, env) for batch in todo]))
            lines, seconds = layout_all(paths, binary)
            results = list(pool.map(check_graph, [(p, l, full) for p, l in zip(paths, lines)], chunksize=16))
            print(f"{name}: {len(paths)} graphs, fdot took {seconds:.2f} s", flush=True)
            all_results += results

    bad = [r for r in all_results if r["problems"]]
    print(f"\nvalidity: {len(bad)} of {len(all_results)} graphs have problems")
    for r in bad[:10]:
        print(f"  {r['path']}: {r['problems'][:3]}")
    good = [r for r in all_results if "ranks_equal" in r]
    ranks = sum(r["ranks_equal"] for r in good)
    lengths = sum(r["rank_length_equal"] for r in good)
    comparable = [r for r in good if r["orders_equal"] is not None]
    orders = sum(r["orders_equal"] for r in comparable)
    print(f"ranks: {ranks} of {len(good)} equal, {lengths} with equal total length")
    print(f"orders of real nodes: {orders} of {len(comparable)} graphs with equal ranks are equal")
    mismatch = [r for r in good if not r["rank_length_equal"]]
    for r in mismatch[:5]:
        print(f"  rank length differs: {r['path']}")
    print(f"coordinates: nodes the same in {sum(r['same_nodes'] for r in good)} of {len(good)}, "
          f"all curves as well in {sum(r['same_nodes'] and r['same_curves'] for r in good)}")
    different = [r for r in good if not (r["same_nodes"] and r["same_curves"]) and r["ranks_equal"]]
    print(f"  with equal ranks: {sum(not r['same_nodes'] for r in different)} differ in nodes, {len(different)} in all")
    for r in sorted(different, key=lambda r: (r["same_nodes"], r["nodes"]))[:options.different]:
        print(f"  different: {r['path']} nodes {r['nodes']} same nodes {r['same_nodes']}")
    useful = [r for r in good if r["edges"] > 0 and r["len_d"] > 0]
    measures = {
        "crossings, +1 each": [(r["cross_f"] + 1) / (r["cross_d"] + 1) for r in useful],
        "edge length": [r["len_f"] / r["len_d"] for r in useful],
        "area": [r["area_f"] / r["area_d"] for r in useful if r["area_d"] > 0],
    }
    print(f"crossings in all: fdot {sum(r['cross_f'] for r in useful)}, dot {sum(r['cross_d'] for r in useful)}")
    failed = bool(bad) or lengths != len(good)
    for label, values in measures.items():
        print(f"{label}: {spread(values)}")
    median_ok = all(abs(statistics.median(v) - 1) <= 0.05 for v in measures.values() if v)
    worst_ok = all(max(v) <= 1.25 for k, v in measures.items() if v and k != "crossings, +1 each")
    for label, values in measures.items():
        over = sum(v > 1.25 for v in values)
        print(f"{label}: {over} of {len(values)} graphs over 25% worse")
    print(f"quality: median within 5%: {median_ok}; no graph over 25% worse in length or area: {worst_ok}")
    if options.list:
        for label, key in (("crossings", lambda r: r["cross_f"] - r["cross_d"]), ("length", lambda r: r["len_f"] / r["len_d"]),
                           ("area", lambda r: r["area_f"] / max(r["area_d"], 1))):
            print(f"worst for {label}:")
            for r in sorted(useful, key=key, reverse=True)[:options.list]:
                print(f"  {r['path']} nodes {r['nodes']} cross {r['cross_f']}/{r['cross_d']} "
                      f"len {r['len_f']:.0f}/{r['len_d']:.0f} area {r['area_f']:.0f}/{r['area_d']:.0f}")
    sys.exit(1 if failed or not median_ok else 0)


if __name__ == "__main__":
    main()
