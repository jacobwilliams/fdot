#!/usr/bin/env python
"""Lay out random graphs with the runtime checks of the debug build, as section 6 of the plan describes.

    pixi run robust [--count N] [--profile debug]

The graphs have deep cycles, long chains, very dense parts, self-loops, parallel edges, labels and
several components. All must end. The last line has a checksum of the output, from a fixed sequence of
graphs, which CI compares between platforms.
"""

import argparse
import hashlib
import random
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TIME_LIMIT = 300


def node_line(rng):
    name = "".join(rng.choice("abcdefghij_klmnop%XYZ") for _ in range(rng.randint(1, 24)))
    link = "page.html" if rng.random() < 0.5 else ""
    outline = "#f0a868" if rng.random() < 0.05 else ""
    return f"N\t{name}\t#646d7e\t{link}\t{outline}"


def edge_line(rng, a, b, labels):
    label = "".join(rng.choice("abcdef ") for _ in range(rng.randint(1, 16))) if rng.random() < labels else ""
    return f"E\t{a}\t{b}\t{label.strip()}\t#8a93a5"


def graph(rng, kind):
    """The lines of a random graph of a kind."""
    n = rng.randint(1, 40)
    labels = rng.choice([0, 0, 0.1, 0.5])
    edges = []
    if kind == "cycles":
        n = rng.randint(2, 120)
        edges = [(k, k % n + 1) for k in range(1, n + 1)]
        edges += [(rng.randint(1, n), rng.randint(1, n)) for _ in range(n)]
    elif kind == "chain":
        n = rng.randint(2, 400)
        edges = [(k, k + 1) for k in range(1, n)]
        edges += [(1, n), (rng.randint(1, n), rng.randint(1, n))]
    elif kind == "dense":
        n = rng.randint(2, 30)
        edges = [(a, b) for a in range(1, n + 1) for b in range(1, n + 1) if rng.random() < 0.6]
    elif kind == "star":
        n = rng.randint(2, 200)
        edges = [(1, k) for k in range(2, n + 1)] + [(k, n) for k in range(2, n)]
    elif kind == "layers":
        width, depth = rng.randint(2, 12), rng.randint(2, 12)
        n = width * depth
        for d in range(depth - 1):
            for w in range(width):
                for _ in range(rng.randint(0, 4)):
                    edges.append((d * width + w + 1, (d + 1) * width + rng.randint(1, width)))
    elif kind == "parts":
        n = rng.randint(2, 60)
        edges = [(a, a + rng.randint(1, 3)) for a in range(1, n - 3) if rng.random() < 0.4]
    else:
        edges = [(rng.randint(1, n), rng.randint(1, n)) for _ in range(rng.randint(0, 3 * n))]
    rng.shuffle(edges)
    direction = rng.choice(["LR", "RL"])
    lines = [f"G\t{direction}\t{rng.choice([0, 1, 1])}"]
    lines += [node_line(rng) for _ in range(n)]
    lines += [edge_line(rng, a, b, labels) for a, b in edges]
    lines.append("L")
    return lines


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=700)
    parser.add_argument("--profile", default="debug")
    options = parser.parse_args()
    rng = random.Random(20261101)
    kinds = ["cycles", "chain", "dense", "star", "layers", "parts", "random"]
    lines = []
    for k in range(options.count):
        lines += graph(rng, kinds[k % len(kinds)])
    text = "\n".join(lines) + "\n"
    subprocess.run(["fpm", "build", "--profile", options.profile], cwd=ROOT, check=True, capture_output=True)
    out = subprocess.run(["fpm", "run", "fdot-layout", "--profile", options.profile, "--runner", "echo"],
                         cwd=ROOT, check=True, capture_output=True, text=True).stdout
    binary = ROOT / out.split()[-1]
    start = time.time()
    result = subprocess.run([str(binary), "svg"], input=text.encode(), capture_output=True, timeout=TIME_LIMIT)
    seconds = time.time() - start
    if result.returncode != 0:
        print(result.stderr.decode(errors="replace")[-3000:])
        sys.exit(f"fdot-layout failed with {result.returncode}")
    output = result.stdout.replace(b"\r\n", b"\n")
    lines = output.split(b"\n")
    failed = sum(line.startswith(b"failed ") for line in lines)
    done = sum(line.strip() == b"%%%" for line in lines)
    print(f"{done} of {options.count} graphs in {seconds:.1f} s, {failed} failed")
    print(f"checksum {hashlib.sha256(output).hexdigest()}")
    sys.exit(1 if failed or done != options.count else 0)


if __name__ == "__main__":
    main()
