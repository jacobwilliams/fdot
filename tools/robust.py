#!/usr/bin/env python
"""Lay out random graphs with the runtime checks of the debug build, as section 6 of the plan describes.

    pixi run robust [--count N] [--profile debug]

The graphs have deep cycles, long chains, very dense parts, self-loops, parallel edges, labels and
several components. All must end. The last line has a checksum of the output, from a fixed sequence of
graphs, which CI compares between platforms.

Graphs run in batches, in parallel, each with its own time limit. A batch that goes over the limit is
named, and its input is kept in `robust-timeout.txt`, so that a hang can be told from a slow machine.
"""

import argparse
import hashlib
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BATCH = 25          # graphs in a process
TIME_LIMIT = 300    # seconds for a batch


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
    parser.add_argument("--count", type=int, default=210)
    parser.add_argument("--profile", default="debug")
    options = parser.parse_args()
    rng = random.Random(20261101)
    kinds = ["cycles", "chain", "dense", "star", "layers", "parts", "random"]
    batches = [[] for _ in range((options.count + BATCH - 1) // BATCH)]
    for k in range(options.count):
        batches[k // BATCH] += graph(rng, kinds[k % len(kinds)])
    subprocess.run(["fpm", "build", "--profile", options.profile], cwd=ROOT, check=True, capture_output=True)
    out = subprocess.run(["fpm", "run", "fdot-layout", "--profile", options.profile, "--runner", "echo"],
                         cwd=ROOT, check=True, capture_output=True, text=True).stdout
    binary = ROOT / out.split()[-1]
    start = time.time()

    def run(numbered):
        number, lines = numbered
        text = ("\n".join(lines) + "\n").encode()
        began = time.time()
        try:
            result = subprocess.run([str(binary), "svg"], input=text, capture_output=True, timeout=TIME_LIMIT)
        except subprocess.TimeoutExpired:
            (ROOT / "robust-timeout.txt").write_bytes(text)
            sys.exit(f"batch {number} of {len(batches)} went over {TIME_LIMIT} s: see robust-timeout.txt")
        print(f"batch {number} of {len(batches)}: {time.time() - began:.1f} s", flush=True)
        if result.returncode != 0:
            print(result.stderr.decode(errors="replace")[-3000:], flush=True)
            sys.exit(f"batch {number}: fdot-layout failed with {result.returncode}")
        return result.stdout

    with ThreadPoolExecutor(os.cpu_count()) as pool:
        output = b"".join(pool.map(run, enumerate(batches, 1)))
    seconds = time.time() - start
    output = output.replace(b"\r\n", b"\n")
    lines = output.split(b"\n")
    failed = sum(line.startswith(b"failed ") for line in lines)
    done = sum(line.strip() == b"%%%" for line in lines)
    print(f"{done} of {options.count} graphs in {seconds:.1f} s, {failed} failed")
    print(f"checksum {hashlib.sha256(output).hexdigest()}")
    sys.exit(1 if failed or done != options.count else 0)


if __name__ == "__main__":
    main()
