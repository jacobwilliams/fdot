#!/usr/bin/env python
"""Export the graphs that fdoc generates for real projects, as the corpus of section 6 of the plan.

    pixi run corpus FDOC_BINARY NAME=PROJECT_FILE ...

Each project is run from the directory of its project file, with fdoc's `graph_dir` setting, which keeps
the `.gv` files. They end up in `corpus/NAME`. The sites go to a temporary directory, and nothing is
written to the projects.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CORPUS = ROOT / "corpus"


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    binary = str(Path(sys.argv[1]).expanduser().resolve())
    for spec in sys.argv[2:]:
        name, project = spec.split("=", 1)
        project = Path(project).expanduser().resolve()
        target = CORPUS / name
        shutil.rmtree(target, ignore_errors=True)
        target.mkdir(parents=True)
        with tempfile.TemporaryDirectory() as site:
            command = [binary, "-q", "-g", "-o", site, "--config", f"graph_dir='{target}'", str(project)]
            result = subprocess.run(command, cwd=project.parent, capture_output=True, text=True)
        for extra in list(target.rglob("*.svg")):
            extra.unlink()
        count = len(list(target.rglob("*.gv")))
        print(f"{name}: {count} graphs (fdoc exited with {result.returncode})", flush=True)
        if result.returncode != 0 or count == 0:
            print(result.stderr[-1000:])


if __name__ == "__main__":
    main()
