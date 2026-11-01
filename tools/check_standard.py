#!/usr/bin/env python
"""Check that fdot is standard Fortran, with no compiler extensions, so that it can build with others.

Each source file is compiled with ``-std=f2018``, in the order of its dependencies. Lines over 132
characters are errors there too. Run in the pixi env: ``pixi run check-standard``.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ORDER = ["fdot_structure", "fdot_text", "fdot_simplex", "fdot_rank", "fdot_order", "fdot_concentrate", "fdot_position",
         "fdot_path", "fdot_route", "fdot_draw", "fdot_svg", "fdot"]


def main():
    problems = 0
    with tempfile.TemporaryDirectory() as scratch:
        files = [ROOT / "src" / f"{name}.f90" for name in ORDER]
        missing = set((ROOT / "src").glob("*.f90")) - set(files)
        if missing:
            sys.exit(f"not in the list of this tool: {sorted(p.name for p in missing)}")
        files += sorted((ROOT / "example").glob("*.f90")) + sorted((ROOT / "app").glob("*.f90"))
        for path in files:
            command = ["gfortran", "-std=f2018", "-fopenmp", "-Wall", "-Wextra", "-Wno-compare-reals", "-fsyntax-only",
                       f"-I{scratch}", f"-J{scratch}", str(path)]
            result = subprocess.run(command, capture_output=True, text=True, cwd=ROOT)
            if result.returncode != 0:
                problems += 1
                print(f"{path.relative_to(ROOT)}:\n{result.stderr[:800]}")
    print(f"{problems} problem file(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
