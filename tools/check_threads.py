#!/usr/bin/env python
"""Check what makes the library safe in threads, and free of C and of the outside world.

The library must have no module variables, no ``save``, no internal or file I/O, no shell commands, no
``bind(C)`` and no function with a deferred-length result. gfortran keeps the length of such a result in
a static ``slen.N`` variable that threads share, so the objects are checked for those as well. Run in the
pixi env: ``pixi run check-threads``.
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BANNED = [
    (re.compile(r"^\s*(read|write|open|close|inquire|print)\b\s*[(*]", re.I), "I/O"),
    (re.compile(r"\bsave\b", re.I), "save"),
    (re.compile(r"\bbind\s*\(\s*c", re.I), "bind(C)"),
    (re.compile(r"\bexecute_command_line\b", re.I), "a shell command"),
    (re.compile(r"\biso_c_binding\b", re.I), "iso_c_binding"),
    (re.compile(r"\bcommon\s*/", re.I), "a common block"),
    (re.compile(r"^\s*character\(len=:\), allocatable\s*::.*\bfunction\b", re.I), "a deferred-length function result"),
]
DECLARATION = re.compile(r"^\s*(integer|real|logical|character|type\s*\()[^!]*::", re.I)


def source_problems(path):
    """Problems in the text of a file."""
    problems = []
    in_contains = False
    in_type = 0
    for number, line in enumerate(path.read_text().splitlines(), 1):
        code = line.split("!")[0]
        for pattern, what in BANNED:
            if pattern.search(code):
                problems.append(f"{path.name}:{number}: {what}: {line.strip()}")
        low = code.strip().lower()
        if low == "contains" and not in_type:
            in_contains = True
        if re.match(r"type\s*(,.*)?::", low) or re.match(r"type\s+\w+\s*$", low):
            in_type += 1
        if low.startswith("end type"):
            in_type -= 1
        if not in_contains and not in_type and DECLARATION.match(code) and "parameter" not in low:
            problems.append(f"{path.name}:{number}: a module variable: {line.strip()}")
    return problems


def object_problems():
    """Static string lengths in the objects of the library."""
    problems = []
    objects = sorted(ROOT.glob("build/gfortran_*/fdot/*.f90.o"))
    if not objects:
        return ["no object files: pixi run build"]
    for path in objects:
        symbols = subprocess.run(["nm", str(path)], check=True, capture_output=True, text=True).stdout
        count = sum("slen." in line for line in symbols.splitlines())
        if count:
            problems.append(f"{path.name}: {count} static string length(s)")
    return problems


def main():
    problems = object_problems()
    for path in sorted((ROOT / "src").glob("*.f90")):
        problems += source_problems(path)
    for problem in problems:
        print(problem)
    print(f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
