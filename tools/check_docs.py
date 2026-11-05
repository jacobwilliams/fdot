#!/usr/bin/env python
"""Check that every routine, argument, type, type component and parameter of fdot is documented.

A routine needs an ``!>`` block before it, and each argument an inline ``!!`` on its declaration (or a
``!!`` line straight after it). A type needs a ``!!`` line straight after its header, and each component,
type-bound procedure included, an inline ``!!``, and so does each named constant. Components and constants
have one name to a declaration, so that each has its own.
Run in the pixi env: ``pixi run check-docs``.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HEADER = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure|module|integer|logical|real|complex|character|type)"
    r"(?:\s*\([^()]*\))?\s+)*"
    r"(subroutine|function)\s+(\w+)\s*(?:\(([^)]*)\))?",
    re.IGNORECASE,
)
END = re.compile(r"^\s*end\s+(subroutine|function)\b", re.IGNORECASE)
TYPE = re.compile(r"^\s*type\s*(?:,[^:]*)?::\s*(\w+)|^\s*type\s+(\w+)\s*$", re.IGNORECASE)
END_TYPE = re.compile(r"^\s*end\s*type\b", re.IGNORECASE)


def declared(line: str):
    """The names declared on ``line``, if it is a declaration"""
    code = line.split("!")[0]
    if "::" not in code:
        return []
    names = re.sub(r"\([^()]*\)|\[[^\]]*\]", "", code.split("::", 1)[1])
    return [n.split("=")[0].strip().lower() for n in names.split(",") if n.strip()]


def check_routines(path: Path, lines):
    problems = []
    for i, line in enumerate(lines):
        match = HEADER.match(line)
        if not match or END.match(line):
            continue
        name = match.group(2)
        where = f"{path.relative_to(ROOT)}:{i + 1}: {name}"
        if i == 0 or not lines[i - 1].lstrip().startswith("!>"):
            problems.append(f"{where}: no !> block before it")
        args = [a.strip().lower() for a in (match.group(3) or "").split(",") if a.strip()]
        for arg in args:
            for j in range(i + 1, len(lines)):
                if END.match(lines[j]):
                    problems.append(f"{where}: argument {arg} not declared")
                    break
                names = declared(lines[j])
                if arg in names:
                    following = j + 1 < len(lines) and lines[j + 1].lstrip().startswith("!!")
                    if len(names) > 1:
                        problems.append(f"{where}: argument {arg} shares a declaration")
                    elif "!!" not in lines[j] and not following:
                        problems.append(f"{where}: argument {arg} is not documented")
                    break
    return problems


def check_types(path: Path, lines):
    problems = []
    i = 0
    while i < len(lines):
        match = TYPE.match(lines[i].split("!")[0])
        if not match:
            i += 1
            continue
        name = match.group(1) or match.group(2)
        where = f"{path.relative_to(ROOT)}:{i + 1}: type {name}"
        if i + 1 >= len(lines) or not lines[i + 1].lstrip().startswith("!!"):
            problems.append(f"{where}: no !! line after its header")
        j = i + 1
        while j < len(lines) and not END_TYPE.match(lines[j]):
            names = declared(lines[j])
            if names:
                if len(names) > 1:
                    problems.append(f"{where}: components {', '.join(names)} share a declaration")
                elif "!!" not in lines[j]:
                    problems.append(f"{where}: component {names[0]} is not documented")
            j += 1
        i = j + 1
    return problems


def check_parameters(path: Path, lines):
    problems = []
    for i, line in enumerate(lines):
        code = line.split("!")[0]
        if not re.search(r",\s*parameter\b", code.split("::")[0], re.IGNORECASE) or "::" not in code:
            continue
        names = declared(line)
        where = f"{path.relative_to(ROOT)}:{i + 1}"
        if len(names) > 1:
            problems.append(f"{where}: parameters {', '.join(names)} share a declaration")
        elif "!!" not in line:
            problems.append(f"{where}: parameter {names[0]} is not documented")
    return problems


def main() -> int:
    problems = []
    for directory in ("app", "src", "test", "example"):
        for path in sorted((ROOT / directory).glob("*.f90")):
            lines = path.read_text(encoding="utf-8").splitlines()
            problems += check_routines(path, lines) + check_types(path, lines) + check_parameters(path, lines)
    print("\n".join(problems))
    print(f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
