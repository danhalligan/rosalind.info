import os
import sys
from typer import Argument, Option, run, Exit
from pathlib import Path
from typing import Optional
from importlib import resources, import_module


def split_path(string):
    head, tail = os.path.split(string)
    if head == "/":
        return [string]
    else:
        return split_path(head) + [tail]


def module_root():
    return resources.files("rosalind")


# guess the location of a problem based on name
def find_location(problem):
    paths = module_root().glob(f"**/{problem}.py")
    if not paths:
        sys.exit(f"Could not find solution for problem {problem}")
    path = list(paths)[0]
    return split_path(path)[-2]


def test_file(problem):
    loc = find_location(problem)
    return module_root().joinpath(
        "resources", "test-data", loc, f"rosalind_{problem}.txt"
    )


def solve(
    problem: str,
    path: Optional[str] = Argument(None),
    test: bool = Option(False, help="Run with test data."),
):
    """Solve a given Rosalind problem"""
    loc = str(find_location(problem))
    if test:
        path = test_file(problem)
    if not path:
        print("An input file is required if not in test mode!")
        raise Exit()
    if not Path(path).is_file():
        print(f"{path} is not a file!")
        raise Exit()
    module = f"rosalind.{loc}.{problem}"
    module = import_module(module)
    getattr(module, "main")(path)


def main():
    run(solve)
