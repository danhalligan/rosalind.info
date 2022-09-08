import typer
import sys
import re
import importlib
from pathlib import Path


# We could keep a lookup of locations, but actually, we can fairly accurately
# guess the location of a problem based on name:
def find_location(problem):
    paths = Path().glob(f"**/{problem}.py")
    if not paths:
        sys.exit(f"Could not find solution for problem {problem}")
    return list(paths)[0]


def solve(problem: str, path: str):
    """Solve a given Rosalind problem"""
    if not Path(path).is_file():
        sys.exit(f"{path} is not an file!")
    loc = str(find_location(problem))
    module = re.sub("/", ".", re.sub(".py", "", loc))
    module = importlib.import_module(module)
    getattr(module, "main")(path)


def main():
    typer.run(solve)
