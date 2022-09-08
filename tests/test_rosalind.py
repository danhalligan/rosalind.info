import pytest
import os
import re
from importlib import import_module
import random

# Gather list of solution scripts
locations = [
    "algorithmic-heights",
    "bioinformatics-armory",
    "bioinformatics-textbook-track",
    "python-village",
]
problems = []
for location in locations:
    for x in os.listdir("rosalind/" + location):
        if re.match("(helpers|solve).+", x):
            continue
        m = re.match(r"([^_]+)\.py$", x)
        if m:
            problems.append([location, m.groups()[0]])


# Run each command with a test file and check our snapshot
# matches the downloaded Sample Output / "expected" version
@pytest.mark.parametrize("problem", problems)
def test_cli_function(capfd, snapshot, problem):
    path = "." + ".".join(problem)
    module = import_module(path, package="rosalind")
    test_file = f"tests/data/{problem[0]}/rosalind_{problem[1]}.txt"
    random.seed(42)
    getattr(module, "main")(test_file)
    out, _ = capfd.readouterr()
    snap_file = f"{problem[0]}/{problem[1]}.txt"
    snapshot.snapshot_dir = "tests/snapshots"
    snapshot.assert_match(out, snap_file)
