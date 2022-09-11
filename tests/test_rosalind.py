import pytest
import os
import re
from importlib import import_module
import random
import yaml
from unittest.mock import patch


# Gather list of solution scripts
locations = [
    "python_village",
    "bioinformatics_stronghold",
    "bioinformatics_armory",
    "bioinformatics_textbook_track",
    "algorithmic_heights",
]
problems = []
for location in locations:
    for x in os.listdir("rosalind/" + location):
        if re.match("(helpers|solve).+", x):
            continue
        m = re.match(r"([^_]+)\.py$", x)
        if m:
            problems.append([location, m.groups()[0]])


# Dictionary of saved uniprot output
uniprot = yaml.safe_load(open("tests/uniprot_output.yaml"))


# Run each command with a test file and check our snapshot
# matches the downloaded Sample Output / "expected" version
@pytest.mark.parametrize("problem", problems)
def test_cli_function(capfd, snapshot, problem):
    path = "." + ".".join(problem)
    module = import_module(path, package="rosalind")
    test_file = f"rosalind/resources/test-data/{problem[0]}/rosalind_{problem[1]}.txt"
    random.seed(42)
    with patch(
        "rosalind.bioinformatics_stronghold.mprt.uniprot_output",
        side_effect=uniprot.get,
    ):
        getattr(module, "main")(test_file)
    out, _ = capfd.readouterr()
    snap_file = f"{problem[0]}/{problem[1]}.txt"
    snapshot.snapshot_dir = "tests/snapshots"
    snapshot.assert_match(out, snap_file)
