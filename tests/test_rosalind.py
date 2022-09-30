import pytest
import os
import re
import random
import yaml
from importlib import import_module
from unittest.mock import patch
from Bio.Seq import Seq  # noqa: F401
from Bio.SeqRecord import SeqRecord  # noqa: F401


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
        if re.match("(helpers).+", x):
            continue
        m = re.match(r"([^_]+)\.py$", x)
        if m:
            problems.append([location, m.groups()[0]])


# Dictionary of saved uniprot output
uniprot = yaml.safe_load(open("tests/uniprot_output.yaml"))


# Evaluate a set of saved SeqRecord objects
rcrds = yaml.safe_load(open("tests/seq_records.yaml"))
rcrds = {k: eval(v) for k, v in rcrds.items()}


def get_records(ids):
    return [rcrds.get(id) for id in ids]


# Run each command with a test file and check our snapshot
# matches the downloaded Sample Output / "expected" version
ids = [x[1] for x in problems]


@pytest.mark.parametrize("problem", problems, ids=ids)
@patch("rosalind.bioinformatics_stronghold.mprt.uniprot_output", uniprot.get)
@patch("rosalind.bioinformatics_armory.frmt.get_records", get_records)
@patch("rosalind.bioinformatics_armory.gbk.gbk", lambda *_: {"Count": 7})
def test_cli_function(capfd, snapshot, problem):
    os.environ["ROSALIND_TEST"] = "1"
    path = "." + ".".join(problem)
    module = import_module(path, package="rosalind")
    test_file = f"rosalind/resources/test-data/{problem[0]}/rosalind_{problem[1]}.txt"
    random.seed(42)
    getattr(module, "main")(test_file)
    out, _ = capfd.readouterr()
    snap_file = f"{problem[0]}/{problem[1]}.txt"
    snapshot.snapshot_dir = "tests/snapshots"
    snapshot.assert_match(out, snap_file)
