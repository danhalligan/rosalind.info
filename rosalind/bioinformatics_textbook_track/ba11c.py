# Convert a Peptide into a Peptide Vector

import os
from .ba4c import mass


# masses with fake masses for tests if we're in test mode
def masses():
    if "ROSALIND_TEST" in os.environ:
        return {"X": 4, "Z": 5}
    else:
        return mass()


def peptide_mass(peptide, masses):
    return sum(masses[x] for x in peptide)


def peptide2vector(peptide, masses):
    vec = [0] * (peptide_mass(peptide, masses) + 1)
    for i in range(len(peptide) + 1):
        vec[peptide_mass(peptide[:i], masses)] = 1
    return vec[1:]


def main(file):
    peptide = open(file).read().rstrip()
    print(*peptide2vector(peptide, masses()))
