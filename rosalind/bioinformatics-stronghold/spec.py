# Inferring Protein from Spectrum

from math import isclose
from rosalind.helpers import aa_mass
from .helpers import Parser
from .prtm import protein_mass


def match_mass(weight, rel_tol=1e-6):
    aam = aa_mass()
    matches = [k for k in aam if isclose(weight, aam[k], rel_tol=rel_tol)]
    return None if len(matches) == 0 else matches[0]


def spectrum(seq):
    """Complete spectrum"""
    r = range(1, len(seq))
    spec = [seq] + [seq[:k] for k in r] + [seq[k:] for k in r]
    return [protein_mass(x) for x in spec]


def main(file):
    weights = [float(x) for x in Parser(file).lines()]
    diff = [j - i for i, j in zip(weights[:-1], weights[1:])]
    print("".join([match_mass(x) for x in diff]))
