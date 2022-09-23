# Convert a Peptide Vector into a Peptide

from .ba11c import masses


def peptide_mass(peptide, masses):
    return sum(masses[x] for x in peptide)


def prefixes2peptide(prefixes, masses):
    m = {v: k for k, v in masses.items()}
    peptide = [m[b - a] for a, b in zip([0] + prefixes, prefixes)]
    return "".join(peptide)


def vector2peptide(vector, masses):
    prefixes = [i + 1 for i, d in enumerate(vector) if d]
    return prefixes2peptide(prefixes, masses)


def main(file):
    vector = list(map(int, open(file).read().split()))
    print(vector2peptide(vector, masses()))
