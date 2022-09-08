# Compute the Score of a Linear Peptide

from .ba4g import linear_score
from .ba4c import mass


def main(file):
    peptide, spectrum = open(file).read().splitlines()
    spectrum = list(map(int, spectrum.split()))
    masses = mass()
    print(linear_score([masses[x] for x in peptide], spectrum))
