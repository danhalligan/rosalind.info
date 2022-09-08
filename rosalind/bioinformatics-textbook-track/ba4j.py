# Generate the Theoretical Spectrum of a Linear Peptide

from .ba4g import linear_spectrum
from .ba4c import mass


def main(file):
    peptide = open(file).read().rstrip()
    masses = mass()
    print(*linear_spectrum([masses[x] for x in peptide]))
