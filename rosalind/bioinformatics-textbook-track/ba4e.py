# Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

from .ba4c import mass, cyclo_spectrum
from .ba1a import substrings


def expand(peptides, masses):
    return [p + [x] for x in masses for p in list(peptides)]


def is_consistent(peptide, target):
    m = [sum(x) for i in range(1, len(peptide) + 1) for x in substrings(peptide, i)]
    return all(x in target for x in m)


def cyclopeptide_sequencing(target):
    masses = set(mass().values())
    peptides = [[]]
    while len(peptides):
        peptides = expand(peptides, masses)
        for peptide in list(peptides):
            if sum(peptide) == target[-1]:
                if list(cyclo_spectrum(peptide)) == target:
                    yield peptide
                peptides.remove(peptide)
            elif not is_consistent(peptide, target):
                peptides.remove(peptide)


def main(file):
    mass = open(file).read().split()
    target = list(map(int, mass))
    res = ["-".join(str(x) for x in p) for p in cyclopeptide_sequencing(target)]
    print(*res)
