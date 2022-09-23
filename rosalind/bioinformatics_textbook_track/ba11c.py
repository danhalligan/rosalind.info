# Convert a Peptide into a Peptide Vector

from .ba4c import mass


# add fake masses for tests.
def modified_mass():
    masses = mass()
    # add fake masses for tests.
    masses["X"] = 4
    masses["Z"] = 5
    return masses


def peptide_mass(peptide, masses):
    return sum(masses[x] for x in peptide)


def peptide2vector(peptide, masses):
    vec = [0] * (peptide_mass(peptide, masses) + 1)
    for i in range(len(peptide) + 1):
        vec[peptide_mass(peptide[:i], masses)] = 1
    return vec[1:]


def main(file):
    peptide = open(file).read().rstrip()
    print(*peptide2vector(peptide, modified_mass()))
