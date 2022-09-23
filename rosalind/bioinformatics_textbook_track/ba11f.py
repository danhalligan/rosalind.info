# Find a Highest-Scoring Peptide in a Proteome against a Spectrum

from .ba11c import masses
from math import inf
from .ba11c import peptide2vector


# This is a brute force scan through all peptides in the proteome
# We can move to the next start position once our peptide candidate mass
# exceeds that implied by the spectral vector (i.e. its length)
def peptide_identification(sv, proteome):
    m = masses()
    best_score = -inf
    best_peptide = ""
    for i in range(len(proteome)):
        for j in range(i + 1, len(proteome)):
            pv = peptide2vector(proteome[i:j], m)
            if len(pv) > len(sv):
                break
            if len(pv) == len(sv):
                score = sum(a * b for a, b in zip(pv, sv))
                if score > best_score:
                    best_score = score
                    best_peptide = proteome[i:j]
    return best_peptide, best_score


def main(file):
    sv, seq = open(file).read().splitlines()
    sv = list(map(int, sv.split()))
    peptide, _ = peptide_identification(sv, seq)
    print(peptide)
