# Calculating Protein Mass

from .helpers import Parser
from rosalind.helpers import aa_mass


def protein_mass(prot):
    aam = aa_mass()
    return sum([aam[x] for x in prot])


def main(file):
    m = protein_mass(Parser(file).line())
    print(round(m, 3))
