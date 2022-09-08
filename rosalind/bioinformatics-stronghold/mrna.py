# Inferring mRNA from Protein

from functools import reduce
from .helpers import Parser
from rosalind.helpers import codons


def mrna(seq, mod=1000000):
    cod = codons()
    seq = seq + "*"
    return reduce(lambda p, c: p * cod[c] % mod, seq, 1)


def main(file):
    print(mrna(Parser(file).line()))
