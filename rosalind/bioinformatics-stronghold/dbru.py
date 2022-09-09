# Constructing a De Bruijn Graph

from .helpers import Parser, Dna


def dbru(seqs, rev=True):
    """Constructing a De Bruijn Graph"""
    seqs = list(seqs)
    if rev:
        # add reverse complement sequences
        seqs = set(seqs).union([Dna(seq).revc().seq for seq in seqs])
    return [(x[:-1], x[1:]) for x in seqs]


def main(file):
    for pair in sorted(dbru(Parser(file).lines())):
        print("(", pair[0], ", ", pair[1], ")", sep="")
