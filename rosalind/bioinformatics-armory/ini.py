# Introduction to the Bioinformatics Armory

from Bio.Seq import Seq


def ini(seq):
    seq = Seq(seq)
    return [seq.count(base) for base in ("A", "C", "G", "T")]


def main(file):
    seq = str(open(file).read())
    print(*ini(seq))
