# Introduction to the Bioinformatics Armory

from Bio.Seq import Seq
import sys


def ini(seq):
    seq = Seq(seq)
    return [seq.count(base) for base in ("A", "C", "G", "T")]


def main(file):
    seq = str(open(file).read())
    print(*ini(seq))


if __name__ == "__main__":
    main(sys.argv[1])
