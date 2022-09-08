# Protein Translation

from Bio.Seq import Seq


def main(file):
    dna, prot = open(file).read().splitlines()
    for i in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15]:
        if prot == Seq(dna).translate(table=i, to_stop=True):
            print(i)
            break
