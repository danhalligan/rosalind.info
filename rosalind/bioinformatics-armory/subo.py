# Suboptimal Local Alignment

from Bio import SeqIO
import sys


def hamm(s1, s2):
    return sum(xi != yi for xi, yi in zip(s1, s2))


def count_hamming(pattern, seq, dist=3):
    count = 0
    np = len(pattern)
    ns = len(seq)
    for i in range(ns - np + 1):
        if hamm(seq[i : i + np], pattern) <= dist:
            count += 1
    return count


def main(file, pattern="GACTCCTTTGTTTGCCTTAAATAGATACATATTT"):
    seqs = [x.seq for x in SeqIO.parse(file, "fasta")]
    print(*[count_hamming(pattern, seq) for seq in seqs])


if __name__ == "__main__":
    pattern = input("Enter to 32-40 bp pattern identified by Lalign: ")
    main(sys.argv[1], pattern)
