# FASTQ format introduction

from Bio import SeqIO
import sys


def main(file):
    for rcrd in SeqIO.parse(file, "fastq"):
        SeqIO.write(rcrd, sys.stdout, "fasta")


if __name__ == "__main__":
    main(sys.argv[1])
