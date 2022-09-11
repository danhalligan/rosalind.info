# Base Filtration by Quality

from Bio import SeqIO
import sys


def main(file):
    handle = open(file)
    q = int(next(handle))
    for x in SeqIO.parse(handle, "fastq"):
        quals = x.letter_annotations["phred_quality"]
        for start in range(len(quals)):
            if quals[start] >= q:
                break
        for end in range(len(quals) - 1, 0, -1):
            if quals[end] >= q:
                break
        SeqIO.write(x[start : end + 1], sys.stdout, "fastq")
