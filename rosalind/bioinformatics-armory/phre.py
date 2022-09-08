# Read Quality Distribution

from Bio import SeqIO
import sys


def main(file):
    handle = open(file)
    n = int(next(handle))
    count = 0
    for rcrd in SeqIO.parse(handle, "fastq"):
        q = rcrd.letter_annotations["phred_quality"]
        count += sum(q) / len(q) < n
    print(count)


if __name__ == "__main__":
    main(sys.argv[1])
