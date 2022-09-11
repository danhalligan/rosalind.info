# Base Quality Distribution

from Bio import SeqIO
from functools import reduce


def main(file):
    handle = open(file)
    q = int(next(handle))
    qarr = [x.letter_annotations["phred_quality"] for x in SeqIO.parse(handle, "fastq")]
    tots = reduce(lambda a, b: [x + y for x, y in zip(a, b)], qarr)
    print(sum(x / len(qarr) < q for x in tots))
