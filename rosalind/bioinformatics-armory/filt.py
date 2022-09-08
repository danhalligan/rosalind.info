# Read Filtration by Quality

from Bio import SeqIO


def ints(x):
    return list(map(int, x.split()))


def good_quality(q, p, seq):
    quals = seq.letter_annotations["phred_quality"]
    return sum(x >= q for x in quals) / len(quals) >= p / 100


def main(file):
    handle = open(file)
    q, p = ints(next(handle))
    print(sum(good_quality(q, p, x) for x in SeqIO.parse(handle, "fastq")))
