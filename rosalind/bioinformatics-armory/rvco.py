# Complementing a Strand of DNA

from Bio import SeqIO


def main(file):
    count = 0
    for rcrd in SeqIO.parse(file, "fasta"):
        count += rcrd.seq == rcrd.reverse_complement().seq
    print(count)
