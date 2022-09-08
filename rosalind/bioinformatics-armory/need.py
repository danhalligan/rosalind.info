# Pairwise Global Alignment

from Bio import SeqIO
from .frmt import get_records
import sys


# Here we just write fasta format to STDOUT, paste these into the Needle
# web interface.
def main(file):
    ids = open(file).read().split()
    records = get_records(ids)
    for rcrd in records:
        SeqIO.write(rcrd, sys.stdout, "fasta")
