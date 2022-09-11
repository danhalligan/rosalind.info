# Data Formats

from .gbk import entrez_email
from Bio import Entrez, SeqIO
import sys


def get_records(ids):
    Entrez.email = entrez_email()
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
    return list(SeqIO.parse(handle, "fasta"))


def frmt(ids):
    records = get_records(ids)
    lengths = [len(x) for x in records]
    return records[lengths.index(min(lengths))]


def main(file):
    ids = open(file).read().split()
    rcrd = frmt(ids)
    SeqIO.write(rcrd, sys.stdout, "fasta")
