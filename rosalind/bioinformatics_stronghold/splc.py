# RNA Splicing

from functools import reduce
from .helpers import Parser, Dna


def main(file):
    def trim(gene, intron):
        s = gene.find(intron)
        return gene[:s] + gene[s + len(intron) :]

    seqs = Parser(file).seqs()
    print(Dna(reduce(trim, seqs)).translate())
