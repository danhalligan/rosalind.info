# Construct a String Spelled by a Gapped Genome Path

from .ba3j import string_from_paired_composition


def main(file):
    ints, *pairs = open(file).read().splitlines()
    k, d = map(int, ints.split())
    pairs = [x.split("|") for x in pairs]
    print(string_from_paired_composition(pairs, k, d))
