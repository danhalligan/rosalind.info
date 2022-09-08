# Generate the k-mer Composition of a String

from .ba1a import substrings


def main(file):
    k, seq = open(file).read().splitlines()
    print(*substrings(seq, int(k)), sep="\n")
