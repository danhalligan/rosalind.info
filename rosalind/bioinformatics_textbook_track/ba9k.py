# Generate the Last-to-First Mapping of a String

from .ba9j import index_seq


def last2first(seq, i):
    first = list(index_seq(sorted(seq)))
    last = list(index_seq(seq))
    return first.index(last[i])


def main(file):
    seq, i = open(file).read().splitlines()
    print(last2first(seq, int(i)))
