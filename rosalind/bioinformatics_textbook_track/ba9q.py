# Construct the Partial Suffix Array of a String

from .ba9g import suffix_array


def partial_suffix_array(seq, k):
    return [(i, x) for i, x in enumerate(suffix_array(seq)) if x % k == 0]


def main(file):
    seq, k = open(file).read().splitlines()
    for x in partial_suffix_array(seq, int(k)):
        print(*x, sep=",")
