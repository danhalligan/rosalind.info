# Reconstruct a String from its Burrows-Wheeler Transform

from collections import defaultdict


# Yields each character with the occurrence number
def index_seq(seq):
    d = defaultdict(int)
    for c in seq:
        yield c, d[c]
        d[c] += 1


def bwtj(seq):
    first = list(index_seq(sorted(seq)))
    last = list(index_seq(seq))
    curr = ("$", 0)
    res = ""
    for _ in range(len(seq)):
        curr = first[last.index(curr)]
        res += curr[0]
    return res


def main(file):
    print(bwtj(open(file).read().rstrip()))
