# Generate the d-Neighborhood of a String

from .ba1g import hamming


def immediate_neighbors(seq):
    for i in range(len(seq)):
        for s in ["A", "T", "G", "C"]:
            if s != seq[i]:
                yield seq[:i] + s + seq[i + 1 :]


def neighbors(seq, d):
    bases = ["A", "T", "G", "C"]
    if d == 0:
        return set([seq])
    if len(seq) == 1:
        return set(bases)
    n = set()
    for x in neighbors(seq[1:], d):
        if hamming(seq[1:], x) < d:
            for b in bases:
                n.add(b + x)
        else:
            n.add(seq[0] + x)
    return n


def main(file):
    seq, d = open(file).read().splitlines()
    print(*sorted(neighbors(seq, int(d))), sep="\n")
