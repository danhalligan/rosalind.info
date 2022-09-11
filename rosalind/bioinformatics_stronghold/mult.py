# Multiple Alignment

from .helpers import Parser
from itertools import product, combinations


def valid_coord(x, pos):
    return all([i >= 0 for i in prev(pos, x)])


def prev(pos, ptr):
    return tuple([p + d for p, d in zip(pos, ptr)])


def insertions(seqs, pos, ptr):
    return [seq[pos[i] - 1] if ptr[i] == -1 else "-" for i, seq in enumerate(seqs)]


# score is obtained as sum over all possible pairs
def score(seqs, pos, ptr):
    a = insertions(seqs, pos, ptr)
    return sum(0 if a == b else -1 for a, b in combinations(a, 2))


def moves(n):
    return list(product([0, -1], repeat=n))[1:]


def mult(seqs):
    m, p = {}, {}
    m[0, 0, 0, 0] = 0
    ranges = [range(0, len(s) + 1) for s in seqs]
    for pos in product(*ranges):
        ptrs = list(filter(lambda x: valid_coord(x, pos), moves(4)))
        if not len(ptrs):
            continue
        sc = [m[prev(pos, x)] + score(seqs, pos, x) for x in ptrs]
        m[pos] = max(sc)
        p[pos] = ptrs[sc.index(max(sc))]

    # traceback to recover alignment
    tot = m[pos]
    aln = ["", "", "", ""]
    while any([x > 0 for x in pos]):
        ptr = p[pos]
        for i, v in enumerate(insertions(seqs, pos, ptr)):
            aln[i] += v
        pos = prev(pos, ptr)

    return tot, *[a[::-1] for a in aln]


def main(file):
    seqs = Parser(file).seqs()
    print(*mult(seqs), sep="\n")
