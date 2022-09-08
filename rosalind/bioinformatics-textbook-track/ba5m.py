# Find a Highest-Scoring Multiple Sequence Alignment

from itertools import product


# Does this pointer/score refer to a non-negative cell in our matrix?
def valid_coord(x, pos):
    return all([i >= 0 for i in prev(pos, x)])


# get previous position given a pointer
def prev(pos, ptr):
    return tuple([p + d for p, d in zip(pos, ptr)])


# given sequences, a position and a pointer, calculate score.
# This represents a scoring function in which the score of an alignment column
# is 1 if all three symbols are identical and 0 otherwise.
def score(seqs, pos, ptr):
    if ptr == (-1, -1, -1):
        bases = [seqs[i][j] for i, j in enumerate(prev(pos, ptr))]
        return all(x == bases[0] for x in bases)
    else:
        return 0


# generate possible previous cells relative to current cell (pointers)
def moves(n):
    return list(product([0, -1], repeat=n))[1:]


def multiple_alignment(seqs):
    m, p = {}, {}
    m[0, 0, 0] = 0
    ranges = [range(0, len(s) + 1) for s in seqs]
    for pos in product(*ranges):
        ptrs = list(filter(lambda x: valid_coord(x, pos), moves(3)))
        if not len(ptrs):
            continue
        sc = [m[prev(pos, x)] + score(seqs, pos, x) for x in ptrs]
        m[pos] = max(sc)
        p[pos] = ptrs[sc.index(max(sc))]

    # traceback to recover alignment
    tot = m[pos]
    aln = ["", "", ""]
    while any([x > 0 for x in pos]):
        ptr = p[pos]
        for i, seq in enumerate(seqs):
            aln[i] += seq[pos[i] - 1] if ptr[i] == -1 else "-"
        pos = prev(pos, ptr)

    return tot, aln[0][::-1], aln[1][::-1], aln[2][::-1]


def main(file):
    seqs = open(file).read().splitlines()
    print(*multiple_alignment(seqs), sep="\n")
