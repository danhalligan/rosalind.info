# Finding Disjoint Motifs in a Gene

from .helpers import Parser


# Generate possible forward moves relative to current position
# if we've started including a pattern, we have to increment sequence "s"
def moves(i, j, k):
    if j > 0 or k > 0:
        return [(1, 0, 1), (1, 1, 0)]
    else:
        return [(1, 0, 0), (1, 0, 1), (1, 1, 0)]


# Calculate the next position given current position and the increment
def npos(pos, move):
    return tuple(a + b for a, b in zip(pos, move))


# Is the proposed position valid, given the increment made and the sequences
# Valid moves are within the sequence lengths and with matching bases
def valid(pos, move, seqs):
    for p, x in zip(pos, seqs):
        if p > len(x):
            return False
    bases = [x[p - 1] for x, p, m in zip(seqs, pos, move) if m == 1 and p > 0]
    return len(set(bases)) == 1


# Queue based iterative approach. Pop a position from the queue and add all
# new valid positions.
# If we reach the end of both our patterns, they can be interleaved.
def itwv(s, t, u):
    pos = (0, 0, 0)
    q = [pos]
    seqs = [s, t, u]

    while q:
        pos = q.pop(0)
        for move in moves(*pos):
            np = npos(pos, move)
            if valid(np, move, seqs):
                if np[1] == len(t) and np[2] == len(u):
                    return 1
                q.append(np)
    return 0


def main(file):
    seq, *patterns = Parser(file).lines()
    res = [[itwv(seq, p1, p2) for p1 in patterns] for p2 in patterns]
    for x in res:
        print(*x)
