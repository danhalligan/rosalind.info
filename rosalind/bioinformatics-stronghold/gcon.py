# Global Alignment with Constant Gap Penalty

from .helpers import Parser
from rosalind.helpers import blosum62


def gcon(s1, s2, score, penalty):
    """Global Alignment with Constant Gap Penalty"""
    # See Biological Sequence Analysis page 29
    # Its a simplification of the more general affine gap penalty
    # with an extension penalty of 0.
    # We now have to keep track of three matrices m, x and y

    m, x, y = {}, {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = penalty
        y[j, 0] = -99
    for i in range(len(s1) + 1):
        m[0, i] = penalty
        x[0, i] = -99

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            x[new] = max([m[j, i + 1] + penalty, x[j, i + 1]])
            y[new] = max([m[j + 1, i] + penalty, y[j + 1, i]])
            m[new] = max([m[j, i] + score[s1[i]][s2[j]], x[new], y[new]])

    return m[len(s2), len(s1)]


def main(file):
    seqs = Parser(file).seqs()
    print(gcon(seqs[0], seqs[1], blosum62(), -5))
