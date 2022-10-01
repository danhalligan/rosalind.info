# Isolating Symbols in Alignments

import numpy as np
from .helpers import Parser


def score(a, b):
    return 1 if a == b else -1


def glob(s1, s2):
    m = np.zeros((len(s1) + 1, len(s2) + 1), int)
    for i in range(len(s1) + 1):
        m[i][0] = -1 * i
    for j in range(len(s2) + 1):
        m[0][j] = -1 * j

    for i in range(len(s1)):
        for j in range(len(s2)):
            m[i + 1][j + 1] = max(
                m[i][j + 1] - 1,
                m[i][j] + score(s1[i], s2[j]),
                m[i + 1][j] - 1,
            )

    return m


def osym(s1, s2):
    pm = glob(s1, s2)
    sm = glob(s1[::-1], s2[::-1])
    scores = [
        pm[i][j] + score(s1[i], s2[j]) + sm[len(s1) - 1 - i][len(s2) - 1 - j]
        for i in range(len(s1))
        for j in range(len(s2))
    ]
    return max(scores), sum(scores)


def main(file):
    print(*osym(*Parser(file).seqs()), sep="\n")
