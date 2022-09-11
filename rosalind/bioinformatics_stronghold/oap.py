# Overlap Alignment

from .helpers import Parser
import numpy as np


# Since the real data is quite big for this challenge, we'll use
# numpy integer arrays for best performance.
def oap(s1, s2, penalty=-2):
    score = np.empty((len(s2) + 1, len(s1) + 1), dtype=int)
    ptr = np.empty((len(s2) + 1, len(s1) + 1), dtype=int)

    for j in range(len(s2) + 1):
        score[j][0] = j * penalty
        ptr[j][0] = 1
    for i in range(len(s1) + 1):
        score[0][i] = 0
        ptr[0][i] = 2

    score[0][0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            opt = [
                score[j][i] + (1 if s1[i] == s2[j] else penalty),
                score[j][i + 1] + penalty,
                score[j + 1][i] + penalty,
            ]
            best = max(opt)
            score[j + 1][i + 1] = best
            ptr[j + 1][i + 1] = opt.index(best)

    sc = [score[j][len(s1)] for j in range(len(s2) + 1)]
    max_score = max(sc)
    j = [j for j, s in enumerate(sc) if s == max_score][-1]
    i = len(s1)
    a1, a2 = "", ""
    while i > 0 and j > 0:
        if ptr[j][i] == 0:
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j, i = j - 1, i - 1
        elif ptr[j][i] == 1:
            a1 += "-"
            a2 += s2[j - 1]
            j = j - 1
        elif ptr[j][i] == 2:
            a1 += s1[i - 1]
            a2 += "-"
            i = i - 1

    return max_score, a1[::-1], a2[::-1]


def main(file):
    s1, s2 = Parser(file).seqs()
    print(*oap(s1, s2, -2), sep="\n")
