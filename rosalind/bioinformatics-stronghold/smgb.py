# Semiglobal Alignment

from .helpers import Parser
import numpy as np

# As with OAP, the real data is quite big, we'll use numpy integer arrays for
# improved performance over other solutions. Even then, this could take
# minutes to run on the "real" dataset.


def smgb(s1, s2):
    m = np.empty((len(s2) + 1, len(s1) + 1), dtype=int)
    p = np.empty((len(s2) + 1, len(s1) + 1), dtype=int)
    for j in range(len(s2) + 1):
        m[j][0] = 0
        p[j][0] = 1
    for i in range(len(s1) + 1):
        m[0][i] = 0
        p[0][i] = 2

    for i in range(len(s1)):
        for j in range(len(s2)):
            opt = [
                m[j][i] + (1 if s1[i] == s2[j] else -1),
                m[j][i + 1] + (-1 if 0 < i + 1 < len(s1) else 0),
                m[j + 1][i] + (-1 if 0 < j + 1 < len(s2) else 0),
            ]
            m[j + 1][i + 1] = max(opt)
            p[j + 1][i + 1] = opt.index(max(opt))

    i, j = len(s1), len(s2)
    max_score = m[j][i]
    a1, a2 = "", ""
    while i > 0 or j > 0:
        if p[j][i] == 0:
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j, i = j - 1, i - 1
        elif p[j][i] == 1:
            a2 += s2[j - 1]
            a1 += "-"
            j = j - 1
        elif p[j][i] == 2:
            a1 += s1[i - 1]
            a2 += "-"
            i = i - 1

    return {"dist": max_score, "a1": a1[::-1], "a2": a2[::-1]}


def main(file):
    seqs = Parser(file).seqs()
    res = smgb(seqs[0], seqs[1])
    print(res["dist"])
    print(res["a1"])
    print(res["a2"])
