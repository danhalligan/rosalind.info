# Local Alignment with Scoring Matrix

from .helpers import Parser
from rosalind.helpers import pam250


def loca(s1, s2, score, penalty):
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = 0
        p[j, 0] = "↑"
    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = "←"

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            opt = [
                m[j, i] + score[s1[i]][s2[j]],
                m[j, i + 1] + penalty,
                m[j + 1, i] + penalty,
                0,
            ]
            m[new] = max(opt)
            p[new] = ["↖", "↑", "←", "↖"][opt.index(max(opt))]

    max_score = max(x for x in m.values())
    j, i = [k for k, v in m.items() if v == max_score][0]
    a1, a2 = "", ""
    while i > 0 or j > 0:
        if m[j, i] == 0:
            break
        if p[j, i] == "↖":
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j, i = j - 1, i - 1
        elif p[j, i] == "←":
            a1 += s1[i - 1]
            i = i - 1
        elif p[j, i] == "↑":
            a2 += s2[j - 1]
            j = j - 1

    return {"dist": max_score, "a1": a1[::-1], "a2": a2[::-1]}


def main(file):
    seqs = Parser(file).seqs()
    res = loca(seqs[0], seqs[1], pam250(), -5)
    print(res["dist"])
    print(res["a1"])
    print(res["a2"])
