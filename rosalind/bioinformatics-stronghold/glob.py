# Global Alignment with Scoring Matrix

from .helpers import Parser
from rosalind.helpers import blosum62


def glob(s1, s2, score, penalty):
    """Global Alignment with Scoring Matrix"""
    m = {}
    for j in range(len(s2) + 1):
        m[j, 0] = penalty * j
    for i in range(len(s1) + 1):
        m[0, i] = penalty * i

    for j in range(len(s2)):
        for i in range(len(s1)):
            pos = [(j + 1, i), (j, i), (j, i + 1)]
            cost = [penalty, score[s1[i]][s2[j]], penalty]
            scores = [m[pos[x]] + cost[x] for x in range(3)]
            m[j + 1, i + 1] = max(scores)

    return m[len(s2), len(s1)]


def main(file):
    seqs = Parser(file).seqs()
    print(glob(seqs[0], seqs[1], blosum62(), -5))
