# Creating a Distance Matrix

from .helpers import Parser
from .hamm import hamm


def pdst(seqs):
    """Creating a Distance Matrix"""
    n = len(seqs)
    d = [[0.0 for x in range(n)] for y in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d[i][j] = d[j][i] = hamm(seqs[i], seqs[j]) / len(seqs[i])
    return d


def main(file):
    for r in pdst(Parser(file).seqs()):
        print(*[f"{x:.5f}" for x in r])
