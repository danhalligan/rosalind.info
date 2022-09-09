# Maximizing the Gap Symbols of an Optimal Alignment

from collections import defaultdict
from .helpers import Parser


def mgap(s1, s2):
    """Maximizing the Gap Symbols of an Optimal Alignment"""
    m = defaultdict(lambda: 0)

    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i] + 1
            else:
                m[j + 1, i + 1] = max([m[j + 1, i], m[j, i + 1]])

    return len(s1) + len(s2) - 2 * m[len(s2), len(s1)]


def main(file):
    seqs = Parser(file).seqs()
    print(mgap(seqs[0], seqs[1]))
