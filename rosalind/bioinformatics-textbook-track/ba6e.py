# Find All Shared k-mers of a Pair of Strings

import re
from .ba1c import revcomp


def overlapping_matches(pattern, seq):
    return re.finditer(rf"(?=({pattern}))", seq)


def ba6e(k, s1, s2):
    for i in range(len(s1) - k + 1):
        seq = s1[i : (i + k)]
        for s in [seq, revcomp(seq)]:
            match = list(overlapping_matches(s, s2))
            for m in match:
                yield (i, m.start())


def main(file):
    k, s1, s2 = open(file).read().splitlines()
    for match in ba6e(int(k), s1, s2):
        print(match)
