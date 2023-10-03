# Find All Approximate Occurrences of a Collection of Patterns in a String

from .ba9q import partial_suffix_array
from .ba9i import bwt
from .ba9m import first_occurrence, count_symbols
from .ba9n import find_location


class BWMatch:
    def __init__(self, seq, k=10):
        self.psa = dict(partial_suffix_array(seq + "$", k))
        self.seq = bwt(seq + "$")
        self.fo = first_occurrence(self.seq)
        self.cs = count_symbols(self.seq)

    def update(self, ptrs, x):
        t, b = ptrs
        return (self.fo[x] + self.cs[t][x], self.fo[x] + self.cs[b + 1][x] - 1)

    # pattern: the pattern to match
    # ptrs: top and bottom pointers within last column
    # mi: counter for number of mismatches
    # mn: total number of acceptable mismatches
    def bwm(self, pattern, ptrs, mi):
        if not pattern:
            return list(range(ptrs[0], ptrs[1] + 1))
        matches = []
        pattern, sym = pattern[:-1], pattern[-1]
        if sym in self.seq[ptrs[0] : ptrs[1] + 1]:
            matches += self.bwm(pattern, self.update(ptrs, sym), mi)
        if mi < self.mn:
            for mm in ["A", "C", "G", "T"]:
                if mm != sym:
                    matches += self.bwm(pattern, self.update(ptrs, mm), mi + 1)
        return matches

    def match_patterns(self, patterns, mn):
        self.mn = mn
        for pattern in patterns:
            for match in self.bwm(pattern, (0, len(self.seq) - 1), 0):
                yield find_location(match, self.psa, self.seq, self.fo, self.cs)


def main(file):
    seq, patterns, mismatches = open(file).read().splitlines()
    patterns = patterns.split()
    mismatches = int(mismatches)
    matcher = BWMatch(seq)
    print(*sorted(matcher.match_patterns(patterns, mismatches)))
