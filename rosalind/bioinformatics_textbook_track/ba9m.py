# Implement BetterBWMatching

from collections import defaultdict
from copy import copy
from itertools import accumulate


# First position of each symbol in the first column. Note the first column
# is always sorted lexicographically.
def first_occurrence(seq):
    letters = sorted(set(seq))
    counts = [0] + list(accumulate(seq.count(x) for x in letters))
    return dict(zip(letters, counts))


def count_symbols(seq):
    count = []
    count.append(defaultdict(int))
    for i, s in enumerate(seq):
        count.append(copy(count[i]))
        count[i + 1][s] += 1
    return count


# Unlike the book, this returns an array of match indexes (not e.g. length of
# that array, or the top/bottom indexes).
def better_bwmatching(FirstOccurrence, LastColumn, Pattern, Countsymbol):
    top, bottom = 0, len(LastColumn) - 1
    while top <= bottom:
        if Pattern:
            Pattern, symbol = Pattern[:-1], Pattern[-1]
            if symbol in LastColumn[top : bottom + 1]:
                top = FirstOccurrence[symbol] + Countsymbol[top][symbol]
                bottom = FirstOccurrence[symbol] + Countsymbol[bottom + 1][symbol] - 1
            else:
                return []
        else:
            return list(range(top, bottom + 1))


def main(file):
    text, patterns = open(file).read().splitlines()
    patterns = patterns.split()
    fo = first_occurrence(text)
    cs = count_symbols(text)
    print(*[len(better_bwmatching(fo, text, p, cs)) for p in patterns])
