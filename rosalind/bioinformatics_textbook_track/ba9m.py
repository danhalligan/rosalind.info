# Implement BetterBWMatching

from collections import defaultdict
from copy import copy


def first_occurrence(seq):
    first = sorted(seq)
    return dict((x, first.index(x)) for x in set(seq))


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
    top = 0
    bottom = len(LastColumn) - 1
    while top <= bottom:
        if Pattern:
            Pattern, symbol = Pattern[:-1], Pattern[-1]
            x = [LastColumn[i] for i in range(top, bottom + 1)]
            if symbol in x:
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
