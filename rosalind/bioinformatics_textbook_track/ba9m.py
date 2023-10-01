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


def better_bwmatching(LastColumn, Pattern):
    FirstOccurrence = first_occurrence(LastColumn)
    Countsymbol = count_symbols(LastColumn)
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
                return None
        else:
            return top, bottom


def count_matches(LastColumn, Pattern):
    m = better_bwmatching(LastColumn, Pattern)
    if m:
        return m[1] - m[0] + 1
    else:
        return 0


def main(file):
    text, patterns = open(file).read().splitlines()
    patterns = patterns.split()
    print(*[count_matches(text, pattern) for pattern in patterns])
