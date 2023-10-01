# Find All Occurrences of a Collection of Patterns in a String

from .ba9q import partial_suffix_array
from .ba9i import bwt
from .ba9m import better_bwmatching


# Given a character that is known to appear in a specific position in a column,
# determine whether it is the first, 2nd, ... occurrence
def getN(ch, row, column):
    n, i = 0, 0
    while i <= row:
        if ch == column[i]:
            n += 1
        i += 1
    return n


# Find a specific occurrence of character in column
def get_char(ch, column, occurrence):
    count = 0
    for pos, x in enumerate(column):
        if ch == x:
            count += 1
            if count == occurrence:
                return pos


# Locate match within text using partial suffix array
def find_position(row, psa, LastColumn):
    psa = dict(psa)
    FirstColumn = sorted(LastColumn)
    steps = 0
    while row not in psa:
        predecessor = LastColumn[row]
        occurrence = getN(predecessor, row, LastColumn)
        row = get_char(predecessor, FirstColumn, occurrence)
        steps += 1
    return steps + psa[row]


def match_positions(text, patterns, k=10):
    psa = partial_suffix_array(text + "$", k)
    LastColumn = bwt(text + "$")
    for Pattern in patterns:
        matches = better_bwmatching(LastColumn, Pattern)
        if matches:
            top, bottom = matches
            yield [find_position(i, psa, LastColumn) for i in range(top, bottom + 1)]


def main(file):
    text, *patterns = open(file).read().splitlines()
    matches = sorted(m for ms in match_positions(text, patterns) for m in ms)
    print(*matches)
