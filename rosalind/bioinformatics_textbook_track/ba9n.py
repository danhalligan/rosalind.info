# Find All Occurrences of a Collection of Patterns in a String

from .ba9q import partial_suffix_array
from .ba9i import bwt
from .ba9m import better_bwmatching, first_occurrence, count_symbols
from .ba9j import index_seq


# Locate match within text using partial suffix array
def find_location(row, psa, first_indexed, last_indexed):
    steps = 0
    while row not in psa:
        row = first_indexed.index(last_indexed[row])
        steps += 1
    return steps + psa[row]


# Use better_bwmatching to obtain match indices, then we use the
# (partial) suffix array to recover the starting position for each
def match_positions(text, patterns, k=10):
    psa = dict(partial_suffix_array(text + "$", k))
    last_column = bwt(text + "$")
    first_indexed = list(index_seq(sorted(last_column)))
    last_indexed = list(index_seq(last_column))
    fo = first_occurrence(last_column)
    cs = count_symbols(last_column)
    for pattern in patterns:
        for match in better_bwmatching(fo, last_column, pattern, cs):
            yield find_location(match, psa, first_indexed, last_indexed)


def main(file):
    text, *patterns = open(file).read().splitlines()
    matches = sorted(match_positions(text, patterns))
    print(*matches)
