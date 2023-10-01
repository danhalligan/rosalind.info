# Find All Occurrences of a Collection of Patterns in a String

from .ba9q import partial_suffix_array
from .ba9i import bwt
from .ba9m import better_bwmatching, first_occurrence, count_symbols


# Locate match within text using partial suffix array
# We don't use our `index_seq` function here (which stores the occurrence of
# each character in a sequence) since it kind of defeats the
# memory saving from using a partial suffix array.
# We can use the fact that the first column is lexicographically sorted to
# easily find the nth occurrence of a character within the first column
# (its just the position of the 1st occurrence + n)!
def find_location(row, psa, last_column, first):
    steps = 0
    while row not in psa:
        predecessor = last_column[row]
        occurrence = last_column[: row + 1].count(predecessor)
        row = first[predecessor] + occurrence - 1
        steps += 1
    return steps + psa[row]


# Use better_bwmatching to obtain match indices, then we use the
# (partial) suffix array to recover the starting position for each
def match_positions(text, patterns, k=10):
    psa = dict(partial_suffix_array(text + "$", k))
    last_column = bwt(text + "$")
    fol = first_occurrence(last_column)
    FirstColumn = sorted(last_column)
    fos = first_occurrence(FirstColumn)
    cs = count_symbols(last_column)
    for pattern in patterns:
        for match in better_bwmatching(fol, last_column, pattern, cs):
            yield find_location(match, psa, last_column, fos)


def main(file):
    text, *patterns = open(file).read().splitlines()
    matches = sorted(match_positions(text, patterns))
    print(*matches)
