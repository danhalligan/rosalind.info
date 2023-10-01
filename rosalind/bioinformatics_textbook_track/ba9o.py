# Find All Approximate Occurrences of a Collection of Patterns in a String

from .ba9q import partial_suffix_array
from .ba9i import bwt
from .ba9m import first_occurrence, count_symbols
from .ba9n import find_location


def update_ptrs(ptrs, fo, cs, sym):
    top, bottom = ptrs
    top = fo[sym] + cs[top][sym]
    bottom = fo[sym] + cs[bottom + 1][sym] - 1
    return (top, bottom)


# Unlike the book, this returns an array of match indexes (not e.g. length of
# that array, or the top/bottom indexes).
def approxbwm(fo, seq, pattern, cs, ptrs, m_count, n):
    if not pattern:
        # no pattern left, so everything between ptrs is a match
        return list(range(ptrs[0], ptrs[1] + 1))
    matches = []
    pattern, sym = pattern[:-1], pattern[-1]
    symbols = [seq[i] for i in range(ptrs[0], ptrs[1] + 1)]
    if sym in symbols:
        # we have a match -- match rest of pattern
        nptrs = update_ptrs(ptrs, fo, cs, sym)
        matches += approxbwm(fo, seq, pattern, cs, nptrs, m_count, n)
    if m_count < n:
        for mm in set(symbols) - set(sym):
            nptrs = update_ptrs(ptrs, fo, cs, mm)
            matches += approxbwm(fo, seq, pattern, cs, nptrs, m_count + 1, n)
    return matches


def approx_match_positions(seq, patterns, mismatches, k=10):
    psa = dict(partial_suffix_array(seq + "$", k))
    seq = bwt(seq + "$")
    fo = first_occurrence(seq)
    cs = count_symbols(seq)
    for pattern in patterns:
        ptrs = (0, len(seq) - 1)
        for match in approxbwm(fo, seq, pattern, cs, ptrs, 0, mismatches):
            yield find_location(match, psa, seq, fo, cs)


def main(file):
    seq, patterns, mismatches = open(file).read().splitlines()
    patterns = patterns.split()
    mismatches = int(mismatches)
    print(*sorted(approx_match_positions(seq, patterns, mismatches)))
