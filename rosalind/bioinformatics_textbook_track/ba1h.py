# Find All Approximate Occurrences of a Pattern in a String

from .ba1a import substrings
from .ba1g import hamming


def find_approx_matches(pattern, text, d):
    for i, x in enumerate(substrings(text, len(pattern))):
        if hamming(x, pattern) <= d:
            yield i


def main(file):
    pattern, text, d = open(file).read().splitlines()
    print(*find_approx_matches(pattern, text, int(d)))
