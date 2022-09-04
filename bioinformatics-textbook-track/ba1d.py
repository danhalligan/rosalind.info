# Find All Occurrences of a Pattern in a String

from .ba1a import substrings


def pattern_find(text, pattern):
    for i, x in enumerate(substrings(text, len(pattern))):
        if x == pattern:
            yield i


def main(file):
    pattern, text = open(file).read().splitlines()
    print(*pattern_find(text, pattern))
