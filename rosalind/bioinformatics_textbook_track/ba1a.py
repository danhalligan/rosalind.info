# Compute the Number of Times a Pattern Appears in a Text


def substrings(text, size):
    for i in range(len(text) - size + 1):
        yield text[i : i + size]


def pattern_count(text, pattern):
    return sum(pattern == x for x in substrings(text, len(pattern)))


def main(file):
    text, pattern = open(file).read().splitlines()
    print(pattern_count(text, pattern))
