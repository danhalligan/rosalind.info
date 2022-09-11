# Find the Most Frequent Words in a String

from collections import defaultdict
from .ba1a import substrings


def count_kmers(text, k):
    d = defaultdict(int)
    for x in substrings(text, k):
        d[x] += 1
    return d


def most_frequent(d):
    return [k for k in d if d[k] == max(d.values())]


def main(file):
    text, k = open(file).read().splitlines()
    print(*most_frequent(count_kmers(text, int(k))))
