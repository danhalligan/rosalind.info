# Find Patterns Forming Clumps in a String

from collections import defaultdict
from .ba1a import substrings


# find (0-based) positions of k-length kmers within text
def find_kmers(text, k):
    x = defaultdict(list)
    for i, t in enumerate(substrings(text, k)):
        x[t] += [i]
    return x


# Do a given array of kmers at positions p form a clump of t within L
# given that each is k-long.
def has_clump(p, L, t, k):
    for i in range(len(p) - t + 1):
        if (p[i + t - 1] + k - p[i]) <= L:
            return True
    return False


def main(file):
    seq, ints = open(file).read().splitlines()
    k, L, t = list(map(int, ints.split()))
    kmers = find_kmers(seq, k)
    print(*[x for x in kmers if has_clump(kmers[x], L, t, k)])
