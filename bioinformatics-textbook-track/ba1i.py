# Find the Most Frequent Words with Mismatches in a String

from .ba1g import hamming
from .ba1b import count_kmers, most_frequent
from itertools import product

# Note, the best kmer might not be observed in our sequence. The simplistic
# method here simply checks all possible kmers (which is ~17M for k = 12)


def generate_kmers(k):
    return ("".join(x) for x in product(["A", "C", "G", "T"], repeat=k))


def count_hamming_kmers(kmers, d, k):
    for x in generate_kmers(k):
        count = sum(kmers[y] for y in kmers if hamming(x, y) <= d)
        if count > 0:
            yield [x, count]


def main(file):
    seq, ints = open(file).read().splitlines()
    k, d = list(map(int, ints.split()))
    kmers = count_kmers(seq, k)
    hkmers = dict(count_hamming_kmers(kmers, d, k))
    print(*most_frequent(hkmers))
