# Find Frequent Words with Mismatches and Reverse Complements

from .ba1c import revcomp
from .ba1g import hamming
from .ba1b import count_kmers, most_frequent
from .ba1i import generate_kmers


def count_hamming_revcomp_kmers(kmers, d, k):
    for x in generate_kmers(k):
        count = sum(kmers[y] for y in kmers if hamming(x, y) <= d)
        count += sum(kmers[y] for y in kmers if hamming(revcomp(x), y) <= d)
        if count > 0:
            yield [x, count]


def main(file):
    seq, ints = open(file).read().splitlines()
    k, d = list(map(int, ints.split()))
    kmers = count_kmers(seq, k)
    hkmers = dict(count_hamming_revcomp_kmers(kmers, d, k))
    print(*most_frequent(hkmers))
