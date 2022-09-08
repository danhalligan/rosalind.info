# Implement MotifEnumeration

from .ba1n import neighbors
from .ba1g import hamming
from .ba1a import substrings


def all_kmers(dna, k):
    return set(x for text in dna for x in substrings(text, k))


def contains_approx_match(pattern, text, d):
    for x in substrings(text, len(pattern)):
        if hamming(x, pattern) <= d:
            return True
    return False


def motif_enumeration(dna, k, d):
    patterns = set()
    for kmer in all_kmers(dna, k):
        for hkmer in neighbors(kmer, d):
            if all(contains_approx_match(hkmer, x, d) for x in dna):
                patterns.add(hkmer)
    return patterns


def main(file):
    ints, *dna = open(file).read().splitlines()
    k, d = map(int, ints.split())
    print(*sorted(motif_enumeration(dna, k, d)))
