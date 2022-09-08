# Find a Median String

from .ba1a import substrings
from .ba1i import generate_kmers
from .ba1g import hamming
import math


def minimum_distance(x, text):
    return min(hamming(w, x) for w in substrings(text, len(x)))


def distance_between_pattern_and_strings(pattern, dna):
    return sum(minimum_distance(pattern, x) for x in dna)


def median_string(dna, k):
    dist = math.inf
    for kmer in generate_kmers(k):
        x = distance_between_pattern_and_strings(kmer, dna)
        if dist > x:
            dist, median = x, kmer
    return median


def main(file):
    k, *dna = open(file).read().splitlines()
    print(median_string(dna, int(k)))
