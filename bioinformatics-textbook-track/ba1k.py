# Generate the Frequency Array of a String

from .ba1a import pattern_count
from .ba1i import generate_kmers


def frequency_array(seq, k):
    return (pattern_count(seq, x) for x in generate_kmers(k))


def main(file):
    seq, k = open(file).read().splitlines()
    print(*frequency_array(seq, int(k)))
