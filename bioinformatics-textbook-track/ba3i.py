# Find a k-Universal Circular String

from itertools import product
from .ba3b import genome_path
from .ba3d import dbru
from .ba3f import eulerian_cycle


def k_universal_binary_string(k):
    kmers = ["".join(x) for x in product(["0", "1"], repeat=k)]
    cycle = eulerian_cycle(dbru(kmers))
    return genome_path(cycle[: -(k - 1)])


def main(file):
    k = int(open(file).read().splitlines()[0])
    print(k_universal_binary_string(k))
