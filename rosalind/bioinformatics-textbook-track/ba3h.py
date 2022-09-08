# Reconstruct a String from its k-mer Composition

from .ba3b import genome_path
from .ba3d import dbru
from .ba3g import eulerian_path


def main(file):
    k, *dna = open(file).read().splitlines()
    print(genome_path(eulerian_path(dbru(dna))))
