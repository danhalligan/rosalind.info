# Implement 2-BreakOnGenome

from .ba6a import read_perm, format_perm
from .ba6h import colored_edges
from .ba6j import two_break_on_genome_graph
from .ba6i import graph2genome


def ints(x):
    return list(map(int, x.split(", ")))


# we don't exactly follow the pseudocode from rosalind here...
def two_break_on_genome(P, i, ip, j, jp):
    genome_graph = colored_edges([P])
    genome_graph = two_break_on_genome_graph(genome_graph, i, ip, j, jp)
    return graph2genome(genome_graph)


def main(file):
    genome, ind = open(file).read().splitlines()
    genome = read_perm(genome)
    new = two_break_on_genome(genome, *ints(ind))
    print(*[format_perm(x) for x in new])
