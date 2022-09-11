# Construct a Trie from a Collection of Patterns

from rosalind.bioinformatics_stronghold.trie import trie


def as_adjacency(graph, nodes=[]):
    node = len(nodes)
    for edge in sorted(graph):
        i = len(nodes)
        nodes.append((node, i + 1, edge))
        nodes + as_adjacency(graph[edge], nodes)
    return nodes


def main(file):
    seqs = open(file).read().splitlines()
    for edge in as_adjacency(trie(seqs)):
        print(f"{edge[0]}->{edge[1]}:{edge[2]}")
