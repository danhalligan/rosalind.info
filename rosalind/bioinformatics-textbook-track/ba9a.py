# Construct a Trie from a Collection of Patterns

from itertools import groupby


def trie(seqs):
    graph = {}
    if len(seqs):
        for base, nseqs in groupby(sorted(seqs), key=lambda s: s[0]):
            graph[base] = trie([seq[1:] for seq in nseqs if len(seq) > 1])
    return graph


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
