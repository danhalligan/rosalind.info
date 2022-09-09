# Encoding Suffix Trees

from functools import cache
from os.path import commonprefix
from .helpers import Parser


def get_edges(graph):
    for k in graph.keys():
        yield k
        yield from get_edges(graph[k])


@cache
def suffix_tree(seq, starts):
    graph = {}
    bases = sorted(set([seq[start] for start in starts]))
    for base in bases:
        matching = [start for start in starts if seq[start] == base]
        seqs = [seq[s:] for s in matching]
        prefix = commonprefix(seqs)
        size = len(prefix)
        new_starts = [start + size for start in matching if start + size < len(seq)]
        graph[prefix] = suffix_tree(seq, tuple(new_starts))
    return graph


def suff(seq):
    return suffix_tree(seq, tuple(range(len(seq))))


def main(file):
    seq = Parser(file).line()
    tree = suff(seq)
    print(*list(get_edges(tree)), sep="\n")
