# Genome Assembly with Perfect Coverage

from .helpers import Parser
from .dbru import dbru


def find_cycle(graph):
    key = sorted(list(graph.keys()))[0]
    visited = set()
    cycle = []
    while key not in visited:
        cycle += [key]
        visited.add(key)
        key = graph[key]
    return cycle


def join_cycle(chain):
    return "".join(x[0] for x in chain)


def pcov(seqs):
    """Genome Assembly with Perfect Coverage"""
    x = find_cycle(dict(dbru(seqs)))
    return join_cycle(sorted(x))


def main(file):
    print(pcov(Parser(file).lines()))
