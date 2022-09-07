# Implement 2-BreakOnGenomeGraph

from .ba6i import parse_edge_string
from copy import copy


def ints(x):
    return list(map(int, x.split(", ")))


def rm_edge(graph, edge):
    if edge[0] in graph:
        graph.pop(edge[0])
    else:
        graph.pop(edge[1])


def two_break_on_genome_graph(graph, i, ip, j, jp):
    new = copy(graph)
    for edge in [(i, ip), (j, jp)]:
        rm_edge(new, edge)
    new[i] = j
    new[ip] = jp
    return new


def main(file):
    s, edges = open(file).read().splitlines()
    graph = parse_edge_string(s)
    edges = ints(edges)
    new = two_break_on_genome_graph(graph, *edges)
    edges = [(k, v) for k, v in new.items()]
    print(*edges, sep=", ")
