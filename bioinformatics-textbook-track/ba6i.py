# Implement GraphToGenome

import re
from .ba6g import cycle2chromosome
from .ba6a import format_perm
from copy import copy


def ints(x):
    return list(map(int, x.split(", ")))


def first_key(g):
    return list(g.keys())[0]


# Find a single cycle from colored edges
def find_node_cycle(graph):
    start = first_key(graph)
    a = start
    component = []
    while graph:
        b = graph.pop(a)
        graph.pop(b)
        n = b + 1 if b % 2 else b - 1
        if n == start:
            return [b] + component + [a]
        component += [a, b]
        a = n


# find cycles in colored edges
# to do this, we first make each edge "undirected"
def find_node_cycles(graph):
    graph = copy(graph)
    for k, v in list(graph.items()):
        graph[v] = k
    while graph:
        yield find_node_cycle(graph)


def graph2genome(genome_graph):
    graph = []
    for nodes in find_node_cycles(genome_graph):
        graph += [cycle2chromosome(nodes)]
    return graph


def parse_edge_string(s):
    g = {}
    for x in re.findall(r"\((.+?)\)", s):
        a, b = ints(x)
        g[a] = b
        # g[b] = a
    return g


def main(file):
    s = open(file).read().rstrip()
    genome_graph = parse_edge_string(s)
    print(*[format_perm(x) for x in graph2genome(genome_graph)], sep="")
