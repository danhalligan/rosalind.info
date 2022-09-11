# Find an Eulerian Path in a Graph

from collections import defaultdict
from copy import deepcopy
from .ba3f import parse_edges, eulerian_cycle


def count_connections(graph):
    counts = defaultdict(lambda: {"in": 0, "out": 0})
    v = sum(graph.values(), [])
    for n in sorted(set(v)):
        counts[n]["in"] = v.count(n)
    for n in graph:
        counts[n]["out"] = len(graph[n])
    return counts


def eulerian_path(graph):
    graph = deepcopy(graph)
    con = count_connections(graph)
    start = [k for k in con if con[k]["in"] > con[k]["out"]][0]
    end = [k for k in con if con[k]["in"] < con[k]["out"]][0]
    graph[start] = [end]
    cycle = eulerian_cycle(graph)[:-1]
    i = cycle.index(end)
    return cycle[i:] + cycle[:i]


def main(file):
    edges = open(file).read().splitlines()
    path = eulerian_path(parse_edges(edges))
    print("->".join(path))
