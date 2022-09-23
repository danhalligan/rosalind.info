# Compute Distances Between Leaves

from re import split
from collections import defaultdict
from math import inf
from heapq import heappush, heappop

# For this solution, since the first n nodes are the leaves
# I run Djiakstra's algorithm starting at the first n nodes,
# and keep distances for the first n termini for each.
# This gives the required matrix, but it's not terribly efficient.

# A smarter solution might be to drop internal nodes and then
# read off the distances between leaves directly


def nodes(graph):
    s = list(graph.keys())
    e = [y["n"] for v in graph.values() for y in v]
    return set(s) | set(e)


# Dijkstra's algorithm to find distance from start to all other nodes
# Assumes nodes are integers starting a 0!
def dij(start, graph):
    d = [inf for i in range(len(nodes(graph)))]
    d[start] = 0
    q = []
    heappush(q, (0, start))
    processed = set()

    while q:
        u = heappop(q)[1]
        processed.add(u)
        for v in graph[u]:
            if v["n"] not in processed:
                d[v["n"]] = min(d[u] + v["w"], d[v["n"]])
                heappush(q, (d[v["n"]], v["n"]))

    return d


def parse_weighted_graph(edges):
    graph = defaultdict(list)

    for edge in edges:
        f, t, w = map(int, split(r"\D+", edge))
        graph[f].append({"n": t, "w": w})

    return graph


def main(file):
    n_leaves, *edges = open(file).read().splitlines()
    n_leaves = int(n_leaves)
    graph = parse_weighted_graph(edges)
    for i in range(n_leaves):
        print(*dij(i, graph)[:n_leaves])
