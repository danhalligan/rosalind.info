# Dijkstra's Algorithm

from helpers import parse_weighed_graph

# from collections import Counter
from math import inf
from heapq import heappush, heappop


def dij(graph):
    d = [inf for i in range(len(graph) + 1)]
    d[1] = 0
    q = []
    heappush(q, (0, 1))
    processed = set()

    while q:
        u = heappop(q)[1]
        processed.add(u)
        for v in graph[u]:
            if v["n"] not in processed:
                d[v["n"]] = min(d[u] + v["w"], d[v["n"]])
                heappush(q, (d[v["n"]], v["n"]))

    return [-1 if x == inf else x for x in d[1:]]


def main(file):
    graph = parse_weighed_graph(file)
    print(*dij(graph))
