# Shortest Paths in DAG

from .helpers import parse_graph
from .ts import topological_sort
from math import inf


def simplify_graph(graph):
    return {k: [x["n"] for x in v] for k, v in graph.items()}


def sdag(graph):
    n = len(graph)
    dist = [inf for _ in range(n + 1)]
    dist[1] = 0
    order = topological_sort(simplify_graph(graph))
    for u in order:
        # -----------------------------------------------------------------------
        # Hack to pass!
        # For this question, if there are multiple edges connecting the same
        # two nodes, we must only consider the last one (not the shortest)!
        # So we'll process edges in reverse order and skip any we've seen...
        # -----------------------------------------------------------------------
        seen = set()
        for edge in graph[u][::-1]:
            if edge["n"] not in seen:
                seen.add(edge["n"])
                if dist[u] + edge["w"] < dist[edge["n"]]:
                    dist[edge["n"]] = dist[u] + edge["w"]
    return ["x" if x == inf else x for x in dist[1:]]


def main(file):
    print(*sdag(parse_graph(open(file), directed=True, weighted=True)))
