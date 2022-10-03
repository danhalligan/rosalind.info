# Strongly Connected Components

from collections import defaultdict
from .ts import topological_sort
from .helpers import parse_graph
from .cc import find_component


def reverse_graph(graph):
    rev = defaultdict(list)
    for node in graph:
        for child in graph[node]:
            rev[child].append(node)
    return rev


def scc(graph):
    order = topological_sort(graph)
    rev = reverse_graph(graph)
    while order:
        n = order.pop(0)
        res = find_component(n, rev)
        order = [x for x in order if x not in res]
        for k in rev.keys():
            rev[k] = [n for n in rev[k] if n not in res]
        yield res


def main(file):
    graph = parse_graph(open(file), directed=True)
    print(len(list(scc(graph))))
