# Hamiltonian Path in DAG

from .ts import topological_sort
from .helpers import parse_graphs

# A Hamiltonian path exists if and only if there are edge between consecutive
# topologically sorted vertices


def hdag(graph):
    x = topological_sort(graph)
    for a, b in zip(x, x[1:]):
        if b not in graph[a]:
            return [-1]
    return [1] + x


def main(file):
    for graph in parse_graphs(open(file), directed=True):
        print(*hdag(graph))
