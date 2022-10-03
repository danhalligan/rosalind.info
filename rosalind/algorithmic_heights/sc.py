# Semi-Connected Graph

from .helpers import parse_graphs
from .scc import scc
from .hdag import hdag

# The idea here is we take strongly connected components
# Then we build a new graph with edges between components
# If there is a hamiltonian path in this new graph, then we are semi connected


def find_comp(n, components):
    for j, comp2 in enumerate(components):
        if n in comp2:
            return j


def sc(graph):
    ngraph = {}
    components = list(scc(graph))
    for i, comp in enumerate(components):
        ngraph[i] = set(
            [
                find_comp(dest, components)
                for node in comp
                for dest in graph[node]
                if dest not in comp
            ]
        )
    return -1 if hdag(ngraph) == [-1] else 1


def main(file):
    graphs = parse_graphs(open(file), directed=True)
    print(*[sc(g) for g in graphs])
