# 2-Satisfiability

# https://en.wikipedia.org/wiki/2-satisfiability#Strongly_connected_components

from .helpers import ints, recursionlimit
from .scc import scc
from .sc import condense
from .ts import topological_sort


def parse_twosat(handle):
    info = next(handle)
    if info == "\n":
        info = next(handle)
    nodes, n_edges = ints(info)
    edges = [next(handle) for _ in range(n_edges)]
    graph = {}

    for n in range(1, nodes + 1):
        graph[n] = list()
        graph[-n] = list()

    for edge in edges:
        f, t = ints(edge)
        graph[-f].append(t)
        graph[-t].append(f)

    return graph


def parse_twosats(handle):
    n = int(next(handle))
    for _ in range(n):
        yield parse_twosat(handle)


def find_comp(n, components):
    for j, comp2 in enumerate(components):
        if n in comp2:
            return j


def twosat(graph):
    components = list(scc(graph))
    ngraph = condense(graph, components)

    # If any component has a literal and its negation,the instance is not
    # satisfiable, so return 0
    for i in topological_sort(ngraph):
        for a in components[i]:
            if -a in components[i]:
                return 0, []

    # If here, then the instance is satisfiable, let's recover assignment
    assignment = []
    for i in topological_sort(ngraph)[::-1]:
        for v in components[i]:
            if v not in assignment and -v not in assignment:
                assignment.append(v)

    return 1, sorted(assignment, key=lambda x: abs(x))


def main(file):
    with recursionlimit(5000):
        graphs = parse_twosats(open(file))
        for graph in graphs:
            x, a = twosat(graph)
            print(x, *a)
