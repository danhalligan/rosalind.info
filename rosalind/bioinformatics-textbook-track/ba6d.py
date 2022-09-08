# Find a Shortest Transformation of One Genome into Another by 2-Breaks

from .ba6c import find_component, parse_genome_graph, breakpoint_graph
from .ba6a import format_perm


def find_components(graph):
    nodes = set(graph.keys())
    while nodes:
        res = find_component(next(iter(nodes)), graph)
        nodes = nodes - res
        yield res


def non_trivial_cycle_nodes(graph):
    for c in find_components(graph):
        if len(c) > 2:
            return list(c)
    return None


def find_genome_component(node, graph):
    q = [node]
    visited = list()
    while q:
        node = q.pop(0)
        visited.append(node)
        for n in graph[node]:
            if -n not in visited:
                q += [-n]
    return visited


def format_genome_graph(g):
    nodes = set(g.keys())
    components = []
    while nodes:
        comp = find_genome_component(next(iter(nodes)), g)
        nodes = nodes - set(comp)
        nodes = nodes - set(-x for x in comp)
        components += [comp]
    x = [format_perm(c) for c in components]
    return "".join(x)


def add_edge(g, i, j):
    g[i] += [j]
    g[j] += [i]


def del_edge(g, i, j):
    g[i].remove(j)
    g[j].remove(i)


def ba6d(P, Q):
    bg = breakpoint_graph(P, Q)
    nodes = non_trivial_cycle_nodes(bg)
    yield format_genome_graph(P)
    while nodes:
        # arbitrary (first) edge from blue (q) edges in non-trivial cycle
        j = nodes[0]
        i2 = Q[nodes[0]][0]

        # edge from ref edges linking to j
        i = P[j][0]

        # edge from ref edges linking to i2
        j2 = P[i2][0]

        # red (p) edges with (i,j) and (i2, j2) removed
        del_edge(P, i, j)
        del_edge(P, i2, j2)

        # red (p) edges with (j, i2) and (j2, i) added
        add_edge(P, j, i2)
        add_edge(P, j2, i)

        yield format_genome_graph(P)
        bg = breakpoint_graph(P, Q)
        nodes = non_trivial_cycle_nodes(bg)


def main(file):
    P, Q = [parse_genome_graph(s) for s in open(file).read().splitlines()]
    for g in ba6d(P, Q):
        print(g)
