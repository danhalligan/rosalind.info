def ints(x):
    return list(map(int, x.split()))


def parse_graph(handle, directed=False):
    info = next(handle)
    nodes, n_edges = ints(info)
    edges = [next(handle) for i in range(n_edges)]
    graph = {}

    for n in range(1, nodes + 1):
        graph[n] = set()

    for edge in edges:
        f, t = ints(edge)
        graph[f].add(t)
        if not directed:
            graph[t].add(f)

    return graph


def parse_weighed_graph(file):
    info, *edges = open(file).read().splitlines()
    nodes, n_edges = ints(info)
    graph = {}

    for n in range(1, nodes + 1):
        graph[n] = list()

    for edge in edges:
        f, t, w = ints(edge)
        graph[f].append({"n": t, "w": w})

    return graph


def parse_graphs(handle, directed=False):
    n = int(next(handle))
    for i in range(n):
        next(handle)
        yield parse_graph(handle, directed=directed)
