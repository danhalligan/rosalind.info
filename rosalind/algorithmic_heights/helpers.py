def ints(x):
    return list(map(int, x.split()))


def parse_graph(handle, directed=False, weighted=False):
    info = next(handle)
    if info == "\n":
        info = next(handle)
    nodes, n_edges = ints(info)
    edges = [next(handle) for _ in range(n_edges)]
    graph = {}

    for n in range(1, nodes + 1):
        graph[n] = list()

    for edge in edges:
        if weighted:
            f, t, w = ints(edge)
            graph[f].append({"n": t, "w": w})
            if not directed:
                graph[t].append({"n": f, "w": w})
        else:
            f, t = ints(edge)
            graph[f].append(t)
            if not directed:
                graph[t].append(f)

    return graph


def parse_graphs(handle, directed=False, weighted=False):
    n = int(next(handle))
    for _ in range(n):
        yield parse_graph(handle, directed=directed, weighted=weighted)
