# Square in a Graph

from .helpers import parse_graphs


def dfs(graph, start, depth):
    def recurse(node, dp, visited):
        if dp == depth:
            yield node
        if dp < depth:
            for new in graph[node]:
                if new not in visited:
                    yield from recurse(new, dp + 1, visited | set([node]))

    return recurse(start, 0, set())


def sq(graph, n):
    for node in graph:
        for tail in dfs(graph, node, n - 1):
            if tail in graph[node]:
                return 1
    return -1


def main(file):
    graphs = parse_graphs(open(file), directed=False)
    print(*[sq(g, 4) for g in graphs])
