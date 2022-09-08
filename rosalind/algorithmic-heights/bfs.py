# Breadth-First Search

from .helpers import parse_graph
import sys


def bfs(graph):
    n = len(graph)
    d = [-1 for i in range(n + 1)]
    d[1] = 0
    q = [1]
    while q:
        u = q.pop(0)
        for v in graph[u]:
            if d[v] == -1:
                q.append(v)
                d[v] = d[u] + 1
    return d[1:]


def main(file):
    graph = parse_graph(open(file), directed=True)
    print(*bfs(graph))


if __name__ == "__main__":
    main(sys.argv[1])
