# Degree Array

from .helpers import parse_graph


def main(file):
    graph = parse_graph(open(file))
    print(*[len(graph[node]) for node in graph.keys()])
