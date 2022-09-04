# Double-Degree Array

from .helpers import parse_graph


def main(file):
    graph = parse_graph(open(file))
    print(*[sum(len(graph[k]) for k in graph[n]) for n in graph.keys()])
