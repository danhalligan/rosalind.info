# Degree Array

from .helpers import parse_graph
import sys


def main(file):
    graph = parse_graph(open(file))
    print(*[len(graph[node]) for node in graph.keys()])


if __name__ == "__main__":
    main(sys.argv[1])
