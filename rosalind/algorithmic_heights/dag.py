# Testing Acyclicity

from .helpers import parse_graphs

# Recursively remove leaf nodes (with no onward connections)
# If we remove all nodes, graph must be acyclic.


def graph_nodes(graph):
    nodes = set()
    for node, val in graph.items():
        nodes.add(node)
        nodes = nodes.union(val)
    return nodes


def remove_leaves(graph):
    nodes = graph_nodes(graph)
    leaves = nodes - graph.keys()
    return {n: set(v) - leaves for n, v in graph.items() if len(set(v) - leaves)}


def dag(graph):
    while graph:
        ngraph = remove_leaves(graph)
        if len(graph) == len(ngraph):
            return -1
        graph = ngraph
    return 1


def main(file):
    gen = parse_graphs(open(file), directed=True)
    print(*[dag(graph) for graph in gen])
