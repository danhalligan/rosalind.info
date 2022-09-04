# Connected Components

from .helpers import parse_graph


def find_component(node, graph):
    def visit(node, visited):
        visited.add(node)
        for n in list(set(graph[node]) - visited):
            visit(n, visited)
        return visited

    return visit(node, set())


def find_components(graph):
    nodes = set(graph.keys())
    components = list()
    while len(nodes):
        res = find_component(next(iter(nodes)), graph)
        nodes = nodes - res
        components.append(res)
    return components


def main(file):
    graph = parse_graph(open(file))
    print(len(find_components(graph)))
