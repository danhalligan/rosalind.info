# Reconstruct a String from its Paired Composition

from .ba3b import genome_path
from .ba3d import dbru
from .ba3g import count_connections


def one_in_out(c):
    return c["in"] == 1 and c["out"] == 1


def maximal_nonbranching_paths(graph):
    con = count_connections(graph)
    paths = []
    for n, c in con.items():
        if not one_in_out(c) and c["out"] > 0:
            for e in graph[n]:
                path = [n, e]
                while one_in_out(con[e]):
                    path += graph[e]
                    e = graph[e][0]
                paths += [path]

    seen = sum(paths, [])
    todo = [x for x in con.keys() if one_in_out(con[x]) and x not in seen]
    while len(todo):
        e = todo.pop(0)
        path = [e]
        while one_in_out(con[e]) and graph[e][0] in todo:
            path += graph[e]
            e = graph[e][0]
            todo.remove(e)
        paths += [path + [path[0]]]
    return paths


def main(file):
    dna = open(file).read().splitlines()
    graph = dbru(dna)
    for path in sorted(maximal_nonbranching_paths(graph)):
        print(genome_path(path))
