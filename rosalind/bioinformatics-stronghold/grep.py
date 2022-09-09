# Genome Assembly with Perfect Coverage and Repeats

from functools import cache
from .helpers import Parser
from .dbru import dbru


def drop_edge(edges, edge):
    g = list(edges)
    g.remove(edge)
    return tuple(g)


@cache
def find_paths(edges, key, assembly):
    if len(edges) == 0:
        yield assembly[: -(len(key) + 1)]
    else:
        opts = [b for a, b in edges if a == key]
        for nkey in opts:
            new = drop_edge(edges, (key, nkey))
            yield from find_paths(new, nkey, assembly + nkey[-1])


def grep(seqs):
    """Genome Assembly with Perfect Coverage and Repeats"""
    db = tuple(dbru(seqs, rev=False))
    res = set(find_paths(db, seqs[0][1:], seqs[0]))
    return sorted(res)


def main(file):
    seqs = Parser(file).lines()
    res = sorted(grep(seqs))
    print("\n".join(res))
