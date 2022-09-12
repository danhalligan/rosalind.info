# Testing Bipartiteness

from .helpers import parse_graphs


# d stores the coloring of each node as 0 or 1
def bip(g):
    d = {1: 0}
    q = [1]
    while q:
        u = q.pop(0)
        for v in g[u]:
            col = (d[u] + 1) % 2
            if v not in d:
                q.append(v)
                d[v] = col
            elif d[v] != col:
                return -1
    return 1


def main(file):
    print(*[bip(g) for g in parse_graphs(open(file))])
