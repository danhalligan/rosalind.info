# Implement AdditivePhylogeny

# Hands down the worst code I've written in a long time.
# TODO: refactor / rewrite!

from collections import defaultdict
from .ba7b import limb, parse_mat
from .ba7a import nodes


def find_path(graph, path, target):
    """Search of tree, returning route to target and cumulative dist"""
    if target in [x[0] for x in path]:
        yield path
    x = path[-1]
    if x[0] in graph:
        for node in graph[x[0]]:
            if node not in path:
                yield from find_path(
                    graph, path + [(node["n"], x[1] + node["w"])], target
                )


def find_leaves(D, n):
    for i in range(len(D)):
        for k in range(i + 1, len(D)):
            if D[i][k] == D[i][n] + D[n][k]:
                return i, n, k


# There may be multiple nodes between i and k, so we need to find the
# path from i to k (there's only 1) and break the appropriate edge.
def add_node(T, i, k, x, n):
    """Add node in graph D, between i and k, dist x from i, labelled n"""
    path = list(find_path(T, [(i, 0)], k))[0]
    for p, node in enumerate(path):
        if node[1] > x:
            break
    p = p - 1
    i, d1 = path[p]
    k, d2 = path[p + 1]

    # delete edge
    T[i] = list(filter(lambda x: x["n"] != k, T[i]))
    T[i].append({"n": n, "w": x - d1})
    T[n].append({"n": k, "w": d2 - x})
    return T


def additive_phylogeny(D, m):
    n = len(D) - 1
    if len(D) == 2:
        g = defaultdict(list)
        g[0].append({"n": 1, "w": D[0][1]})
        return g

    limb_len = limb(D, n)
    for j in range(len(D)):
        if j != n:
            D[j][n] = D[j][n] - limb_len
            D[n][j] = D[j][n]

    # three leaves such that Di,k = Di,n + Dn,k
    (i, n, k) = find_leaves(D, n)
    x = D[i][n]

    # remove row n and column n from D
    D = [x[:n] + x[n + 1 :] for x in D[:n] + D[n + 1 :]]

    T = additive_phylogeny(D, m)

    # label for new internal node
    v = max(max(nodes(T)), m - 1) + 1

    # break an internal edge adding a new node (possibly)
    # and add the new leaf node
    T = add_node(T, i, k, x, v)
    T[v].append({"n": n, "w": limb_len})
    return T


def as_edges(graph):
    edges = []
    for k in sorted(graph):
        for v in graph[k]:
            edges += [f"{k}->{v['n']}:{v['w']}"]
            edges += [f"{v['n']}->{k}:{v['w']}"]
    return sorted(edges)


def main(file):
    n, *D = open(file).read().splitlines()
    graph = additive_phylogeny(parse_mat(D), int(n))
    for edge in as_edges(graph):
        print(edge)
