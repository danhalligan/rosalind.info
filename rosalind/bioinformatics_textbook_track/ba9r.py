# Construct a Suffix Tree from a Suffix Array


class Node:
    def __init__(self, parent=None, label=""):
        self.parent = parent
        self.label = label


def desc(tree, node):
    """
    length of the concatenation of all path labels from the root to node
    """
    n = len(tree[node].label)
    while tree[node].parent:
        node = tree[node].parent
        n += len(tree[node].label)
    return n


def ba9r(txt, sa, lcp):
    tree = {-1: Node()}

    for i in range(len(txt)):
        v = i - 1
        while tree[v].parent and desc(tree, v) > lcp[i]:
            v = tree[v].parent
        dv = desc(tree, v)

        if lcp[i] == dv:
            tree[i] = Node(v, txt[sa[i] + lcp[i] :])
        else:
            w = i - 1
            while tree[w].parent and tree[w].parent != v:
                w = tree[w].parent

            tree[f"y{i}"] = Node(v, txt[sa[i - 1] + dv : sa[i - 1] + lcp[i]])
            tree[w] = Node(f"y{i}", txt[sa[i - 1] + lcp[i] : sa[i - 1] + desc(tree, w)])
            tree[i] = Node(f"y{i}", txt[sa[i] + lcp[i] :])

    del tree[-1]
    return tree


def main(file):
    txt, sa, lcp = open(file).read().splitlines()
    sa = [int(x) for x in sa.split(", ")]
    lcp = [int(x) for x in lcp.split(", ")]

    tree = ba9r(txt, sa, lcp)
    labels = [tree[key].label for key in tree.keys()]
    print(*sorted(labels), sep="\n")
