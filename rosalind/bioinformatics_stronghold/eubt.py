# Enumerating Unrooted Binary Trees

from copy import deepcopy


class Node:
    def __init__(self, name=None):
        self.name = name

    def leaf(self):
        return self.name is not None


class Tree:
    def __init__(self, edges=[]):
        self.edges = edges

    def nodes(self):
        nodes = [n for edge in self.edges for n in edge]
        return list(dict.fromkeys(nodes))

    def internal_nodes(self):
        return [node for node in self.nodes() if not node.leaf()]

    def leaf_edges(self):
        return [edge for edge in self.edges if edge[0].leaf() or edge[1].leaf()]

    def adjacent_edges(self, node):
        return [edge for edge in self.edges if node in edge]

    def adjacent_nodes(self, node):
        return [n for edge in self.adjacent_edges(node) for n in edge if n is not node]

    def newick(self):
        tree = deepcopy(self)
        n = tree.nodes()
        if len(n) == 1:
            return f"{n[0].name};"
        elif len(n) == 2:
            return f"({n[0].name},{n[1].name});"
        elif len(tree.nodes()) > 2:
            for node in tree.internal_nodes():
                adjacent_nodes = tree.adjacent_nodes(node)
                adjacent_leaves = [n for n in adjacent_nodes if n.leaf()]

                if len(adjacent_leaves) >= 2:
                    leaf1, leaf2 = adjacent_leaves[0:2]
                    for edge in tree.adjacent_edges(node):
                        if leaf1 in edge or leaf2 in edge:
                            tree.edges.remove(edge)
                    node.name = f"({leaf1.name},{leaf2.name})"
                    return tree.newick()

    def insert(self, i, leaf):
        tree = deepcopy(self)
        n1, n2 = tree.edges.pop(i)
        ni = Node()
        tree.edges += [(n1, ni), (n2, ni), (leaf, ni)]
        return tree


def enumerate_trees(leaves):
    if len(leaves) == 2:
        n1, n2 = leaves
        yield Tree([(Node(n1), Node(n2))])
    elif len(leaves) > 2:
        leaf = Node(leaves[-1])
        for t in enumerate_trees(leaves[:-1]):
            for i in range(len(t.edges)):
                yield t.insert(i, leaf)


def main(file):
    sp = open(file).read().split()
    for tree in enumerate_trees(sp):
        print(tree.newick())
