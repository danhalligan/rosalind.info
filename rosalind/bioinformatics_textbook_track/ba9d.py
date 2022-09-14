# Find the Longest Repeat in a String

from rosalind.bioinformatics_stronghold.suff import suff


# We'll do DFS over tree, ending iteration in the node above a leaf
# By doing this, edges returned are shared and therefore repeats
def internal_edges(tree):
    for n1 in tree.keys():
        if not len(tree[n1]):
            yield ""
        for n2 in internal_edges(tree[n1]):
            yield n1 + n2


def longest_shared_substring(tree):
    return max(internal_edges(tree), key=lambda x: len(x))


def main(file):
    seq = open(file).read().rstrip()
    tree = suff(seq)
    print(longest_shared_substring(tree))
