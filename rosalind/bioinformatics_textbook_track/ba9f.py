# Find the Shortest Non-Shared Substring of Two Strings

# We have to be careful here as the title is misleading. We actually need
# The shortest substring of **Text1** that does not appear in **Text2**
#
# The logic here is that we can find the shortest a non-shared sequence by
# finding the sequence from the root up to a non-purple node.
# For the shortest non-shared sequence, we don't need to consider nodes
# deeper in the tree and we also only need the first character of the final
# edge to the non-purple node (since as soon as we add this character, we
# have a subsequence unique to one of the sequences).
#
# Since we're only interested in subsequences present in Text1, we condition
# on our final node being red.


from rosalind.bioinformatics_stronghold.suff import suff
from .ba9e import as_graph, color_tree, leaf_colors


def non_purple_edges(tree, colors):
    def _non_purple_edges(node, seq, path):
        if colors[node] == "purple":
            for child in tree[node]:
                yield from _non_purple_edges(
                    child["n"], seq + [child["l"]], path + [node]
                )
        else:
            if colors[node] == "red":
                seq[-1] = seq[-1][0]
                seq = "".join(seq)
                if "#" not in seq and "$" not in seq:
                    yield seq

    yield from _non_purple_edges(0, [], [])


def shortest_nonshared_substring(seq1, seq2):
    tree = suff(seq1 + "$" + seq2 + "#")
    tree = as_graph(tree)
    colors = color_tree(tree, leaf_colors(tree))
    return min(non_purple_edges(tree, colors), key=lambda x: len(x))


def main(file):
    seq1, seq2 = open(file).read().splitlines()
    print(shortest_nonshared_substring(seq1, seq2))
