# Implement TrieMatching

from .ba9a import trie, edge_match


def prefix_trie_matching(text, graph):
    v = 0
    for symbol in text:
        m = edge_match(graph, v, symbol)
        if m:
            if not graph[m]:
                return True
            else:
                v = m
        else:
            return False


def trie_matching(text, graph):
    for i in range(len(text)):
        if prefix_trie_matching(text[i:], graph):
            yield i


def main(file):
    text, *seqs = open(file).read().splitlines()
    graph = trie(seqs)
    print(*trie_matching(text, graph))
