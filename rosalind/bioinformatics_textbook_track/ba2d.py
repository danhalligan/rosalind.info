# Implement GreedyMotifSearch

from .ba1a import substrings
from .ba2c import pmpkmer
from collections import Counter


def create_profile(seqs, pc=0):
    bases = ["A", "C", "G", "T"]
    profile = [[] for b in bases]
    for i, b in enumerate(bases):
        profile[i] = [
            (sum(x[j] == b for x in seqs) + pc) / len(seqs) for j in range(len(seqs[0]))
        ]
    return profile


def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        bases = [x[i] for x in motifs]
        c = Counter(bases).most_common()[0][0]
        score += sum(x != c for x in bases)
    return score


def greedy_motif_search(dna, k, pc=0):
    best_motifs = [x[0:k] for x in dna]
    for kmer in substrings(dna[0], k):
        motifs = [kmer]
        for i in range(1, len(dna)):
            motifs += [pmpkmer(dna[i], k, create_profile(motifs, pc=pc))]
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def main(file):
    ints, *dna = open(file).read().splitlines()
    k, t = map(int, ints.split())
    print(*greedy_motif_search(dna, k), sep="\n")
