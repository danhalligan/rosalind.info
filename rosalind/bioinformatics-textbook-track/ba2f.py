# Implement RandomizedMotifSearch

from .ba2c import pmpkmer
from .ba2d import create_profile, score
from random import randint
import math


def rand_substring(seq, k):
    i = randint(0, len(seq) - k)
    return seq[i : (i + k)]


def motifs(profile, dna):
    k = len(profile[0])
    return [pmpkmer(seq, k, profile) for seq in dna]


def randomised_motif_search(dna, k):
    m = [rand_substring(x, k) for x in dna]
    s = math.inf
    while score(m) < s:
        s = score(m)
        m = motifs(create_profile(m, pc=1), dna)
    return s, m


def minimise(fun, n, *args):
    """
    Minimise function fun by calling n times with *args.
    """
    bs, bx = fun(*args)
    for i in range(n - 1):
        s, x = fun(*args)
        if s < bs:
            bs, bx = s, x
    return bx


def main(file):
    ints, *dna = open(file).read().splitlines()
    k, t = map(int, ints.split())
    print(*minimise(randomised_motif_search, 1000, dna, k), sep="\n")
