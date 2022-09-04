# Implement GibbsSampler

import math
import random
from .ba1a import substrings
from .ba2d import create_profile, score
from .ba2f import minimise, rand_substring


def profile_probs(text, k, profile):
    base = {"A": 0, "C": 1, "G": 2, "T": 3}
    return [
        math.prod(profile[base[x[j]]][j] for j in range(k)) for x in substrings(text, k)
    ]


def pick_kmer(seq, k, p):
    pr = profile_probs(seq, k, p)
    c = random.choices(range(len(pr)), pr)[0]
    return seq[c : (c + k)]


def gibbs_sampler(dna, k, n):
    m = [rand_substring(x, k) for x in dna]
    best = m
    for j in range(n):
        i = random.randint(0, len(dna) - 1)
        m[i] = pick_kmer(dna[i], k, create_profile(m[:i] + m[i + 1 :], pc=1))
        if score(m) < score(best):
            best = m
    return score(best), best


def main(file):
    ints, *dna = open(file).read().splitlines()
    k, t, n = map(int, ints.split())
    print(*minimise(gibbs_sampler, 20, dna, k, n), sep="\n")
