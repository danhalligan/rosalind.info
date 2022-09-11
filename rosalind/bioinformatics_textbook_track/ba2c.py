# Find a Profile-most Probable k-mer in a String

from .ba1a import substrings
import math


def pmpkmer(text, k, profile):
    base = {"A": 0, "C": 1, "G": 2, "T": 3}
    prob = -1
    for x in substrings(text, k):
        p = math.prod(profile[base[x[j]]][j] for j in range(k))
        if p > prob:
            prob, best = p, x
    return best


def main(file):
    text, k, *profile = open(file).read().splitlines()
    profile = [list(map(float, x.split())) for x in profile]
    print(pmpkmer(text, int(k), profile))
