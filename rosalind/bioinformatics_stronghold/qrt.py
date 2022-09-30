# Quartets

from itertools import combinations, product


def qrt(s):
    taxa = [[i for i, x in enumerate(s) if x == v] for v in ["0", "1"]]
    for x in product(combinations(taxa[0], 2), combinations(taxa[1], 2)):
        yield tuple(sorted(x))


def format_pair(pair, names):
    return "{" + ", ".join([names[x] for x in pair]) + "}"


def main(file):
    names, *splits = open(file).read().splitlines()
    names = names.split()
    for pairs in set().union(*[set(qrt(s)) for s in splits]):
        print(*[format_pair(pair, names) for pair in pairs])
