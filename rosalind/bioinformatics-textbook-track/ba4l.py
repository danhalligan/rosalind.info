# Trim a Peptide Leaderboard

from .ba4c import mass
from .ba4g import linear_score


def cut(peptides, spectrum, n):
    masses = mass()
    specs = [[masses[x] for x in p] for p in peptides]
    sc = [linear_score(p, spectrum) for p in specs]
    lim = sorted(sc, reverse=True)[n - 1]
    return [p for p, sc in zip(peptides, sc) if sc >= lim]


def main(file):
    lb, spectrum, n = open(file).read().splitlines()
    lb = lb.split()
    spectrum = list(map(int, spectrum.split()))
    print(*cut(lb, spectrum, int(n)))
