# Implement LeaderboardCyclopeptideSequencing

from .ba4c import mass, cyclo_spectrum
from .ba1a import substrings
from .ba4f import score
from .ba4e import expand


def linear_spectrum(peptide):
    spec = [0]
    for i in range(1, len(peptide) + 1):
        for x in substrings(peptide, i):
            spec.append(sum(x))
    return spec


def linear_score(peptide, spectrum):
    return score(linear_spectrum(peptide), spectrum)


def cyclo_score(peptide, spectrum):
    return score(cyclo_spectrum(peptide), spectrum)


def cut(peptides, spectrum, n):
    if len(peptides) < n:
        return peptides
    sc = [linear_score(p, spectrum) for p in peptides]
    lim = sorted(sc, reverse=True)[n - 1]
    return [p for p, sc in zip(peptides, sc) if sc >= lim]


def leaderboard_cyclopeptide_sequencing(spec, n, masses):
    lb = [[]]
    leader = []
    while len(lb):
        lb = expand(lb, masses)
        for pep in lb.copy():
            if sum(pep) == spec[-1]:
                if cyclo_score(pep, spec) > cyclo_score(leader, spec):
                    leader = pep
            elif sum(pep) > spec[-1]:
                lb.remove(pep)
        lb = cut(lb, spec, n)
    return leader


def main(file):
    n, m = open(file).read().splitlines()
    spec = list(map(int, m.split()))
    masses = set(mass().values())
    p = leaderboard_cyclopeptide_sequencing(spec, int(n), masses)
    print(*p, sep="-")
