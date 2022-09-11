# Compute the Score of a Cyclic Peptide Against a Spectrum

from .ba4c import mass, cyclo_spectrum


def score(theoretical, expected):
    spec, score = theoretical[:], 0
    if spec:
        for m in expected:
            if m in spec:
                score += 1
                spec.remove(m)
    return score


def main(file):
    peptide, spectrum = open(file).read().splitlines()
    spectrum = list(map(int, spectrum.split()))
    m = mass()
    spec = cyclo_spectrum([m[x] for x in peptide])
    print(score(spec, spectrum))
