# Generate the Theoretical Spectrum of a Cyclic Peptide

import yaml


def mass():
    with open("resources/mass.yaml", "r") as stream:
        return yaml.safe_load(stream)


def cyclo_spectrum(peptide):
    spec = [0, sum(peptide)]
    for i in range(1, len(peptide)):
        for j in range(len(peptide)):
            spec += [sum((peptide[j:] + peptide[:j])[0:i])]
    return sorted(spec)


def main(file):
    peptide = open(file).read().splitlines()[0]
    m = mass()
    print(*cyclo_spectrum([m[x] for x in peptide]))
