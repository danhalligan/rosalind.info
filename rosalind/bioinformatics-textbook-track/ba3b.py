# Reconstruct a String from its Genome Path


def genome_path(dna):
    return dna[0] + "".join(x[-1] for x in dna[1:])


def main(file):
    dna = open(file).read().splitlines()
    print(genome_path(dna))
