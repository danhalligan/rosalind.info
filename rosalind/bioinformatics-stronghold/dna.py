# Counting DNA Nucleotides

from .helpers import Parser


def main(file):
    print(*Parser(file).dna().table().values())
