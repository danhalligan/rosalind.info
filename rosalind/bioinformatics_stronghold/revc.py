# Complementing a Strand of DNA

from .helpers import Parser


def main(file):
    print(Parser(file).dna().revc())
