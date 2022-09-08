# Transcribing DNA into RNA

from .helpers import Parser


def main(file):
    print(Parser(file).dna().rna())
