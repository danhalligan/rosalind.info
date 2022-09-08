# Translating RNA into Protein

from .helpers import Parser


def main(file):
    print(Parser(file).rna().translate())
