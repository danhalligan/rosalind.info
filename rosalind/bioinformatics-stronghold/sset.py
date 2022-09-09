# Counting Subsets

from .helpers import Parser


def main(file):
    n = Parser(file).ints()[0]
    print(2**n % 1000000)
