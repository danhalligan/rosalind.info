# Sex-Linked Inheritance

from .helpers import Parser


def main(file):
    arr = Parser(file).floats()
    print(*[round(2 * v * (1 - v), 3) for v in arr])
