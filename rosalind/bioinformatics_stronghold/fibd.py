# Mortal Fibonacci Rabbits

from .helpers import Parser


def fibd(n, m):
    v = [1] + (m - 1) * [0]
    for i in range(2, n + 1):
        v = [sum(v[1:])] + v[:-1]
    return sum(v)


def main(file):
    print(fibd(*Parser(file).ints()))
