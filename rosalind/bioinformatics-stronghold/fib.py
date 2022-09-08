# Rabbits and Recurrence Relations

from .helpers import Parser


def fib(n, k):
    a, b = 1, 1
    for _ in range(2, n):
        a, b = b, k * a + b
    return b


def main(file):
    print(fib(*Parser(file).ints()))
