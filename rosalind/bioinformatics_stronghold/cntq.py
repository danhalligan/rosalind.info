# Counting Quartets

from math import comb


def main(file):
    n = int(open(file).readline())
    print(comb(n, 4) % 1000000)
