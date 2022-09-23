# Compute the Size of a Spectral Dictionary


from functools import cache
from .ba11c import masses


def dict_size(sv, T, max_score):
    @cache
    def size(i, t):
        if i == 0 and t == 0:
            return 1
        if t < 0 or i <= 0:
            return 0
        return sum(size(i - x, t - sv[i]) for x in mass)

    mass = masses().values()
    n = len(sv)
    sv = [0] + sv
    return sum(size(n, x) for x in range(T, max_score + 1))


def main(file):
    sv, T, max_score = open(file).read().splitlines()
    sv = list(map(int, sv.split()))
    print(dict_size(sv, int(T), int(max_score)))
