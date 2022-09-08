# Implement ChromosomeToCycle

from .ba6a import read_perm


def chromosome2cycle(perm):
    nodes = []
    for i in perm:
        if i > 0:
            nodes += [2 * i - 1, 2 * i]
        else:
            nodes += [-2 * i, -2 * i - 1]
    return nodes


def ba6f(s):
    return chromosome2cycle(read_perm(s))


def format_cycle(cycle):
    return "(" + " ".join(map(str, cycle)) + ")"


def main(file):
    s = open(file).read().rstrip()
    print(format_cycle(ba6f(s)))
