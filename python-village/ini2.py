# Variables and Some Arithmetic


def main(file):
    print(sum([int(x) ** 2 for x in open(file).read().split()]))
