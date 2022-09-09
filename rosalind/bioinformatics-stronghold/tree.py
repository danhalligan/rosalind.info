# Completing a Tree

from .helpers import Parser


# Nb. A connected tree of n nodes will always contain n-1 edges
def main(file):
    data = Parser(file).lines()
    n_nodes = int(data[0])
    print(n_nodes - len(data[1:]) - 1)
