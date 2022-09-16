# Implement TreeColoring

from collections import defaultdict
from .ba9e import color_tree


def parse_input(file):
    f = open(file)
    g = defaultdict(list)
    for line in f:
        line = line.rstrip()
        if line == "-":
            break
        x, nodes = line.split(" -> ")
        if nodes == "{}":
            continue
        for n in nodes.split(","):
            g[int(x)].append({"n": int(n)})

    colors = defaultdict(lambda: None)
    for line in f:
        line = line.rstrip()
        x, c = line.split(": ")
        colors[int(x)] = c
    return g, colors


def main(file):
    g, colors = parse_input(file)
    colors = color_tree(g, colors)
    for k in sorted(colors.keys()):
        print(f"{k}: {colors[k]}")
