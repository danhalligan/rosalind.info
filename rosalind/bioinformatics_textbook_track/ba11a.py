# Construct the Graph of a Spectrum

from collections import defaultdict
from .ba4c import mass


def spectrum_graph(x):
    m = {v: k for k, v in mass().items()}
    g = defaultdict(list)
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            d = x[j] - x[i]
            if d in m:
                g[x[i]].append({"n": x[j], "l": m[d]})
    return g


def main(file):
    x = [0] + list(map(int, open(file).read().split()))
    for k, v in spectrum_graph(x).items():
        for v2 in v:
            print(f"{k}->{v2['n']}:{v2['l']}")
