# Construct the Graph of a Spectrum

from .ba4c import mass


def main(file):
    m = {v: k for k, v in mass().items()}
    x = [0] + list(map(int, open(file).read().split()))
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            d = x[j] - x[i]
            if d in m:
                print(f"{x[i]}->{x[j]}:{m[d]}")
