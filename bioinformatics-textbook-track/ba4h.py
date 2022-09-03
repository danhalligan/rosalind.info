# Generate the Convolution of a Spectrum

from collections import defaultdict, Counter


def spectrum_convolution(spec):
    dict = defaultdict(int)
    for i in range(len(spec)):
        for j in range(i + 1, len(spec)):
            m = abs(spec[i] - spec[j])
            if m > 0:
                dict[abs(spec[i] - spec[j])] += 1
    return Counter(dict).most_common()


def main(file):
    mass = open(file).read().split()
    mass = list(map(int, mass))
    res = spectrum_convolution(mass)
    print(*sum([[k] * v for k, v in res], []))
