# Find a Position in a Genome Minimizing the Skew


def find_minima(seq):
    skew = [0]
    delta = {"G": 1, "C": -1, "A": 0, "T": 0}
    for i in range(len(seq)):
        skew.append(skew[i] + delta[seq[i]])
    m = min(skew)
    return (i for i, x in enumerate(skew) if x == m)


def main(file):
    seq = open(file).read().splitlines()[0]
    print(*find_minima(seq))
