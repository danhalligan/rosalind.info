# Compute the Hamming Distance Between Two Strings


def hamming(s1, s2):
    return sum(s1[i] != s2[i] for i in range(len(s1)))


def main(file):
    s1, s2 = open(file).read().splitlines()
    print(hamming(s1, s2))
