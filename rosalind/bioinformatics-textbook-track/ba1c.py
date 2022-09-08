# Find the Reverse Complement of a String


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def main(file):
    seq = open(file).read().splitlines()[0]
    print(revcomp(seq))
