# Creating a Character Table from Genetic Strings


def main(file):
    seqs = open(file).read().splitlines()
    for pos in zip(*seqs):
        pos = [int(x == pos[0]) for x in pos]
        if 1 < sum(pos) < (len(seqs) - 1):
            print(*pos, sep="")
