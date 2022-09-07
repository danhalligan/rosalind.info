# Compute Limb Lengths in a Tree


def parse_mat(lines):
    """Parse integer matrix from set of lines"""
    return [[int(x) for x in y.split()] for y in lines]


def limb(j, d):
    """Calculate limb length j for distance matrix d"""
    limb_len = 100000
    n = len(d)
    for k in range(int(n)):
        for i in range(int(n)):
            if j != k and i != j:
                limb_len = min((d[i][j] + d[j][k] - d[i][k]) // 2, limb_len)
    return limb_len


def main(file):
    _, j, *arr = open(file).read().splitlines()
    print(limb(int(j), parse_mat(arr)))
