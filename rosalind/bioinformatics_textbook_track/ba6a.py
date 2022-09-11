# Implement GreedySorting


def reversal(perm, a, b):
    perm[a : b + 1] = [-x for x in perm[a : b + 1][::-1]]


def locate(x, perm):
    return [abs(y) for y in perm].index(x)


def format_perm(perm):
    return "(" + " ".join([f"{x:+}" for x in perm]) + ")"


def read_perm(s):
    return list(map(int, s[1:-1].split()))


def greedy_sorting(perm):
    x = 1
    while x <= len(perm):
        if perm[x - 1] == x:
            x += 1
        else:
            reversal(perm, x - 1, locate(x, perm))
            yield perm


def main(file):
    s = open(file).read().rstrip()
    for step in greedy_sorting(read_perm(s)):
        print(format_perm(step))
