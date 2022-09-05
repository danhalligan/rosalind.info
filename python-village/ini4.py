# Conditions and Loops

import sys


def main(file):
    a, b = map(int, open(file).read().split())
    print(sum([x for x in range(a, b) if x % 2 == 1]))


if __name__ == "__main__":
    main(sys.argv[1])
