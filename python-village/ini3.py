# Strings and Lists

import sys


def main(file):
    str, slices = open(file).read().splitlines()
    a, b, c, d = map(int, slices.split())
    print(str[a : b + 1], str[c : d + 1])


if __name__ == "__main__":
    main(sys.argv[1])
