import sys

from collections import Counter


def main(file):
    for k, v in Counter(open(file).read().split()).items():
        print(k, v)


if __name__ == "__main__":
    main(sys.argv[1])
