import sys


def main(file):
    print(sum([int(x) ** 2 for x in open(file).read().split()]))


if __name__ == "__main__":
    main(sys.argv[1])
