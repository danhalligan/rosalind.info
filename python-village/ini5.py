import sys


def main(file):
    for i, line in enumerate(open(file).readlines()):
        if i % 2 == 1:
            print(line, end="")


if __name__ == "__main__":
    main(sys.argv[1])
if __name__ == "__main__":
    main(sys.argv[1])
