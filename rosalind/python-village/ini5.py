# Working with Files


def main(file):
    for i, line in enumerate(open(file).readlines()):
        if i % 2 == 1:
            print(line, end="")
