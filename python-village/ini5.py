def main(file):
    lines = open(file).read().splitlines()
    print(*[lines[i] for i in range(1, len(lines), 2)], sep="\n")
