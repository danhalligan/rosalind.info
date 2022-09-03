# Implement PatternToNumber


def symbol_to_number(x):
    return ["A", "C", "G", "T"].index(x)


def pattern_to_number(seq):
    if len(seq) == 0:
        return 0
    else:
        return 4 * pattern_to_number(seq[:-1]) + symbol_to_number(seq[-1])


def main(file):
    seq = open(file).read().splitlines()[0]
    print(pattern_to_number(seq))
