# Implement NumberToPattern


def number_to_symbol(x):
    return ["A", "C", "G", "T"][x]


def number_to_pattern(i, k):
    if k == 1:
        return number_to_symbol(i)
    else:
        pi, r = divmod(i, 4)
        return number_to_pattern(pi, k - 1) + number_to_symbol(r)


def main(file):
    index, k = open(file).read().splitlines()
    print(number_to_pattern(int(index), int(k)))
