# Dictionaries

from collections import Counter


def main(file):
    for k, v in Counter(open(file).read().split()).items():
        print(k, v)
