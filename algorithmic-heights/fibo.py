# Fibonacci Numbers


def fibo(n):
    x, y = 0, 1
    for i in range(n):
        x, y = y, x + y
    return x


def main(file):
    n = int(open(file).read())
    print(fibo(n))
