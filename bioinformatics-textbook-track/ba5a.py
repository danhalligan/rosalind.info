# Find the Minimum Number of Coins Needed to Make Change


def dp_change(money, coins):
    mc = [0] * (money + 1)
    for m in range(1, money + 1):
        mc[m] = min(mc[m - coin] + 1 for coin in coins if m >= coin)
    return mc[money]


def main(file):
    money, coins = open(file).read().splitlines()
    coins = list(map(int, coins.split(",")))
    print(dp_change(int(money), coins))
