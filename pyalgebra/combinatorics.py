import functools
import itertools
from fractions import Fraction as Frac


def factorial(n):
    assert type(n) == int and n >= 0
    ans = 1
    for i in range(n):
        ans *= (i + 1)
    return ans

@functools.cache
def binomial(x, n): #x choose n
    assert type(n) == type(x) == int
    if n > x or 0 > n:
        return 0
    if n == 0 or n == x:
        return 1
    else:
        return binomial(x - 1, n) + binomial(x - 1, n - 1)


#ordered length n subsets of X with repeats: lists
#counted by: x^n
def sequences(n, X):
    assert type(n) == int
    X = set(X)
    if n == 0:
        yield []
    else:
        for first in X:
            for rest in sequences(n - 1, X):
                yield [first] + rest

#ordered length n subsets of X with no repeats: lists
#counted by: x! / (x-n)!
def permutations(n, X):
    assert type(n) == int
    X = set(X)
    if n == 0:
        yield []
    else:
        X = list(X)
        for idx in range(len(X)):
            first = X[idx]
            Y = X[:idx] + X[idx+1:]
            for rest in permutations(n - 1, Y):
                yield [first] + rest

#unordered length n subsets of X with no repeats: sets
#counted by: x choose n
def combinations(n, X):
    assert type(n) == int
    X = set(X)
    if n == 0:
        yield set()
    else:
        X = list(X)
        for idx in range(len(X)):
            first = X[idx]
            Y = X[idx+1:]
            for rest in combinations(n - 1, Y):
                yield {first} | rest

#unordered length n subsets of X with repeats: counts
#counted by: x + n - 1 choose n
def multisubsets(n, X):
    assert type(n) == int
    X = list(set(X))
    symbols = list(range(n + len(X) - 1))
    for bars in combinations(len(X) - 1, symbols):
        bars = sorted([-1] + list(bars) + [n + len(X) - 1])
        counts = [bars[i + 1] - bars[i] - 1 for i in range(len(bars) - 1)]        
        assert sum(counts) == n
        assert len(counts) == len(X)
        yield {x : counts[idx] for idx, x in enumerate(X)}

#partitions of set N into x subsets
#counted by: Striling number of the second kind {n, x}
def set_partition_eq(N, x):
    N = list(N)
    assert type(x) == int and x >= 0
    if x == 0 and len(N) == 0:
        yield []
    elif x == 0 or len(N) == 0:
        pass
    else:
        n0 = N[0]
        N_rest = N[1:]
        #all partitions come from:
        #1) {n0} & partitions of N_rest into x - 1 parts
        #2) partitions of N_rest into x parts, where n0 is added to a member of the partition

        for part in set_partition_eq(N_rest, x - 1):
            part.append([n0])
            yield part

        for rest_part in set_partition_eq(N_rest, x):
            for idx in range(len(rest_part)):
                part = [subset[:] for subset in rest_part]
                part[idx].append(n0)
                yield part

#partitions of set N into <= x subsets
def set_partition_le(N, x):
    for y in range(x + 1):
        yield from set_partition_eq(N, y)

#counted by: Bell number Bn
def set_partitions(N):
    return set_partition_le(N, len(N))

#ordered partitions of set N into x subsets
#counted by: x!{n stirling2 x}
def set_compositions(N, x):
    for part in set_partition_eq(N, x):
        for perm in permutations(x, range(x)):
            yield [part[i] for i in perm]
                    


#partitions of int n into x parts: lists
#counted by: partition number p_x(n)
def int_partition_eq(n, x):
    def int_partition_eq_max(n, x, m):
        if x == 0 and n == 0:
            yield []
        elif x <= 0 or n <= 0:
            pass
        else:
            for first in reversed(range(1, min(n - x + 1, m) + 1)): #reversed so that first elem decreases
                for rest in int_partition_eq_max(n - first, x - 1, first):
                    yield [first] + rest
    return int_partition_eq_max(n, x, n)

#partitions of int n into <= x parts: sets
def int_partition_le(n, x):
    for y in range(x + 1):
        yield from int_partition_eq(n, y)

def int_partitions(n):
    return int_partition_le(n, n)

#ordered partitions of int n into x parts: lists
#counted by: n - 1 choose n - x
def int_compositions(n, x):
    if x == 0 and n == 0:
        yield []
    elif x <= 0 or n <= 0:
        pass
    else:
        for first in range(1, n - x + 2):
            for rest in int_compositions(n - first, x - 1):
                yield [first] + rest

#ordered partitions of int n <= x parts: lists
#counted by: sum n - 1 choose n - x
def int_compositions_le(n, x):
    for y in range(1, x + 1):
        yield from int_compositions(n, y)










































