import math
from algebra import basic
import fractions
import itertools
import sympy
import functools


def factors(n):
    assert n >= 1
    return sympy.factorint(n)

def is_prime(n):
    return factors(n) == {n : 1}





def gcd(x, y):
    while y != 0:
        x, y = y, x % y
    return x


def xgcd(x, y):
    pa, a, pb, b = 1, 0, 0, 1
    while y != 0:
        q, r = divmod(x, y)
        a, pa = pa - q * a, a
        b, pb = pb - q * b, b
        x, y = y, r
    #g, a, b
    return x, pa, pb

def coprime(n):
    return set(x for x in range(1, n + 1) if gcd(x, n) == 1)

def phi(n):
    return len(coprime(n))

def prim_roots(n):
    mults = coprime(n)
    roots = set(x for x in mults if set(pow(x, p, n) for p in range(n)) == mults)
    assert phi(len(mults)) == len(roots) or len(roots) == 0
    return roots



def jacobi(n, k):
    assert(k > 0 and k % 2 == 1)
    n = n % k
    t = 1
    while n != 0:
        while n % 2 == 0:
            n = n // 2
            r = k % 8
            if r == 3 or r == 5:
                t = -t
        n, k = k, n
        if n % 4 == 3 and k % 4 == 3:
            t = -t
        n = n % k
    if k == 1:
        return t
    else:
        return 0

def find_jacobi(k, j):
    for x in range(k):
        if jacobi(x, k) == j:
            return x
    raise Exception()

def isqrt(n):
    assert n >= 0
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

#find sqrt mod p for prime p
def psqrt(n, p):
    assert is_prime(p)
    assert jacobi(n, p) == 1

    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1

    z = find_jacobi(p, -1)

    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)
    while True:
##        assert pow(c, 2 ** (M - 1), p) == p - 1
##        assert pow(t, 2 ** (M - 1), p) == 1
##        assert pow(R, 2, p) == (t * n) % p
        if t == 0:
            return 0
        elif t == 1:
            return R
        tp = t
        for i in range(M):
            if tp == 1:
                break
            tp = pow(tp, 2, p)
        else:
            raise Exception()
##        assert pow(t, 2 ** i, p) == 1
##        assert pow(t, 2 ** (i - 1), p) != 1
        b = pow(c, 2 ** (M - i - 1), p)
        M = i
        c = pow(b, 2, p)
        t = (t * pow(b, 2, p)) % p
        R = (R * b) % p
        





def two_squares(n):
    Gaussian = basic.ComplexOver(basic.Int)
    
    if n == 0:
        return 0, 0
    powers = factors(n)
    two_pow = powers.get(2, 0)
    one_four = {p : powers[p] for p in powers if p % 4 == 1}
    three_four = {p : powers[p] for p in powers if p % 4 == 3}

    one_primes = [p for p in one_four]

    for prime, power in three_four.items():
        if power % 2 == 1:
            return

    # p = 2 then 2 = (1 + i)(1 - i)
    # p = 1 mod 4 then p = aa* for some irreducible a in Z[i]
    # p = 3 mod 4 then p is irredcible in Z[i]

    one_four_factors = {}
    for p in one_primes:
        #prime = 1 mod 4
        #so -1 is a QR mod p
        #so there exists n such that p | n^2 + 1

        n = psqrt(-1, p)
        
        #now, p | (n + i) * (n - i)

        p_gaussian = Gaussian(p, 0)
        a = Gaussian(n, 1)
        z = Gaussian.gcd(a, p_gaussian)
        assert z * z.conjugate() == p_gaussian
        one_four_factors[p] = z
    
    const = Gaussian.product([Gaussian(1, 1) for _ in range(two_pow)] + [Gaussian(prime, 0) ** (power // 2) for prime, power in three_four.items()])
    done = set([])

    for powers in itertools.product(*[range(one_four[p] + 1) for p in one_primes]):
        one_facts = [[(one_four_factors[p] if powers[idx] <= i else one_four_factors[p].conjugate()) for i in range(one_four[p])] for idx, p in enumerate(one_primes)]

        ans = const
        for fact in one_facts:
            ans *= Gaussian.product(fact)

        ans = tuple(sorted([abs(int(ans.a)), abs(int(ans.b))]))
        if not ans in done:
            done.add(ans)
            yield int(ans[0]), int(ans[1])

    total = 1
    for p, k in one_four.items():
        total *= (k + 1)
    if total % 2 == 1:
        total += 1
    total //= 2

    assert total == len(done)
                



def four_squares(n):
    if n == 0:
        return 0, 0, 0, 0
    H = basic.QuaternionOver(basic.Int)
    fs = factors(n)

    two_pow = fs.get(2, 0)
    one_pows = {p : fs[p] for p in fs if p % 4 == 1}
    three_pows = {p : fs[p] for p in fs if p % 4 == 3}

    ans = H(1, 1, 0, 0) ** two_pow
    for p in one_pows:
        #always exists for p = 1 mod 4
        a, b = next(two_squares(p))
        ans *= H(a, b, 0, 0) ** one_pows[p]

    def descent(p, x, y, z, w):
        rp = x ** 2 + y ** 2 + z ** 2 + w ** 2
        assert rp % p == 0
        r = rp // p
        if r == 1:
            return x, y, z, w
        elif r % 2 == 0:
            #in this case, rp is even, so wlog x=y mod 2 or z=w mod 2
            #then (x+y)/2, (x-y)/2, (z+w)/2, (z-w)/2 do the trick
            if (x - z) % 2 == 0:
                y, z = z, y
            elif (x - w) % 2 == 0:
                y, z, w = w, y, z
            assert (x - y) % 2 == 0
            assert (w - z) % 2 == 0
            a, b, c, d = (x + y) // 2, (x - y) // 2, (z + w) // 2, (z - w) // 2
            return descent(p, a, b, c, d)
        else:
            #reduce x, y, z, w mod r to between -r/2 < a, b, c, d < r/2
            r_off = (r - 1) // 2
            a = (x + r_off) % r - r_off
            b = (y + r_off) % r - r_off
            c = (z + r_off) % r - r_off
            d = (w + r_off) % r - r_off
            #now a^2 + b^2 + c^2 + d^2 = rr' for some r' != 0
            #(a^2 + b^2 + c^2 + d^2)(x^2 + y^2 + z^2 + w^2) = r^2r'p
            quart = H(a, b, c, d) * H(x, y, z, w).conjugate()
            x, y, z, w = quart.a, quart.b, quart.c, quart.d
            assert x % r == 0
            assert y % r == 0
            assert z % r == 0
            assert w % r == 0
            x, y, z, w = x // r, y // r, z // r, w // r
            return descent(p, x, y, z, w)
        
    for p in three_pows:
        # -1 is NOT a quadratic residue mod p
        # it is therefore possible to find a in Z such that a is QR & a + 1 is NQR
        # then a = x^2 and a + 1 = -y^2 for some x, y in Z
        # now x^2 + y^2 + 1 = 0 mod p
        # so x^2 + y^2 + 1 = rp for some r in Z
        # can then apply desent on this expression until we get a sum of 4 squares equal to p

        #first find a
        q = (p - 1) // 2
        for a in range(1, p - 1):
            if pow(a, q, p) == 1 and pow(a + 1, q, p) == p - 1:
                break
        else:
            raise Exception("Maths / Python is bork")

        x = psqrt(a, p)
        y = psqrt(-a-1, p)

        assert (pow(x, 2, p) + pow(y, 2, p) + 1) % p == 0
        #now perform descent on x, y, 1, 0 untill we get x^2 + y^2 + z^2 + w^2 = p
        z = H(*descent(p, x, y, 1, 0))
        assert z.norm() == p
        ans *= z ** three_pows[p]

    return tuple(sorted([abs(int(ans.a)), abs(int(ans.b)), abs(int(ans.c)), abs(int(ans.d))]))
        





def pells(d):
    #sols to x^2 - d * y^2
    assert d >= 1
    assert isqrt(d) ** 2 != d
    from algebra import algebraic
    yield 1, 0
    rd = algebraic.Algebraic.root(d, 2)
    #find fundamental solution
    for a, p, q in rd.continued_fraction():
        if p ** 2 - d * q ** 2 == 1:
            break
    #generate there rest from there
    x, y = p, q
    while True:
        yield x, y
        x, y = x * p + y * q * d, x * q + y * p #(x + yr(d))(p + qr(d))
    






def run():
##    from algebra import algebraic
##    
##    x = algebraic.Algebraic.root(37, 2)
##    for a, p, q in x.continued_fraction():
##        print(a, p, q)

##    for x, y in pells(37):
##        print(x, y)

    
##    for n in range(100):
##        print(n, four_squares(n))        

    print(jacobi(1, 7))
    print(jacobi(2, 7))
    print(jacobi(3, 7))
    print(jacobi(4, 7))
    print(jacobi(5, 7))
    print(jacobi(6, 7))
    


















    

