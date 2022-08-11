import itertools
from fractions import Fraction as Frac


def product(nums):
    ans = 1
    for n in nums:
        ans *= n
    return ans

def factors(n):
    import sympy
    fs = sympy.factorint(n)
    unit = 1
    if -1 in fs:
        unit = -1
        assert fs[-1] == 1
        del fs[-1]
    return unit, fs

def is_prime(n):
    assert type(n) == int and n >= 2
    _, f = factors(n)
    if len(f) == 1:
        p = next(iter(f.keys()))
        if f[p] == 1:
            return True
    return False



#p adic norm of x
def pnorm(p, x):
    if x == 0:
        return 0
    x = Frac(x)
    n, d = x.numerator, x.denominator
    k = 0
    while n % p == 0:
        k += 1
        n //= p
    while d % p == 0:
        k -= 1
        d //= p
    return Frac(1, p) ** k


#p adic expansion of p-adic integer x
def pexpand(p, x, maxdigs = 10):
    def get_digs(x, maxdigs):
        assert pnorm(p, x) <= 1
        if maxdigs == 0:
            return        
        x = Frac(x)
        n, d = x.numerator, x.denominator

        #x - m = (n - md) / d should have |x - m| < 1 (iff n == md mod 3), use xgcd for this
        g, e, f = xgcd(d, p)
        assert g == 1
        m = (e * n) % p       
        yield m % p
        if x - m != 0:
            yield from get_digs(Frac(x - m, p), maxdigs - 1)

    assert is_prime(p)
    assert type(x) == Frac

    digits = list(reversed(list(get_digs(x, maxdigs))))
    if p < 10:
        ans = "".join(str(m) for m in digits)
    else:
        ans = ",".join(str(m) for m in digits)
    if len(digits) == maxdigs:
        ans = "..." + ans
    return ans
    







def gcd(x, y):
    while y != 0:
        x, y = y, x % y
    return x

def gcd_list(elems):
    elems = list(elems)
    assert len(elems) >= 1
    if len(elems) == 1:
        return elems[0]
    elif len(elems) == 2:
        return gcd(*elems)
    else:
        i = len(elems) // 2
        return gcd(gcd_list(elems[:i]), gcd_list(elems[i:]))

def lcm(x, y):
    return (x * y) // gcd(x, y)

def lcm_list(elems):
    assert len(elems) >= 1
    if len(elems) == 1:
        return elems[0]
    elif len(elems) == 2:
        return lcm(*elems)
    else:
        i = len(elems) // 2
        return lcm(lcm_list(elems[:i]), lcm_list(elems[i:]))

def xgcd(x, y):
    pa, a, pb, b = 1, 0, 0, 1
    while y != 0:
        q, r = divmod(x, y)
        a, pa = pa - q * a, a
        b, pb = pb - q * b, b
        x, y = y, r
    return x, pa, pb

def xgcd_list(elems):
    if all(elem == 0 for elem in elems):
        return 0, [0 for elem in elems]
    elems = list(elems)
    assert len(elems) >= 1
    if len(elems) == 1:
        return elems[0], [1]
    elif len(elems) == 2:
        g, x, y = xgcd(elems[0], elems[1])
        return g, [x, y]
    else:
        n = len(elems) // 2
        g1, coeffs1 = xgcd_list(elems[:n])
        g2, coeffs2 = xgcd_list(elems[n:])
        g, x, y = xgcd(g1, g2)
        #x*g1 + y*g2 = g
        #x*(coeffs1[i]*elems[i]) + y*(coeffs2[i]*elems[i+n]) = g
        return g, [x * c for c in coeffs1] + [y * c for c in coeffs2]



def divisors(n):
    assert type(n) == int and n >= 1
    _, fs = factors(n)
    for prime_power_factors in itertools.product(*[[prime ** k for k in range(power + 1)] for prime, power in fs.items()]):
        yield product(prime_power_factors)


def phi(n):
    assert type(n) == int and n >= 1
    _, f = factors(n)
    return product((p - 1) * p ** (k - 1) for p, k in f.items())




def phi_inverse(n):
    #yield all k such that phi(k) == n

    def phi_factor(n, ex = set([])):
        if n == 1:
            yield []
            
        for d in divisors(n):
            p = d + 1
            if is_prime(p):
                if not p in ex:
                    #possible first bits
                    small_n = n // (p - 1)
                    m_pow = 0
                    while small_n % p ** (m_pow + 1) == 0:
                        m_pow += 1
                        
                    for k in range(m_pow + 1):
                        first = (p - 1) * p ** k                        
                        if first != 1:                    
                            for phi_f in phi_factor(n // first, ex | {p}):
                                yield [(p, k)] + phi_f

    found = set([])
    for phi_f in phi_factor(n):
        phi_inv = product(p ** (k + 1) for p, k in phi_f)
        if not phi_inv in found:
            found.add(phi_inv)
            yield phi_inv
            
        if phi_inv % 2 == 1:
            phi_inv *= 2
            if not phi_inv in found:
                found.add(phi_inv)
                yield phi_inv





























        

