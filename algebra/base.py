import warnings
from fractions import Fraction as Frac
import functools
import itertools

def cached_classmethod(func):
    @functools.cache
    def do_func(*args, **kwargs):
        return func(*args, **kwargs)
    @classmethod
    def new_func(*args, **kwargs):
        return do_func(*args, **kwargs)
    return new_func

class NumType(type):   
    def __init__(self, name, bases, dct):
        super().__init__(name, bases, dct)

    def __new__(*args, **kwargs):
        x = type.__new__(*args, **kwargs)
        x.init_cls()
        return x

    def __str__(self):
        ans = self.typestr()
        if ans is None:
            return super().__str__()
        else:
            return ans

class MathsSet(metaclass = NumType):
    @classmethod
    def init_cls(cls):
        #use to create classes specific to each ring.
        #for example to create a unique Factorization class for each new type of ring
        pass

    @classmethod
    def test_values(cls):
        warnings.warn(f"No test values provided for {cls}", RuntimeWarning)
        return []
    
    @classmethod
    def test_axioms(cls, test):
        pass
        
    @classmethod
    def typestr(cls):
        return None

    @classmethod
    def convert(cls, x):
        #try to convert an arbitrary object x to one of type cls
        #should only be done via cannonical maps
        if type(x) == cls:
            return x
        raise NotImplementedError()




class AbGroup(MathsSet):
    @classmethod
    def test_axioms(cls, test):
        values = cls.test_values()
        for a in values:
            test.assertEqual(a, a + cls.zero())
        for a in values:
            for b in values:
                test.assertEqual(a + (-b), a - b)
                test.assertEqual(a + b, b + a)
    
    @classmethod
    def zero(cls):
        raise NotImplementedError()
        
    @classmethod
    def sum(cls, vals):
        if len(vals) == 0:
            return cls.zero()
        else:
            ans = vals[0]
            for val in vals[1:]:
                ans += val
            return ans
        
    def hash(self):
        raise NotImplementedError()
    def equal(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
    def add(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
    def neg(self):
        raise NotImplementedError()
    
    def __hash__(self):
        return self.hash()
    def __eq__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.equal(other)
    def __add__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.add(other)
    def __radd__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return other.add(self)
    def __sub__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.add(other.neg())
    def __rsub__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.neg().add(other)
    def __neg__(self):
        return self.neg()


class NotDivisibleError(ZeroDivisionError):
    pass



class Ring(AbGroup):
    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        values = cls.test_values()
        for a in values:
            test.assertEqual(a, a * 1)
            test.assertEqual(a ** 3, a * a * a)
            with test.assertRaises(ZeroDivisionError):
                a / cls.int(0)
        for a in values:
            for b in values:
                test.assertEqual(a * b, b * a)
        for a in values:
            for b in values:
                for c in values:
                    test.assertEqual(a * (b + c), a * b + a * c)
                    test.assertEqual((a + b) * c, a * c + b * c)

                
    @classmethod
    def int(cls, n):
        raise NotImplementedError()

    @classmethod
    def zero(cls):
        return cls.int(0)

    @classmethod
    def product(cls, vals):
        ans = cls.int(1)
        for val in vals:
            ans *= val
        return ans
    
    @classmethod
    def convert(cls, x):
        try:
            return super().convert(x)
        except NotImplementedError:
            pass
        if type(x) == int:
            return cls.int(x)
        else:
            raise NotImplementedError()
        
    def mul(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
    #called by the syntax | or implemented in terms of // and % in an ED
    def exactdiv(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
        raise NotDivisibleError(f"Cant divide {self} by {other}")

    def recip(self):
        return 1 / self
    
    def is_unit(self):
        if self == 0:
            return False
        try:
            1 / self
            return True
        except NotDivisibleError:
            return False

    def __mul__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.mul(other)
    def __rmul__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return other.mul(self)

    #implement a/b as exact division
    #note that this is not compatible with builtin integer division, so builtin integers should
    #not be used as Ring elements which rely on the implementation of Ring. e.g. as matrix entries
    #they should instead first be converted to ZZ type
    def __truediv__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.exactdiv(other)
    def __rtruediv__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return other.exactdiv(self)

    def __pow__(self, other):
        cls = type(self)
        if type(other) == int:
            if other < 0:
                return (self ** (-other)).recip()
            elif other == 0:
                return cls.int(1)
            elif other == 1:
                return self
            elif other == 2:
                return self * self
            else:
                q, r = divmod(other, 2)
                return (self ** q) ** 2 * self ** r
        return NotImplemented


class IntegralDomain(Ring):
    @classmethod
    def are_associate(cls, a, b):
        try:
            a / b
            b / a
        except NotDivisibleError:
            return False
        else:
            return True
    
    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        values = cls.test_values()
        for a in values:
            for b in values:
                if a != 0 and b != 0:
                    test.assertNotEqual(a * b, 0)





class UniqueFactorizationDomain(IntegralDomain):
    @classmethod
    def init_cls(cls):
        super().init_cls()
        ring = cls
        class Factorization():
            @classmethod
            def match(cls, facts):
                facts = list(facts)
                for a in facts:
                    assert isinstance(a, cls)
                units = [a.unit for a in facts]
                powers = {}
                for i, a in enumerate(facts):
                    for g, q in a.powers.items():
                        for f in powers:
                            if ring.are_associate(f, g):
                                #we want to replace g^q with f^q
                                #there will be a unit factor induced by this
                                #the unit factor is (g/f)^q
                                powers[f][i] = q
                                units[i] *= (g / f) ** q
                                break
                        else:
                            powers[g] = [0] * len(facts)
                            powers[g][i] = q
                return units, powers

            #it is nice how product, gcd, lcm corespond to sum, min, max.
            @classmethod
            def product(cls, facts):
                units, powers = cls.match(facts)
                return cls(ring.product(units), {f : sum(ps) for f, ps in powers.items()})
            @classmethod
            def gcd(cls, facts):
                units, powers = cls.match(facts)
                return cls(ring.int(1), {f : min(ps) for f, ps in powers.items()})
            @classmethod
            def lcm(cls, facts):
                units, powers = cls.match(facts)
                return cls(ring.int(1), {f : max(ps) for f, ps in powers.items()})

            def __init__(self, unit, powers):
                unit = ring.convert(unit)
                assert type(powers) == dict
                for p in powers.values():
                    assert type(p) == int and p >= 0
                powers = {ring.convert(f) : p for f, p in powers.items() if p != 0}
                self.unit = unit
                self.powers = powers
            def __eq__(self, other):
                if (cls := type(self)) == type(other):
                    return self.expand() == other.expand()
                return False
            def __repr__(self):
                return f"{ring}_Factorization({self.unit}; " + " * ".join([f"{f}^{p}" for f, p in self.powers.items()]) + ")"
            def __mul__(self, other):
                if (cls := type(self)) == type(other):
                    return cls.product([self, other])
                return NotImplemented
            
            def expand(self):
                return self.unit * ring.product([f ** p for f, p in self.powers.items()])

            #list all irreducible factors repeated by multiplicity
            def list(self):
                factors = []
                for f, p in self.powers.items():
                    for _ in range(p):
                        factors.append(f)
                return factors
            #yield all divisors up to associates
            def divisors(self):
                fs = list(self.powers.keys())
                powers = [list(range(self.powers[f] + 1)) for f in fs]
                for rep in itertools.product(*powers):
                    yield type(self)(1, dict(zip(fs, rep)))
            
        cls.Factorization = Factorization
        
    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        if cls.can_factor():
            values = cls.test_values()
            for a in values:
                if a == 0:
                    test.assertEqual(a.factor(), None)
                else:
                    factored = a.factor()
                    test.assertEqual(factored.expand(), a)
                    fs = list(factored.powers.keys())
                    for i in range(len(fs)):
                        for j in range(i + 1, len(fs)):
                            test.assertFalse(cls.are_associate(fs[i], fs[j]))

            for a in values:
                for b in values:
                    if a != 0 and b != 0:
                        test.assertEqual(a.factor() * b.factor(), (a * b).factor())

    @classmethod
    def can_factor(cls):
        return False

    def factor(self):
        raise NotImplementedError()

    def is_irreducible(self):
        fp = self.factor().powers
        if len(fp) == 1:
            if fp[list(fp.keys())[0]] == 1:
                return True
        return False


class PrincipalIdealDomain(UniqueFactorizationDomain):
    pass


class EuclideanDomain(PrincipalIdealDomain):
    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        test.assertEqual(cls.int(0).norm(), None)
        values = cls.test_values()
        for a in values:
            if a != 0:
                test.assertTrue(type(a.norm()) == int)
                test.assertTrue(a.norm() >= 0)
        for a in values:
            for b in values:
                if b != 0:
                    q, r = divmod(a, b)
                    if r != 0:
                        test.assertTrue(r.norm() < b.norm())
        for a in values:
            for b in values:
                for c in values:
                    if a != 0 and b != 0 and c != 0:
                        g = cls.gcd_list([a, b, c])
                        test.assertTrue(a % g == b % g == c % g == 0)                            
                        m = cls.lcm_list([a, b, c])
                        test.assertTrue(m % a == m % b == m % c == 0)
                        g, coeffs = cls.xgcd_list([a, b, c])
                        x, y, z = coeffs
                        test.assertEqual(x * a + y * b + z * c, g)

                        #does euclidean gcd/lcm agree with the gcd/lcm computed via counting irreducible powers
                        if cls.can_factor():
                            test.assertTrue(cls.are_associate(cls.gcd_list([a, b, c]), cls.Factorization.gcd([a.factor(), b.factor(), c.factor()]).expand()))
                            test.assertTrue(cls.are_associate(cls.lcm_list([a, b, c]), cls.Factorization.lcm([a.factor(), b.factor(), c.factor()]).expand()))

    @classmethod
    def gcd(cls, x, y):
        if x == 0 or y == 0:
            raise ZeroDivisionError()
        while y != cls.int(0):
            x, y = y, x % y
        return x

    @classmethod
    def gcd_list(cls, elems):
        assert len(elems) >= 1
        if len(elems) == 1:
            g = elems[0]
        elif len(elems) == 2:
            g = cls.gcd(*elems)
        else:
            i = len(elems) // 2
            g = cls.gcd(cls.gcd_list(elems[:i]), cls.gcd_list(elems[i:]))
        return g

    @classmethod
    def lcm(cls, x, y):
        if x == 0 or y == 0:
            raise ZeroDivisionError()
        x = cls.convert(x)
        y = cls.convert(y)
        g = cls.gcd(x, y)
        m = x * y
        q, r = divmod(m, g)
        assert r == 0
        return q

    @classmethod
    def lcm_list(cls, elems):
        assert len(elems) >= 1
        if len(elems) == 1:
            return elems[0]
        elif len(elems) == 2:
            return cls.lcm(*elems)
        else:
            i = len(elems) // 2
            return cls.lcm(cls.lcm_list(elems[:i]), cls.lcm_list(elems[i:]))

    @classmethod
    def xgcd(cls, x, y):
        pa, a, pb, b = cls.int(1), cls.int(0), cls.int(0), cls.int(1)
        while y != cls.int(0):
            q, r = divmod(x, y)
            a, pa = pa - q * a, a
            b, pb = pb - q * b, b
            x, y = y, r
        return x, pa, pb

    @classmethod
    def xgcd_list(cls, elems):
        if all(elem == 0 for elem in elems):
            return cls.int(0), [cls.int(0) for elem in elems]
        elems = list(elems)
        assert len(elems) >= 1
        if len(elems) == 1:
            return elems[0], [cls.int(1)]
        elif len(elems) == 2:
            g, x, y = cls.xgcd(elems[0], elems[1])
            return g, [x, y]
        else:
            n = len(elems) // 2
            g1, coeffs1 = cls.xgcd_list(elems[:n])
            g2, coeffs2 = cls.xgcd_list(elems[n:])
            g, x, y = cls.xgcd(g1, g2)
            return g, [x * c for c in coeffs1] + [y * c for c in coeffs2]

    def norm(self) -> int:
        raise NotImplementedError()
    #implement either divmod OR floordiv and mod
    def divmod(self, other):
        assert (cls := type(self)) == type(other)
        return self.floordiv(other), self.mod(other)
    def floordiv(self, other):
        assert (cls := type(self)) == type(other)
        return self.divmod(other)[0]
    def mod(self, other):
        assert (cls := type(self)) == type(other)
        return self.divmod(other)[1]
    
    def exactdiv(self, other):
        assert (cls := type(self)) == type(other)
        q, r = divmod(self, other)
        if r != 0:
            raise NotDivisibleError()
        else:
            return q

    def __floordiv__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.floordiv(other)
    def __rfloordiv__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return other.floordiv(self)
    def __mod__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.mod(other)
    def __rmod__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return other.mod(self)
    def __divmod__(self, other):
        try:
            other = type(self).convert(other)
        except NotImplementedError:
            return NotImplemented
        else:
            return self.divmod(other)
    
    @cached_classmethod
    def QuotientRing(ring, n):
        n = ring.convert(n)
        
        class QuotientRing(Ring):
            @classmethod
            def test_values(cls):
                return list(set(cls(x) for x in ring.test_values()))

            @classmethod
            def test_axioms(cls, test):
                super().test_axioms(test)
                values = [z.rep for z in cls.test_values()]
                for a in values:
                    test.assertEqual(cls(-a), -cls(a))
                for a in values:
                    for b in values:
                        test.assertEqual(a, a + cls(n) * b)
        
            @classmethod
            def typestr(cls):
                return f"{ring}/<{n}>"
            @classmethod
            def int(cls, x):
                assert type(x) == int
                return cls(ring.int(x))
            @classmethod
            def convert(cls, x):
                try:
                    return super().convert(x)
                except NotImplementedError:
                    pass
                if isinstance(x, ring):
                    return cls(x)
                else:
                    return cls.convert(ring.convert(x))

            def __init__(self, rep):
                rep = ring.convert(rep)
                self.rep = rep % n
                
            def __str__(self):
                return str(self.rep)
            def __repr__(self):
                return f"Mod({self.rep}, {n})"

            def hash(self):
                #rep is not unique in general
                #is unique for integers & polynomials over fields tho
                return 0
            def equal(self, other):
                assert (cls := type(self)) == type(other)
                return self.rep - other.rep % n == 0
            def add(self, other):
                assert (cls := type(self)) == type(other)
                return cls(self.rep + other.rep)
            def neg(self):
                return type(self)(-self.rep)
            def mul(self, other):
                assert (cls := type(self)) == type(other)
                return cls(self.rep * other.rep)

            def exactdiv(self, other):
                assert (cls := type(self)) == type(other)
                g, x, y = ring.xgcd(n, other.rep)
                #g = x*n + y*other
                if g != ring.int(1):
                    raise NotDivisibleError()
                #1 = y*other mod n
                return self * cls(y)

        if issubclass(ring, UniqueFactorizationDomain):
            if n.is_irreducible():
                class QuotientRing(QuotientRing, Field):
                    pass

        return QuotientRing


class Field(EuclideanDomain):
    @classmethod
    def convert(cls, x):
        try:
            return super().convert(x)
        except NotImplementedError:
            pass
        if type(x) == Frac:
            return cls.int(x.numerator) / cls.int(x.denominator)
        else:
            raise NotImplementedError()

    @classmethod
    def can_factor(cls):
        return True

    #should implement either recip or exactdiv to get a field
    def exactdiv(self, other):
        assert (cls := type(self)) == type(other)
        if other == 0:
            raise ZeroDivisionError()
        else:
            return self * other.recip()
        
    def norm(self):
        if self == 0:
            return None
        else:
            return 0
    def floordiv(self, other):
        return self / other
    def mod(self, other):
        assert (cls := type(self)) == type(other)
        return 0
    def factor(self):
        if self == 0:
            return None
        else:
            return self.Factorization(self, {})



class ZZ(EuclideanDomain):
    @classmethod
    def typestr(cls):
        return "ℤ"
    
    @classmethod
    def test_values(cls):
        return [cls(n) for n in [-48, -31, -1, 0, 1, 2, 6, 10, 13, 15, 168]]

    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        values = [z.rep for z in cls.test_values()]
        for a in values:
            test.assertEqual(cls(-a), -cls(a))
        for a in values:
            for b in values:
                test.assertEqual(cls(a + b), cls(a) + cls(b))
                test.assertEqual(cls(a - b), cls(a) - cls(b))
                test.assertEqual(cls(a * b), cls(a) * cls(b))

    @classmethod
    def can_factor(cls):
        return True
        
    @classmethod
    def int(cls, n):
        return cls(n)
    
    def __init__(self, rep):
        assert type(rep) == int
        self.rep = rep

    def __str__(self):
        return str(self.rep)
    def __repr__(self):
        return f"{type(self).typestr()}{repr(self.rep)}"

    def hash(self):
        return hash(self.rep)
    def equal(self, other):
        assert (cls := type(self)) == type(other)
        return self.rep == other.rep
    def add(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep + other.rep)
    def neg(self):
        return type(self)(-self.rep)
    def mul(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep * other.rep)
    def norm(self) -> int:
        if self == 0:
            return None
        else:
            return abs(self.rep)
    def floordiv(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep // other.rep)
    def mod(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep % other.rep)

    def factor(self):
        if self == 0:
            return None
        import sympy
        fs = sympy.factorint(self.rep)
        unit = 1
        if -1 in fs:
            unit = -1
            del fs[-1]
        powers = {type(self)(f) : p for f, p in fs.items()}
        unit = type(self)(unit)
        return type(self).Factorization(unit, powers)
    



class Gaussian(EuclideanDomain):
    @classmethod
    def typestr(cls):
        return "ℤi"
    
    @classmethod
    def test_values(cls):
        return [cls(1, 2), cls(3, 4), cls(-2, 5), cls(0, 0), cls(1, 0), cls(0, 1), cls(11, 0), cls(13, 0)]

    @classmethod
    def can_factor(cls):
        return True
        
    @classmethod
    def int(cls, n):
        return cls(n, 0)
    
    def __init__(self, a, b):
        assert type(a) == int
        assert type(b) == int
        self.a = a
        self.b = b

    def __str__(self):
        return f"{self.a} + {self.b}i"
    def __repr__(self):
        return f"{type(self).typestr()}({repr(self.a)}, {repr(self.b)})"

    def hash(self):
        return hash((self.a, self.b))
    def equal(self, other):
        assert (cls := type(self)) == type(other)
        return self.a == other.a and self.b == other.b
    def add(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.a + other.a, self.b + other.b)
    def neg(self):
        return type(self)(-self.a, -self.b)
    def mul(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.a * other.a - self.b * other.b, self.a * other.b + self.b * other.a)
    def norm(self) -> int:
        if self == 0:
            return None
        else:
            return self.a ** 2 + self.b ** 2
    def divmod(self, other):
        assert (cls := type(self)) == type(other)
        if other == 0:
            raise ZeroDivisionError()
        d = other.norm()
        x_prime = Frac(int(self.a) * int(other.a) + int(self.b) * int(other.b), int(d))
        y_prime = Frac(int(self.b) * int(other.a) - int(self.a) * int(other.b), int(d))
        x, y = round(x_prime), round(y_prime)
        q = cls(x, y)
        r = self - q * other
        return q, r

    def factor(self):
        if self == 0:
            return None
        
        import sympy
        cls = type(self)
        int_fs = sympy.factorint(self.norm())
        assert not -1 in int_fs

        def psqrt(n, p):
            assert sympy.ntheory.isprime(p)
            assert sympy.ntheory.legendre_symbol(n, p) == 1

            Q = p - 1
            S = 0
            while Q % 2 == 0:
                Q //= 2
                S += 1

            for z in range(p):
                if sympy.ntheory.legendre_symbol(z, p) == -1:
                    break
            else:
                assert False

            M = S
            c = pow(z, Q, p)
            t = pow(n, Q, p)
            R = pow(n, (Q + 1) // 2, p)
            while True:
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
                b = pow(c, 2 ** (M - i - 1), p)
                M = i
                c = pow(b, 2, p)
                t = (t * pow(b, 2, p)) % p
                R = (R * b) % p

        def find_two_sq_sum(p):
            #p is a prime congruent to 1 mod 4
            #find a, b such that a^2 + b^2 = p
            m = psqrt(-1, p)
            #we have p dividing m^2+1 = (m+i)(m-i)
            g = cls.gcd(cls(p, 0), cls(m, 1))
            a, b = g.a, g.b
            assert a ** 2 + b ** 2 == p
            return a, b
        
        powers = {}
        for p in int_fs:
            if p == 2:
                powers[cls(1, 1)] = int_fs[2]
            elif p % 4 == 3:
                assert int_fs[p] % 2 == 0
                powers[cls(p, 0)] = int_fs[p] // 2

        omf_part = (self // cls.Factorization(cls.int(1), powers).expand()) #one mod four part
        for p in int_fs:
            if p % 4 == 1:
                a, b = find_two_sq_sum(p)
                for d in [cls(a, b), cls(a, -b)]:
                    powers[d] = 0
                    while omf_part % d == 0:
                        powers[d] += 1
                        omf_part //= d
        i = cls(0, 1)
        assert omf_part in [1, i, -1, -i]
        unit = self / cls.Factorization(cls.int(1), powers).expand()
        return cls.Factorization(unit, powers)


class QQ(Field):
    @classmethod
    def typestr(cls):
        return "ℚ"
    
    @classmethod
    def test_values(cls):
        return [QQ(0, 1), QQ(1, 1), QQ(1, 2), QQ(-2, 3), QQ(-3, 4), QQ(5, 7)]

    @classmethod
    def test_axioms(cls, test):
        super().test_axioms(test)
        values = [z.rep for z in cls.test_values()]
        for a in values:
            test.assertEqual(cls(-a), -cls(a))
        for a in values:
            for b in values:
                test.assertEqual(cls(a + b), cls(a) + cls(b))
                test.assertEqual(cls(a - b), cls(a) - cls(b))
                test.assertEqual(cls(a * b), cls(a) * cls(b))
        
    @classmethod
    def int(cls, n):
        return cls(n, 1)
    
    def __init__(self, *rep):
        if len(rep) == 2:
            rep = Frac(rep[0], rep[1]),
        rep = rep[0]
        assert type(rep) == Frac
        self.rep = rep

    def __str__(self):
        return str(self.rep)
    def __repr__(self):
        return f"{type(self).typestr()}({self.rep.numerator}/{self.rep.denominator})"

    def hash(self):
        return hash(self.rep)
    def equal(self, other):
        assert (cls := type(self)) == type(other)
        return self.rep == other.rep
    def add(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep + other.rep)
    def neg(self):
        return type(self)(-self.rep)
    def mul(self, other):
        assert (cls := type(self)) == type(other)
        return cls(self.rep * other.rep)
    def recip(self):
        return type(self)(1 / self.rep)









    

