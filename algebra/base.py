import warnings
from fractions import Fraction as Frac
import functools

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

    def __str__(self):
        ans = self.typestr()
        if ans is None:
            return super().__str__()
        else:
            return ans

class MathsSet(metaclass = NumType):
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
                        
    @cached_classmethod
    def FactorizationCls(ring):
        if not ring.can_factor():
            raise NotImplementedError() #set cls.can_factor to return true if factorization is possible
        
        class Factorization():
            @classmethod
            def gcd(cls, items):
                keys, items = cls.make_compatible(items)
                return cls(1, {key : min(factorization[key] for factorization in items) for key in keys})

            @classmethod
            def lcm(cls, items):
                keys, items = cls.make_compatible(items)
                return cls(1, {key : max(factorization[key] for factorization in items) for key in keys})


            def __init__(self, unit, powers):
                unit = ring.convert(unit)
                assert type(powers) == dict
                powers = {ring.convert(f) : p for f, p in powers.items()}
                for p in powers.values():
                    assert type(p) == int
                self.unit = unit
                self.powers = powers

            def __eq__(self, other):
                if (cls := type(self)) == type(other):
                    return self.expand() == other.expand()
                return False
                
            def __repr__(self):
                return f"{ring}_Factorization({self.unit} " + " ".join([f"{f}^{p}" for f, p in self.powers.items()]) + ")"
            
            def expand(self):
                return self.unit * ring.product([f ** p for f, p in self.powers.items()])
            
            def __mul__(self, other):
                if (cls := type(self)) == type(other):
                    unit = self.unit * other.unit
                    powers = {f : p for f, p in self.powers.items()}
                    for g, q in other.powers.items():
                        for f in powers:
                            if ring.are_associate(f, g):
                                #we want to replace g^q with f^q
                                #there will be a unit factor induced by this
                                #the unit factor is (g/f)^q
                                powers[f] += q
                                unit *= (g / f) ** q
                                break
                        else:
                            powers[g] = q
                    return cls(unit, powers)
                return NotImplemented
        return Factorization

    @classmethod
    def Factorization(cls, *args, **kwargs):
        return cls.FactorizationCls()(*args, **kwargs)

    def factor(self):
        raise NotImplementedError()

    def is_irreducible(self):
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
    def floordiv(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
    def mod(self, other):
        assert (cls := type(self)) == type(other)
        raise NotImplementedError()
    
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
        return self.__floordiv__(other), self.__mod__(other)
    
    @cached_classmethod
    def QuotientRing(ring, n):
        n = ring.convert(n)
        
        class QuotientRing(Ring):
            @classmethod
            def test_values(cls):
                return [cls(x) for x in ring.test_values()]

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
        if self.rep == 0:
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
    
##    def __add__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(self.rep + other.rep)
##        return super().__add__(other)
##    def __radd__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(other.rep + self.rep)
##        return super().__radd__(other)
##    def __sub__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(self.rep - other.rep)
##        return super().__sub__(other)
##    def __rsub__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(other.rep - self.rep)
##        return super().__rsub__(other)
##    def __mul__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(self.rep * other.rep)
##        return super().__mul__(other)
##    def __rmul__(self, other):
##        if (cls := type(self)) == type(other):
##            return cls(other.rep * self.rep)
##        return super().__rmul__(other)
##    def __divmod__(self, other):
##        if (cls := type(self)) == type(other):
##            q, r = divmod(self.rep, other.rep)
##            return cls(q), cls(r)
##        return super().__floordiv__(other), super().__mod__(other)

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


##
##@functools.cache
##def Modulo(n):
##    ring = type(n)
##    assert issubclass(ring, EuclideanDomain)
##    
##    
##        
##    QuotientRing.ring = ring
##    QuotientRing.n = n
##
##    if issubclass(ring, UniqueFactorizationDomain):
##        if n.is_irreducible():
##            class QuotientRing(QuotientRing, Field):
##                pass
##
##    return QuotientRing










    

