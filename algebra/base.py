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
        raise NotImplementedError()
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


class ZZ(Ring):
    @classmethod
    def test_values(cls):
        return [cls(n) for n in [-48, -31, -1, 0, 1, 2, 13, 168]]

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
        return cls(n)
    
    def __init__(self, rep):
        assert type(rep) == int
        self.rep = rep

    def hash(self):
        return hash(self.rep)
    def equal(self, other):
        assert (cls := type(self)) == type(other)
        return self.rep == other.rep
    def add(self, other):
        return type(self)(self.rep + other.rep)
    def neg(self):
        return type(self)(-self.rep)
    def mul(self, other):
        assert (cls := type(self)) == type(other)
        return type(self)(self.rep * other.rep)
    def exactdiv(self, other):
        assert (cls := type(self)) == type(other)
        q, r = divmod(self.rep, other.rep)
        if r == 0:
            return cls(q)
        else:
            raise NotDivisibleError()













    

