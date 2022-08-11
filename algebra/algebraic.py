from algebra import base
from algebra import polynomials
from fractions import Fraction as Frac
import sympy
import pyalgebra
import functools

ZZ = base.ZZ
QQ = base.ZZ.FractionField
PolyZZ = polynomials.PolyOver(ZZ)
PolyQQ = polynomials.PolyOver(QQ)

#eventual plan: replace the sympy stuff here with my own algorithms which work over arbitrary polynomial rings
def poly_to_sympy(poly, x):
    assert isinstance(poly, PolyZZ)
    poly = [int(c) for c in poly.rep]
    s_poly = sum([c * x ** k for k, c in enumerate(poly)], sympy.Integer(0))
    s_poly = sympy.Poly(s_poly, x, domain = "ZZ")
    return s_poly

def sympy_to_poly(poly, x):
    assert poly.gen == x
    return PolyZZ([int(c) for c in reversed(poly.all_coeffs())])
    
    
def root_sum_poly(p, q):
    x, z = sympy.symbols("x, z")
    p = poly_to_sympy(p, x)
    q = poly_to_sympy(q, x)
    r = sympy.Poly(q.compose(sympy.Poly(z - x)), x, domain = "ZZ[z]")
    sum_poly = sympy.Poly(sympy.resultant(p, r)).sqf_part()
    return sympy_to_poly(sum_poly, z)

def root_prod_poly(p, q):
    x, t = sympy.symbols("x, t")
    p = poly_to_sympy(p, t)
    q = poly_to_sympy(q, x)  
    r = sympy.Poly(q.homogenize(t), t, domain = "ZZ[x]")
    #x ** q.degree() * q(t * x ** -1)
    prod_poly = sympy.Poly(sympy.resultant(p, r), x, domain = "ZZ").sqf_part()
    return sympy_to_poly(prod_poly, x)

def eval_at_frac(poly, x):
    assert isinstance(poly, PolyZZ)
    assert type(x) == Frac
    ans = Frac(0, 1)
    x_pow = Frac(1, 1)
    for p, c in enumerate(poly.rep):
        c = int(c)
        ans += c * x_pow
        x_pow *= x
    return ans



#the idea now is to define standalone classes _RealRep and _ComplexRep representing algebraic numbers
#then define a class called Algebraic which uses an instance of either Frac, _RealRep, or _ComplexRep as its representative

class _RealRep():
    def __init__(self, poly, a, b):
        assert isinstance(poly, PolyZZ)
        assert poly.degree() >= 2
        assert poly.is_irreducible()
        assert type(a) == Frac
        assert type(b) == Frac
        assert a < b
        at_a = eval_at_frac(poly, a)
        at_b = eval_at_frac(poly, b)
        assert at_a != 0
        assert at_b != 0
        assert (at_a < 0) != (at_b < 0)
        self.poly = poly
        self.at_a = at_a
        self.at_b = at_b
        self.a = a
        self.b = b

    def __str__(self):
        dp = 3
        while self.b - self.a > 10 ** (-dp-1):
            self.refine()
        return "≈" + str(round(float((self.a + self.b) / 2), dp))
    
    def __add__(self, other):
        if type(other) == _RealRep:
            polys = list(root_sum_poly(self.poly, other.poly).factor().powers.keys())
            while True:
                counts = [pyalgebra.polynomials.count_real_roots([int(c) for c in poly.rep], self.a + other.a, self.b + other.b) for poly in polys]
                non_zero = set([i for i, c in enumerate(counts) if c >= 1])
                counts = [c for i, c in enumerate(counts) if i in non_zero]
                polys = [poly for i, poly in enumerate(polys) if i in non_zero]
                assert len(polys) != 0
                if len(polys) == 1:
                    if counts[0] == 1:
                        return _RealRep(polys[0], self.a + other.a, self.b + other.b)
                self.refine()
                other.refine()
        if type(other) == Frac:            
            _, poly = self.poly(PolyZZ.int(other.denominator) * PolyZZ.var() - PolyZZ.int(other.numerator)).factor_primitive()
            return _RealRep(poly, self.a + other, self.b + other)
        return NotImplemented
    
    def __radd__(self, other):
        if type(other) == Frac:
            return self + other
        return NotImplemented
    
    def __neg__(self):            
        return _RealRep(self.poly(-PolyZZ.var()), -self.b, -self.a)
    
    def __mul__(self, other):
        if type(other) == _RealRep:
            polys = list(root_prod_poly(self.poly, other.poly).factor().powers.keys())
            while True:
                points = (self.a * other.a, self.a * other.b, self.b * other.a, self.b * other.b)
                counts = [pyalgebra.polynomials.count_real_roots([int(c) for c in poly.rep], min(points), max(points)) for poly in polys]
                non_zero = set([i for i, c in enumerate(counts) if c >= 1])
                counts = [c for i, c in enumerate(counts) if i in non_zero]
                polys = [poly for i, poly in enumerate(polys) if i in non_zero]
                assert len(polys) != 0
                if len(polys) == 1:
                    if counts[0] == 1:
                        return _RealRep(polys[0], self.a + other.a, self.b + other.b)
                self.refine()
                other.refine()
                
                try:
                    ans = Real(poly, min(points), max(points))
                except BadRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
##        if type(other) == QQ:
##            if other == 0:
##                return other
##            elif other < 0:
##                return (-other) * (-self)
##            else:
##                assert other > 0
##                other = Frac(int(other.n), int(other.d))
##                poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [0, 1 / other]))
##                return Real(PolyZZ(poly), self.a * other, self.b * other)
        return NotImplemented
    
##    def __rmul__(self, other):
####        if type(other) == QQ:
####            return self * other
##        return NotImplemented

    def refine(self):
        #update our internal variables to get a better approximation
        m = Frac(self.a + self.b, 2)
        at_m = eval_at_frac(self.poly, m)
        assert at_m != 0 #contradics irreducibility of self.poly which has degree >= 2
        if (at_m < 0) == (self.at_a < 0):
            self.a = m
            self.at_a = at_m
        else:
            self.b = m
            self.at_b = at_m
    

class _ComplexRep():
    def __init__(self):
        pass



if False:
    @functools.lru_cache()
    def primitive_factors(poly):
        assert type(poly) in {PolyZZ, PolyQQ}
        mult, factor_powers = poly.factor_powers()
        return {f.primitive() : k for f, k in factor_powers.items()}

    #should only be generated by the real_roots function
    #this means two real roots are equal iff they are actually the same instance
    class _Real():
        def __str__(self):
            dp = 3
            while self.b - self.a > 10 ** (-dp-1):
                self.refine()
            return "≈" + str(round(float((self.a + self.b) / 2), dp))
        
        def __add__(self, other):
            if type(other) == _Real:
                poly = PolyZZ(root_sum_poly(self.poly, other.poly))
                while True:
                    try:
                        ans = Real(poly, self.a + other.a, self.b + other.b)
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
            if type(other) == QQ:
                other = Frac(int(other.n), int(other.d))
                poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [-other, 1]))
                return Real(PolyZZ(poly), self.a + other, self.b + other)
            return NotImplemented
        def __radd__(self, other):
            if type(other) == QQ:
                return self + other
            return NotImplemented
        def __neg__(self):
            poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [0, -1]))
            return Real(PolyZZ(poly), -self.b, -self.a)
        def __mul__(self, other):
            if type(other) == _Real:
                poly = PolyZZ(root_prod_poly(self.poly, other.poly))
                while True:
                    points = (self.a * other.a, self.a * other.b, self.b * other.a, self.b * other.b)
                    try:
                        ans = Real(poly, min(points), max(points))
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
            if type(other) == QQ:
                if other == 0:
                    return other
                elif other < 0:
                    return (-other) * (-self)
                else:
                    assert other > 0
                    other = Frac(int(other.n), int(other.d))
                    poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [0, 1 / other]))
                    return Real(PolyZZ(poly), self.a * other, self.b * other)
            return NotImplemented
        def __rmul__(self, other):
            if type(other) == QQ:
                return self * other
            return NotImplemented
        def recip(self):
            poly = PolyZZ(tuple(reversed(self.poly)))
            while self.a == 0 or self.b == 0 or (self.a < 0) != (self.b < 0):
                self.refine()
            if self.a < 0:
                assert self.b < 0
                return -(-self).recip()
            assert self.a > 0 and self.b > 0
            while True:
                try:
                    ans = Real(poly, 1/self.b, 1/self.a)
                except BadRoot:
                    self.refine()
                else:
                    return ans
        
        def __lt__(self, other):
            raise NotImplementedError()
        def __int__(self):
            raise Exception("Real root is not an integer")
        def __float__(self):
            while self.b - self.a > 2 ** -64:
                self.refine()
            return float((self.a + self.b) / 2)
        def floor(self):
            raise NotImplementedError()

        def refine(self):
            m = Frac(self.a + self.b, 2)
            at_m = pyalgebra.polynomials.Poly.evaluate(self.poly, m)
            assert at_m != 0 #contradics irreducibility of self.poly which has degree >= 2
            if (at_m < 0) == (self.at_a < 0):
                self.a = m
                self.at_a = at_m
            else:
                self.b = m
                self.at_b = at_m


    class BadRoot(Exception):
        pass
        
    def Real(poly, a, b):
        def overlap(r):
            if type(r) == QQ:
                r = Frac(int(r.n), int(r.d))
                return a < r < b
            else:
                assert type(r) == _Real
                return a < r.b and r.a < b
        def contains(r):
            if type(r) == QQ:
                r = Frac(int(r.n), int(r.d))
                return a < r < b
            else:
                assert type(r) == _Real
                return a < r.a < b and a < r.b < b
        
        if type(a) == int:
            a = Frac(a)
        if type(b) == int:
            b = Frac(b)
        assert type(poly) in {PolyZZ, PolyQQ}
        assert type(a) == Frac
        assert type(b) == Frac
        assert a < b    
        roots = list(rat_roots(poly)) + list(real_roots(poly))
        while True:            
            roots = [r for r in roots if overlap(r)]
            contained = [r for r in roots if contains(r)]
            if len(contained) >= 2:
                raise BadRoot(f"{poly} has multiple real roots between {a} and {b}")
            if len(roots) == 0:
                raise BadRoot(f"{poly} has no real roots between {a} and {b}")
            if len(roots) == 1:
                return roots[0]
            for r in roots:
                if type(r) == _Real:
                    r.refine()
        
        

            
    #should only be generated by the imag_roots function
    #this means two complex roots are equal iff they are actually the same instance
    class _Complex():
        def __init__(self, poly, a, b, c, d):
            _validate_poly(poly)
            assert type(a) == Frac
            assert type(b) == Frac
            assert type(c) == Frac
            assert type(d) == Frac
            assert a < b
            assert c < d
            assert pyalgebra.polynomials.count_complex_roots(poly, a, b, c, d) == 1
            self.a = a
            self.b = b
            self.c = c
            self.d = d
            self.poly = poly

        def __str__(self):
            dp = 3
            while self.error > 10 ** (-dp-1):
                self.refine()
            return "≈" + str(complex(round(float((self.a + self.b) / 2), dp),
                                     round(float((self.c + self.d) / 2), dp))).replace("j", "i").replace("(", "").replace(")", "")
            
        @property
        def error(self):
            return max(self.b - self.a, self.d - self.c)

        def __add__(self, other):
            if type(other) == _Complex:
                poly = PolyZZ(root_sum_poly(self.poly, other.poly))
                while True:
                    try:
                        ans = Complex(poly, self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
                                
            elif type(other) == _Real:
                poly = PolyZZ(root_sum_poly(self.poly, other.poly))
                while True:
                    try:
                        ans = Complex(poly, self.a + other.a, self.b + other.b, self.c, self.d)
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
                
            if type(other) == QQ:
                other = Frac(int(other.n), int(other.d))
                poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [-other, 1]))
                return type(self)(tuple(poly), self.a + other, self.b + other, self.c, self.d)
            return NotImplemented
        def __radd__(self, other):
            if type(other) in {_Real, QQ}:
                return self + other
            return NotImplemented
        def __neg__(self):
            poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [0, -1]))
            return Complex(PolyZZ(poly), -self.b, -self.a, -self.d, -self.c)
        def __mul__(self, other):
    ##        def c_mul(z, w): #z and w are pairs of fractions representing complex numbers
    ##            return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

            if type(other) == _Complex:
                poly = PolyZZ(root_prod_poly(self.poly, other.poly))
                while True:
                    points_re = []
                    points_im = []
                    for w in [(other.a, other.c), (other.a, other.d), (other.b, other.c), (other.b, other.d)]:
                        for z in [(self.a, self.c), (self.a, self.d), (self.b, self.c), (self.b, self.d)]:
                            points_re.append(z[0] * w[0] - z[1] * w[1])
                            points_im.append(z[0] * w[1] + z[1] * w[0])
                    try:
                        ans = Complex(poly, min(points_re), max(points_re), min(points_im), max(points_im))
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
                
            if type(other) == _Real:
                poly = PolyZZ(root_prod_poly(self.poly, other.poly))
                while True:
                    points_re = []
                    points_im = []
                    for x in [other.a, other.b]:
                        for z in [(self.a, self.c), (self.a, self.d), (self.b, self.c), (self.b, self.d)]:
                            points_re.append(x * z[0])
                            points_im.append(x * z[1])
                    try:
                        ans = Complex(poly, min(points_re), max(points_re), min(points_im), max(points_im))
                    except BadRoot:
                        self.refine()
                        other.refine()
                    else:
                        return ans
                    
            if type(other) == QQ:
                if other == 0:
                    return other
                elif other < 0:
                    return (-other) * (-self)
                else:
                    assert other > 0
                    other = Frac(int(other.n), int(other.d))
                    poly = pyalgebra.polynomials.Poly.primitive(pyalgebra.polynomials.Poly.compose(self.poly, [0, 1 / other]))
                    return Complex(PolyZZ(poly), self.a * other, self.b * other, self.c * other, self.d * other)
            return NotImplemented
        def __rmul__(self, other):
            if type(other) in {_Real, QQ}:
                return self * other
            return NotImplemented
        def recip(self):
            poly = PolyZZ(tuple(reversed(self.poly)))

            #z = the root represented by self
            #a = the center of the approximation
            #eps = the radius of the approximation

            #result used to get region which contains the reciprocal:
            #if |z - a| < eps < a then |1/z - 1/a| < eps / 2|a|^2

            while True:            
                eps = max(self.b - self.a, self.d - self.c) / 2
                #now |1/z-1/a| <= eps / 2|a|^2
                
                w_re = Frac(self.a + self.b, 2)
                w_im = Frac(self.c + self.d, 2)
                w_mag_sq = w_re ** 2 + w_im ** 2

                delta = eps / (2 * w_mag_sq)
                w_recip_re = w_re / w_mag_sq
                w_recip_im = -w_im / w_mag_sq
                
                try:
                    ans = Complex(poly, w_recip_re - delta, w_recip_re + delta, w_recip_im - delta, w_recip_im + delta)
                except BadRoot as e:
                    self.refine()
                else:
                    return ans

            

            
        
        def __lt__(self, other):
            raise Exception("Cant compare ordering on complex numbers")
        def __int__(self):
            raise Exception("Complex root is not an integer")
        def __float__(self):
            raise Exception("Complex root is not real")
        def __complex__(self):
            while self.error > 2 ** -64:
                self.refine()
            return complex((self.a + self.b) / 2, (self.c + self.d) / 2)
        def floor(self):
            raise Exception("Can't take the floor of a complex root")

        def refine(self):
            def pick(box1, box2):
                if pyalgebra.polynomials.count_complex_roots(self.poly, *box1) == 1:
                    return box1
                else:
                    return box2
            a, b, c, d = self.a, self.b, self.c, self.d
            if b - a > d - c:
                m = Frac(a + b, 2)
                try:
                    box = pick((a, m, c, d), (m, b, c, d))
                except pyalgebra.polynomials.BoundaryRoot:
                    m = Frac(2 * a + b, 3)
                    box = pick((a, m, c, d), (m, b, c, d))
            else:
                m = Frac(c + d, 2)
                try:
                    box = pick((a, b, c, m), (a, b, m, d))
                except pyalgebra.polynomials.BoundaryRoot:
                    m = Frac(2 * c + d, 3)
                    box = pick((a, b, c, m), (a, b, m, d))
            self.a, self.b, self.c, self.d = box


    def Complex(poly, a, b, c, d):
        def overlap(r):
            if type(r) == QQ:
                r = Frac(int(r.n), int(r.d))
                return a < r < b and c < 0 < d
            elif type(r) == _Real:
                return a < r.b and r.a < b and c < 0 < d
            else:
                assert type(r) == _Complex
                return a < r.b and r.a < b and c < r.d and r.c < d
        def contains(r):
            if type(r) == QQ:
                r = Frac(int(r.n), int(r.d))
                return a < r < b and c < 0 < d
            elif type(r) == _Real:
                return a < r.a < b and a < r.b < b and c < 0 < d and c < 0 < d
            else:
                assert type(r) == _Complex
                return a < r.a < b and a < r.b < b and c < r.c < d and c < r.d < d
            
        if type(a) == int: a = Frac(a)
        if type(b) == int: b = Frac(b)
        if type(c) == int: c = Frac(c)
        if type(d) == int: d = Frac(d)
        assert type(poly) in {PolyZZ, PolyQQ}
        assert type(a) == Frac
        assert type(b) == Frac
        assert type(c) == Frac
        assert type(d) == Frac
        assert a < b
        assert c < d
        roots = list(rat_roots(poly)) + list(real_roots(poly)) + list(imag_roots(poly))
        while True:
            roots = [r for r in roots if overlap(r)]
            contained = [r for r in roots if contains(r)]
            if len(contained) >= 2:
                raise BadRoot(f"{poly} has multiple imaginary roots in {(a, b, c, d)}")
            if len(roots) == 0:
                raise BadRoot(f"{poly} has no imaginary roots in {(a, b, c, d)}")
            if len(roots) == 1:
                return roots[0]
            for r in roots:
                if type(r) in {_Real, _Complex}:
                    r.refine()



    class Algebraic(base.Field, base.Real):
        @classmethod
        def int(cls, n):
            return cls(QQ.int(n))

        def __init__(self, rep):
            assert type(rep) in [QQ, _Real, _Complex]
            self.rep = rep
        def __str__(self):
            return str(self.rep)

        def is_rat(self):
            return type(self.rep) == QQ
        def is_real(self):
            return type(self.rep) == _Real
        def is_complex(self):
            return type(self.rep) == _Complex
        
        def hash(self):
            return hash(self.rep)
        def equal(self, other):
            if (cls := type(self)) == type(other):
                if type(self.rep) == type(other.rep):
                    return self.rep == other.rep
            return False
            
        def add(self, other):
            assert (cls := type(self)) == type(other)
            return cls(self.rep + other.rep)
        def neg(self):
            return type(self)(-self.rep)
        def mul(self, other):
            assert (cls := type(self)) == type(other)
            return cls(self.rep * other.rep)
        def recip(self):
            return type(self)(self.rep.recip())
        def lt(self, other):
            assert (cls := type(self)) == type(other)
            return self.rep < other.rep
        def __int__(self):
            return int(self.rep)
        def __float__(self):
            return float(self.rep)
        def __complex__(self):
            if type(self.rep) in [QQ, _Real]:
                return complex(float(self.rep), 0)
            else:
                return complex(self.rep)
        
        def floor(self):
            return self.rep.floor()
        def min_poly(self):
            if self.is_rat():
                return (PolyQQ.var() - self.rep).primitive()
            elif self.is_real():
                return PolyZZ(self.rep.poly)
            elif self.is_complex():
                return PolyZZ(self.rep.poly)
        def degree(self):
            return self.min_poly().degree()





    def split_into_factors(find_roots):
        def new_f(poly):
            assert type(poly) in {PolyZZ, PolyQQ}
            poly = poly.primitive() #unique associate
            assert type(poly) == PolyZZ
            for p, k in primitive_factors(poly).items():
                for _ in range(k):
                    yield from find_roots(p)
        return new_f

    @split_into_factors
    @lambda f : lambda p : tuple(f(p))
    def rat_roots(poly):
        if poly.degree() == 1:
            a, b = poly[0], poly[1]
            yield QQ(-a, b)

    @split_into_factors
    @functools.cache
    @lambda f : lambda p : tuple(f(p))
    def real_roots(poly):
        if poly.degree() >= 2:
            poly_tuple = tuple(int(x) for x in poly.rep)
            for a, b in  pyalgebra.polynomials.isolate_real_roots(poly_tuple):
                yield _Real(poly_tuple, Frac(a), Frac(b))
        
    @split_into_factors
    @functools.cache
    @lambda f : lambda p : tuple(f(p))
    def imag_roots(poly):
        if poly.degree() >= 2:
            poly_tuple = tuple(int(x) for x in poly.rep)
            for a, b, c, d in pyalgebra.polynomials.isolate_imag_roots(poly_tuple):
                yield _Complex(poly_tuple, Frac(a), Frac(b), Frac(c), Frac(d))
        
        

    def roots(poly, *, rat = False, real = False, imag = False):
        if rat == real == imag == False:
            rat, real, imag = True, True, True
            
        if rat:
            for x in rat_roots(poly):
                yield Algebraic(x)
        if real:
            for x in real_roots(poly):
                yield Algebraic(x)
        if imag:
            for x in imag_roots(poly):
                yield Algebraic(x)
                



    def test():
        import numpy
        x = PolyZZ.var()
        
        p = x ** 7 + x ** 6 - 1
        rs = list(roots(p))

        for r in rs:
            print(r)

        print(Algebraic.sum(rs))


        print("done")

    ##    a = Algebraic.int(5) / Algebraic.int(6)
    ##    b = Algebraic.int(7) / Algebraic.int(8)
    ##
    ##    print(a)
    ##    print(b)
    ##    print(a + b)
    ##    print(a * b)
    











def roots(f):
    assert isinstance(f, PolyZZ)
    if f == 0:
        raise Exception()
    
    fs = f.factor().powers
    print(fs)



































