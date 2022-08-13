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

class NoUniqueRoot(Exception):
    pass

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

    def __hash__(self):
        return hash(self.poly.degree())

    def __eq__(self, other):
        if type(other) == _RealRep:
            if self.b < other.a or other.b < self.a:
                return False
            return self - other == 0
        return False
    
    def __add__(self, other):
        if type(other) == _RealRep:
            poly = root_sum_poly(self.poly, other.poly)
            while True:
                try:
                    ans = RealRep(poly, self.a + other.a, self.b + other.b)
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
        if type(other) == Frac:
            _, poly = self.poly(PolyQQ([-QQ(other.numerator, other.denominator), 1])).factor_primitive_field()
            return _RealRep(poly, self.a + other, self.b + other)
        return NotImplemented
    def __radd__(self, other):
        if type(other) == Frac:
            return self + other
        return NotImplemented
    def __sub__(self, other):
        return self + (-other)
    def __rsub__(self, other):
        return (-self) + other
    def __neg__(self):
        return _RealRep(self.poly(-PolyZZ.var()), -self.b, -self.a)
    def __mul__(self, other):
        if type(other) == _RealRep:
            poly = root_prod_poly(self.poly, other.poly)
            while True:
                points = (self.a * other.a, self.a * other.b, self.b * other.a, self.b * other.b)
                try:
                    ans = RealRep(poly, min(points), max(points))
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
        if type(other) == Frac:
            if other == 0:
                return other
            elif other < 0:
                return (-other) * (-self)
            else:
                assert other > 0
                _, poly = self.poly(PolyQQ([0, QQ.int(other.denominator) / QQ.int(other.numerator)])).factor_primitive_field()
                return _RealRep(poly, self.a * other, self.b * other)
        return NotImplemented
    def __rmul__(self, other):
        if type(other) == Frac:
            return self * other
        return NotImplemented
    def __truediv__(self, other):
        if type(other) == Frac:
            return self * (1 / other)
        elif type(other) in {_RealRep}:
            return self * other.recip()
        return NotImplemented
    def __rtruediv__(self, other):
        if type(other) in {Frac, _RealRep}:
            return self.recip() * other
        return NotImplemented

    def recip(self):
        poly = PolyZZ(tuple(reversed(self.poly.rep)))
        while self.a == 0 or self.b == 0 or (self.a < 0) != (self.b < 0):
            self.refine()
        if self.a < 0:
            assert self.b < 0
            return -(-self).recip()
        assert self.a > 0 and self.b > 0
        while True:
            try:
                ans = RealRep(poly, 1/self.b, 1/self.a)
            except NoUniqueRoot:
                self.refine()
            else:
                return ans

    def __int__(self):
        raise Exception("Real root is not an integer")
    def __float__(self):
        while self.b - self.a > 2 ** -64:
            self.refine()
        return float((self.a + self.b) / 2)
    
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

def RealRep(poly, a, b):
    ra = pyalgebra.polynomials.rational_real_root(a)
    rb = pyalgebra.polynomials.rational_real_root(b)
    polys = list(poly.factor().powers.keys())
    roots = [pyalgebra.polynomials.real_roots([int(c) for c in poly.rep], a, b) for poly in polys]
    prs = list(zip(polys, roots))
    prs = [(poly, [r for r in roots if ra < r < rb]) for poly, roots in prs] #delete roots which lie outside the search range
    prs = [(poly, roots) for poly, roots in prs if len(roots) != 0] #delete polynomials with no possible roots
    if len(prs) == 0:
        raise Exception() #if there are _no_ roots then something has gone really wrong
    #if there is a unique root in the range, return it
    if len(prs) == 1:
        poly, roots = prs[0]
        if len(roots) == 1:
            root = roots[0]
            if poly.degree() == 1:
                r, s = poly.rep
                rat = -Frac(int(r), int(s))
                assert root.a <= rat <= root.b
                return rat
            else:
                return _RealRep(poly, root.a, root.b)
    raise NoUniqueRoot()
    


class _ComplexRep():
    def __init__(self, poly, a, b, c, d):
        assert poly.degree() >= 2
        assert poly.is_irreducible()
        assert type(a) == Frac
        assert type(b) == Frac
        assert type(c) == Frac
        assert type(d) == Frac
        assert a < b
        assert c < d
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.poly = poly

    @property
    def error(self):
        return max(self.b - self.a, self.d - self.c)

    def __hash__(self):
        return hash(self.poly.degree())
    
    def __eq__(self, other):
        if type(other) == _ComplexRep:
            if self.b < other.a or other.b < self.a or self.d < other.c or other.d < self.c:
                return False
            return self - other == 0
        return False

    def __str__(self):
        dp = 3
        while self.error > 10 ** (-dp-1):
            self.refine()
        return "≈" + str(complex(round(float((self.a + self.b) / 2), dp),
                                 round(float((self.c + self.d) / 2), dp))).replace("j", "i").replace("(", "").replace(")", "")

    def __int__(self):
        raise Exception("Complex root is not an integer")
    def __float__(self):
        raise Exception("Complex root is not real")
    def __complex__(self):
        while self.error > 2 ** -64:
            self.refine()
        return complex((self.a + self.b) / 2, (self.c + self.d) / 2)

    def __add__(self, other):
        if type(other) == _ComplexRep:
            poly = root_sum_poly(self.poly, other.poly)
            while True:
                try:
                    ans = ComplexRep(poly, self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
                            
        elif type(other) == _RealRep:
            poly = root_sum_poly(self.poly, other.poly)
            while True:
                try:
                    ans = ComplexRep(poly, self.a + other.a, self.b + other.b, self.c, self.d)
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
            
        if type(other) == Frac:
            _, poly = self.poly(PolyQQ([-QQ(other.numerator, other.denominator), 1])).factor_primitive_field()
            return _ComplexRep(poly, self.a + other, self.b + other, self.c, self.d)
        return NotImplemented
    def __radd__(self, other):
        if type(other) in {Frac, _RealRep}:
            return self + other
        return NotImplemented
    def __sub__(self, other):
        return self + (-other)
    def __rsub__(self, other):
        return (-self) + other
    def __neg__(self):
        return _ComplexRep(self.poly(-PolyZZ.var()), -self.b, -self.a, -self.d, -self.c)
    def __mul__(self, other):
        if type(other) == _ComplexRep:
            poly = root_prod_poly(self.poly, other.poly)
            while True:
                points_re = []
                points_im = []
                for w in [(other.a, other.c), (other.a, other.d), (other.b, other.c), (other.b, other.d)]:
                    for z in [(self.a, self.c), (self.a, self.d), (self.b, self.c), (self.b, self.d)]:
                        points_re.append(z[0] * w[0] - z[1] * w[1])
                        points_im.append(z[0] * w[1] + z[1] * w[0])
                try:
                    ans = ComplexRep(poly, min(points_re), max(points_re), min(points_im), max(points_im))
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
            
        if type(other) == _RealRep:
            poly = root_prod_poly(self.poly, other.poly)
            while True:
                points_re = []
                points_im = []
                for x in [other.a, other.b]:
                    for z in [(self.a, self.c), (self.a, self.d), (self.b, self.c), (self.b, self.d)]:
                        points_re.append(x * z[0])
                        points_im.append(x * z[1])
                try:
                    ans = ComplexRep(poly, min(points_re), max(points_re), min(points_im), max(points_im))
                except NoUniqueRoot:
                    self.refine()
                    other.refine()
                else:
                    return ans
                
        if type(other) == Frac:
            if other == 0:
                return other
            elif other < 0:
                return (-other) * (-self)
            else:
                assert other > 0
                _, poly = self.poly(PolyQQ([0, QQ.int(other.denominator) / QQ.int(other.numerator)])).factor_primitive_field()
                return _ComplexRep(poly, self.a * other, self.b * other, self.c * other, self.d * other)
        return NotImplemented
    def __rmul__(self, other):
        if type(other) in {Frac, _RealRep}:
            return self * other
        return NotImplemented
    def __truediv__(self, other):
        if type(other) == Frac:
            return self * (1 / other)
        elif type(other) in {_RealRep, _ComplexRep}:
            return self * other.recip()
        return NotImplemented
    def __rtruediv__(self, other):
        if type(other) in {Frac, _RealRep, _ComplexRep}:
            return self.recip() * other
        return NotImplemented

    def recip(self):
        poly = PolyZZ(reversed(self.poly.rep))

        #z = the root represented by self
        #a = the center of the approximation
        #eps = the radius of the approximation

        #result used to get region which contains the reciprocal:
        #if eps < |a|/2 is such that |z - a| < eps then
        #|1/z - 1/a| = |z - a| / |a||z| < |z - a| / |a|^2-|a|eps < 2 eps / |a|^2

        while True:            
            eps = (self.b - self.a + self.d - self.c) / 2 #max(self.b - self.a, self.d - self.c) / 2
            #now |1/z-1/a| <= eps / 2|a|^2
            
            w_re = Frac(self.a + self.b, 2)
            w_im = Frac(self.c + self.d, 2)
            w_mag_sq = w_re ** 2 + w_im ** 2

            if eps >= w_mag_sq / 2:
                self.refine()
                continue

            delta = 2 * eps / w_mag_sq
            w_recip_re = w_re / w_mag_sq
            w_recip_im = -w_im / w_mag_sq
            
            try:
                ans = ComplexRep(poly, w_recip_re - delta, w_recip_re + delta, w_recip_im - delta, w_recip_im + delta)
            except NoUniqueRoot as e:
                self.refine()
            else:
                return ans

    def refine(self):
        def pick(box1, box2):
            if pyalgebra.polynomials.count_complex_roots([int(c) for c in self.poly.rep], *box1) == 1:
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


def ComplexRep(poly, a, b, c, d):
    polys = list(poly.factor().powers.keys())
    try:
        counts = [pyalgebra.polynomials.count_complex_roots([int(c) for c in poly.rep], a, b, c, d) for poly in polys]
    except pyalgebra.polynomials.BoundaryRoot:
        raise NoUniqueRoot()
    count = sum(counts)
    if count == 0:
        raise Exception()
    elif count >= 2:
        raise NoUniqueRoot()

    for i in range(len(polys)):
        if counts[i] == 1:
            break
    else:
        assert False
    poly = polys[i] #this is the irreducible part of the input poly containing the root we are after
    #it remains to figure out which root of this polynomial we are after
    coeffs = [int(c) for c in poly.rep]

    #deal with rational roots
    if poly.degree() == 1:
        return -Frac(int(coeffs[0]), int(coeffs[1]))    
    
    roots = []
    for ar, br in pyalgebra.polynomials.isolate_real_roots(coeffs, a, b):
        roots.append(_RealRep(poly, ar, br))
    for ar, br, cr, dr in pyalgebra.polynomials.isolate_imag_roots(coeffs):
        roots.append(_ComplexRep(poly, ar, br, cr, dr))

    def overlap(r):
        if type(r) == _RealRep:
            return a < r.b and r.a < b and c < 0 < d
        else:
            assert type(r) == _ComplexRep
            return a < r.b and r.a < b and c < r.d and r.c < d

    while True:
        roots = [r for r in roots if overlap(r)]
        if len(roots) <= 1:
            break
        for r in roots:
            r.refine()
    assert len(roots) == 1
    return roots[0]
                     
def all_roots_rep(poly):
    assert isinstance(poly, PolyZZ)
    factors = poly.factor()
    for sub_poly, k in factors.powers.items():
        coeffs = [int(c) for c in sub_poly.rep]
        if sub_poly.degree() == 1:
            yield -Frac(int(coeffs[0]), int(coeffs[1]))
        else:
            for a, b in pyalgebra.polynomials.isolate_real_roots(coeffs):
                yield _RealRep(sub_poly, a, b)
            for a, b, c, d in pyalgebra.polynomials.isolate_imag_roots(coeffs):
                yield _ComplexRep(sub_poly, a, b, c, d)

Algebraic = base.ZZ.FractionField.AlgebraicClosure



























