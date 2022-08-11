from algebra import base
from algebra import matricies
import functools
import heapq
import itertools
import random


@functools.cache
def PolyOver(ring, var = "λ"):
    assert issubclass(ring, base.Ring)

    class Poly(base.Ring):
        @classmethod
        def test_values(cls):
            values = ring.test_values()
            return [cls([7 * c, 2]) for c in values[:2]] + [cls([2, -1, 3 * c]) for c in values[:2]]
    
        @classmethod
        def typestr(cls):
            return f"{ring}[{var}]"
        @classmethod
        def int(cls, n):
            assert type(n) == int
            return cls([n])
        @classmethod
        def convert(cls, x):
            try:
                return super().convert(x)
            except NotImplementedError:
                pass
            if isinstance(x, ring):
                return cls([x])
            return cls.convert(ring.convert(x))
        @classmethod
        def var(cls):
            return cls([0, 1])

        @classmethod
        def lagrange_interp(cls, points):
            points = tuple((ring.convert(inp), ring.convert(out)) for inp, out in points)
            n = len(points)
            Mat = matricies.MatrixOver(ring)
            M = Mat(n, n, [[points[r][0] ** c for c in range(n)] for r in range(n)])
            v = Mat(n, 1, [[points[r][1]] for r in range(n)])
            #need to solve Mx = v where x is the vector of coefficients
            x = M.col_solve(v)
            assert x.rows == n and x.cols == 1
            ans = cls([x[i, 0] for i in range(n)])
            for inp, out in points:
                assert ans(inp) == out

            return ans
        def __init__(self, rep):
            rep = tuple(ring.convert(x) for x in rep)
            while (False if len(rep) == 0 else rep[-1] == 0): #remove trailing zeros
                rep = rep[:-1]
            self.rep = rep
            
        def __str__(self):
            bo = "("
            bc = ")"
            
            def strterm(c, p):
                if p == 0:
                    if str(c) == "0":
                        return ""
                    ans = base.add_mulbrac(str(c))
                    if len(ans) > 0:
                        if not ans[0] in {"+", "-"}:
                            ans = "+" + ans
                    return ans
                if p == 1:
                    var_term_str = var
                else:
                    var_term_str = f"{var}^{p}"
                if c == 0:
                    return ""
                elif str(c) == "1":
                    return "+" + var_term_str
                elif str(c) == "-1":
                    return "-" + var_term_str
                else:
                    return strterm(c, 0) + var_term_str
                    
            if self == 0:
                return "0"
            ans = "".join(strterm(c, p) for p, c in enumerate(self.rep))
            assert len(ans) != 0
            if ans[0] == "+":
                ans = ans[1:]
            return ans
        
        def __repr__(self):
            return f"{type(self)}({str(self)})"

        def __getitem__(self, p):
            assert type(p) == int
            assert p >= 0
            if 0 <= p < len(self.rep):
                return self.rep[p]
            else:
                return ring.int(0)
                
        def hash(self):
            return hash(self.rep)
        def equal(self, other):
            assert (cls := type(self)) == type(other)
            return self.rep == other.rep
        def add(self, other):
            assert (cls := type(self)) == type(other)
            rep = [ring.int(0) for _ in range(max(len(self.rep), len(other.rep)))]
            for i, x in enumerate(self.rep):
                rep[i] += x
            for i, x in enumerate(other.rep):
                rep[i] += x
            return cls(rep)
        def neg(self):
            return type(self)(tuple(-x for x in self.rep))
        def mul(self, other):
            assert (cls := type(self)) == type(other)
            rep = [ring.int(0) for _ in range(len(self.rep) + len(other.rep))]
            for i, x in enumerate(self.rep):
                for j, y in enumerate(other.rep):
                    rep[i + j] += x * y
            return cls(rep)
        def degree(self):
            #zero polynomial -> degree -1
            if self == 0:
                return None
            else:
                return len(self.rep) - 1
        def __call__(self, obj):
            ans = self.rep[0]
            obj_pow = obj
            for c in self.rep[1:]:
                ans += c * obj_pow
                obj_pow *= obj
            return ans
        def derivative(self):
            return type(self)(i * self.rep[i] for i in range(1, len(self.rep)))


    if issubclass(ring, base.IntegralDomain):
        class Poly(Poly, base.IntegralDomain):
            #we implement divmod early as the algorithm can be used for exact division too
            def divmod(self, other):
                assert (cls := type(self)) == type(other)
                if other == 0:
                    raise ZeroDivisionError
                quo = cls([])
                rem = self
                while True:
                    if rem == 0:
                        break
                    elif rem.degree() < other.degree():
                        break
                    new_quo = cls([ring.int(0)] * (rem.degree() - other.degree()) + [rem[rem.degree()] / other[other.degree()]])
                    rem = rem - new_quo * other
                    quo = quo + new_quo
                return quo, rem
            def exactdiv(self, other):
                assert (cls := type(self)) == type(other)
                q, r = self.divmod(other) #this only raises NotDivisibleError in the base ring if they are not divisible here too
                if r == 0:
                    return q
                else:
                    raise base.NotDivisibleError()

    if issubclass(ring, base.UniqueFactorizationDomain):
        class Poly(Poly, base.UniqueFactorizationDomain):
            @classmethod
            def test_axioms(cls, test):
                super().test_axioms(test)
                for f in cls.test_values():
                    m, p = f.factor_primitive()
                    test.assertEqual(m * p, f)
                        
            @classmethod
            def can_factor(cls):
                if issubclass(ring, base.Field):
                    if ring.is_fraction_field():
                        return PolyOver(ring.fraction_ring).can_factor()
                return ring.can_factor() and ring.has_finite_units()

            def factor_primitive(self):
                m = ring.gcd_list(self.rep)
                return m, self / m
            
            def factor(self):
                import sympy
                if ring == base.ZZ:
                    x = sympy.symbols("x")
                    s_poly = sum((x ** p * sympy.Integer(c.rep) for p, c in enumerate(self.rep)), sympy.Integer(0))
                    s_poly = sympy.poly(s_poly, x, domain = sympy.ZZ)
                    mult, factors = s_poly.factor_list()
                    mult = ring(int(mult))
                    mult = mult.factor()
                    mult = type(self).Factorization(self.convert(mult.unit), {self.convert(f) : p for f, p in mult.powers.items()})
                    
                    factors = {type(self)(tuple(reversed(tuple(ring(int(coeff)) for coeff in poly.all_coeffs())))) : power for poly, power in factors}
                    factors = type(self).Factorization(1, factors)
                    return mult * factors

                #kronekers method over general UFDs
                #need to add squarefree factorization to speed this up
                #also finding rational roots

                print(self)
                
                m, f = self.factor_primitive()

                def factor_primitive_kronekers(f):
                    assert f.degree() != 0
                    if f.degree() == 1:
                        return type(self).Factorization(type(self).int(1), {f : 1})

##                    g = f.sqfree_part()
##                    print(f, g)
                    
                    def possible_factors(f):
                        #we are reduced to factoring m in ring and f in cls
                        #self = m * f(x)
                        k = f.degree() // 2
                        #if f(x) factors as g(x)h(x) then wolg we need only check g(x) up to degree k
                        def yield_points():
                            yield 0
                            n = 0
                            while True:
                                n += 1
                                yield n
                                yield -n
                        yield_points = yield_points()
                        extra = 10 * (k + 1) #can try tweaking this
                        assert extra >= 0
                        
                        points = [next(yield_points) for _ in range(k + 1 + extra)]                
                        f_at_points = []
                        for x in points:
                            fx = f(ring.int(x))
                            if fx == 0:
                                #found an irreducible factor already
                                yield type(self).var() - type(self).int(x)
                                assert False #we know this factor works so should never be here
                            f_at_points.append(fx)
                        points = [ring.int(x) for x in points]
                        
                        f_at_points_divs = [list(f(x).factor().divisors()) for x in points]

                        idxs = heapq.nsmallest(k+1, range(k + 1 + extra), key = lambda i : len(f_at_points_divs[i]))
                        points = [points[i] for i in idxs]
                        f_at_points_divs = [f_at_points_divs[i] for i in idxs]
                        f_at_points_divs = [sum([[u * d.expand() for d in divs] for u in ring.all_units()], []) for divs in f_at_points_divs]
                        all_divs = list(itertools.product(*f_at_points_divs))
                        random.shuffle(all_divs)
                        for divs in all_divs:
                            #divs ranges over a finite set of possible values of g(x) at points
                            #we can reconstruct g from this using lagrange interpolation
                            try:
                                g = type(self).lagrange_interp(zip(points, divs))
                            except matricies.NoSolution:
                                pass
                            else:
                                if g.degree() > 0:
                                    yield g

                    for g in possible_factors(f):
                        try:
                            new_f = f / g
                        except base.NotDivisibleError:
                            pass
                        else:
                            return factor_primitive_kronekers(g) * factor_primitive_kronekers(new_f)
                            
                    return type(self).Factorization(type(self).int(1), {f : 1})

                m_factored = m.factor()
                m_factored = type(self).Factorization(type(self).convert(m_factored.unit), {type(self).convert(irr) : p for irr, p in m_factored.powers.items()})
                return m_factored * factor_primitive_kronekers(f)
                

    if issubclass(ring, base.Field):
        class Poly(Poly, base.EuclideanDomain):
            @classmethod
            def test_axioms(cls, test):
                super().test_axioms(test)
                if ring.is_fraction_field(): 
                    for f in cls.test_values():
                        m, p = f.factor_primitive_field()
                        test.assertEqual(cls.convert(m) * cls.convert(p), f)

            @classmethod
            def convert(cls, x):
                try:
                    return super().convert(x)
                except NotImplementedError:
                    pass
                if ring.is_fraction_field():
                    PolyOverRing = PolyOver(ring.fraction_ring)
                    if isinstance(x, PolyOverRing):
                        return cls([ring.convert(c) for c in x.rep])
                return cls.convert(ring.convert(x))

            def norm(self):
                return self.degree()
            #divmod is defined over an integral domain

            @classmethod
            def can_factor(cls):
                if ring.is_fraction_field():
                    return PolyOver(ring.fraction_ring).can_factor()
                return super().can_factor()

            #note the difference with factor_primitive
            #here the primitive part is a polynomial over the base ring of the fractions
            def factor_primitive_field(self):
                if ring.is_fraction_field():
                    PolyOverRing = PolyOver(ring.fraction_ring)
                    d = ring.fraction_ring.lcm_list([c.d for c in self.rep])
                    prim = self * d #over the field
                    prim = PolyOverRing(c.n / c.d for c in prim.rep) #over the base ring
                    n, prim = prim.factor_primitive() #factor_primitive over base ring
                    return ring(n, d), prim
                raise Exception("Can't factor_primitive_field over a field which is not a field of fractions")

            def factor(self):
                if ring.is_fraction_field():
                    PolyOverRing = PolyOver(ring.fraction_ring)                    
                    m, prim = self.factor_primitive_field()
                    prim_factored = prim.factor()
                    return type(self).Factorization(type(self).convert(prim_factored.unit) * type(self).convert(m), {type(self).convert(irr) : p for irr, p in prim_factored.powers.items()})                    
                raise NotImplementedError()

            def sqfree_part(self):
                return self // type(self).gcd(self, self.derivative())

    return Poly






if False:
    class PolyType():
        pass


    @functools.cache
    def PolyOver(ring):
        class Poly(PolyType, basic.Euclidean):
            @basic.cached_classmethod
            def cyclotomic(cls, n):
                if ring is basic.ZZ:
                    n = basic.ZZ.convert(n)
                    # x ** n - 1 = prod_d|n cyc(d)
                    return (cls.var() ** n.rep - 1) // (cls.product(cls.cyclotomic(d) for d in n.divisors() if d != n))
                else:
                    return cls([int(x) for x in PolyOver(basic.ZZ).cyclotomic(n)])
                
            def roots(self):
                #return all roots in the base ring
                rs = []
                for f, p in self.factor().powers().items():
                    if f.degree() == 1:
                        a, b = f.rep
                        #a + bx = 0
                        if b.is_unit():
                            for _ in range(p):
                                rs.append(-a / b)
                return rs


            def monic(self):
                return self * self[self.degree()].recip()
            def derivative(self):
                return type(self)(i * self.rep[i] for i in range(1, len(self.rep)))
            def sqfree_part(self):
                return self // type(self).gcd(self, self.derivative())
            def reversed(self):
                return type(self)(reversed(self.rep))

            def over_frac_field(self):
                frac_field = basic.FractionsOver(ring)
                return PolyOver(frac_field)([frac_field.convert(coeff) for coeff in self.rep])
            
        Poly.ring = ring

        if ring is basic.ZZ:
            class Poly(Poly):
                def factor(self):
                    import sympy
                    x = sympy.symbols("x")
                    s_poly = sum((x ** p * sympy.Integer(c.rep) for p, c in enumerate(self.rep)), sympy.Integer(0))
                    s_poly = sympy.poly(s_poly, x, domain = sympy.ZZ)
                    mult, factors = s_poly.factor_list()
                    mult = ring(int(mult))
                    mult = mult.factor().convert(type(self))
                    factors = {type(self)(tuple(reversed(tuple(ring(int(coeff)) for coeff in poly.all_coeffs())))) : power for poly, power in factors}
                    factors = basic.FactorsOver(type(self))(1, factors)
                    return mult * factors


        if issubclass(ring, basic.Quotient):
            if ring.ring is basic.ZZ:
                if ring.n.is_irreducible():
                    prime = int(ring.n)
                    class Poly(Poly):
                        def factor(self):
                            import sympy
                            x = sympy.symbols("x")
                            s_poly = sum((x ** p * sympy.Integer(c.rep) for p, c in enumerate(self.rep)), sympy.Integer(0))
                            s_poly = sympy.poly(s_poly, x, modulus = prime)                        
                            mult, factors = s_poly.factor_list()
                            mult = ring(int(mult))                        
                            mult = mult.factor().convert(type(self))
                            factors = {type(self)(tuple(reversed(tuple(ring(int(coeff)) for coeff in poly.all_coeffs())))) : power for poly, power in factors}
                            factors = basic.FactorsOver(type(self))(1, factors)
                            return mult * factors

                

        if issubclass(ring, basic.Factorable):
            class Poly(Poly):
                def factor_primitive(self):
                    if self == 0:
                        return ring.int(1), self
                    g = ring.gcd_list([c for c in self.rep if c != 0])
                    poly = self / g
                    u = poly.rep[-1].assoc_recip()
                    return g * u.recip(), poly * u

                def primitive(self):
                    return self.factor_primitive()[1]


        if issubclass(ring, basic.FieldOfFractions):
            class Poly(Poly):
                def factor_primitive(self):
                    if self == 0:
                        return ring.int(1), self
                    ds = [coeff.d for coeff in self.rep]
                    m = ring.ring.lcm_list(ds)
                    poly = PolyOver(ring.ring)([(coeff * m).n for coeff in self.rep])
                    g, poly = poly.factor_primitive()
                    return ring(g, m), poly

                def primitive(self):
                    return self.factor_primitive()[1]
                    
                def factor(self):
                    mult, prim = self.factor_primitive()
                    return basic.FactorsOver(type(self))(mult, {key.over_frac_field() : item for key, item in prim.factor().items()})
            
        return Poly





    def test():
        D = basic.Modulo(basic.ZZ(2))
        Po = PolyOver(D)
        x = Po.var()

        p = x ** 2 + 1 - x ** 2 - 1

        c, fs = p.factor_list()
        print(c)
        for f in fs:
            print(f)









        
