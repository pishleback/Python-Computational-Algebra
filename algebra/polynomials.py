from algebra import base
import functools




@functools.cache
def PolyOver(ring, var = "λ"):
    assert issubclass(ring, base.Ring)

    class Poly(base.Ring):
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

        def factor(self):
            raise NotImplementedError()

    return Poly






if False:
    class PolyType():
        pass


    @functools.cache
    def PolyOver(ring):
        class Poly(PolyType, basic.Euclidean):
            @classmethod
            def lagrange_interp(cls, points):
                points = tuple((cls.convert(inp), cls.convert(out)) for inp, out in points)
                n = len(points)
                return cls.sum([cls.convert(points[i][1]) * cls.ring.product([points[i][0] - points[j][0] for j in range(n) if j != i]).recip() * cls.product([cls.var() - points[j][0] for j in range(n) if j != i]) for i in range(n)])

            @classmethod
            def lagrange_interp_mult(cls, points):
                #lagrange interp scaled by a constant so that coefficients are in the same base ring
                points = tuple((cls.convert(inp), cls.convert(out)) for inp, out in points)
                n = len(points)
                return cls.sum([cls.convert(points[i][1]) * cls.ring.product([cls.ring.product([points[k][0] - points[j][0] for j in range(n) if j != k]) for k in range(n) if k != i]) * cls.product([cls.var() - points[j][0] for j in range(n) if j != i]) for i in range(n)])

            @basic.cached_classmethod
            def cyclotomic(cls, n):
                if ring is basic.ZZ:
                    n = basic.ZZ.convert(n)
                    # x ** n - 1 = prod_d|n cyc(d)
                    return (cls.var() ** n.rep - 1) // (cls.product(cls.cyclotomic(d) for d in n.divisors() if d != n))
                else:
                    return cls([int(x) for x in PolyOver(basic.ZZ).cyclotomic(n)])

            #put leading term into standard associate form
            def assoc_recip(self):
                if self == 0:
                    raise ZeroDivisionError()
                else:
                    return self.rep[-1].assoc_recip()
                
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

            #euclidean stuff
            def norm(self):
                return self.degree()
            def divmod(self, other):
                assert (cls := type(self)) == type(other)
                assert other != 0
                quo = cls([])
                rem = self
                while rem.degree() - other.degree() >= 0 and rem != 0:
                    new_quo = cls([ring.int(0)] * (rem.degree() - other.degree()) + [rem[rem.degree()] / other[other.degree()]])
                    rem = rem - new_quo * other
                    quo = quo + new_quo
                return quo, rem
            def div(self, other):
                assert (cls := type(self)) == type(other)
                return self.divmod(other)[0]
            def mod(self, other):
                assert (cls := type(self)) == type(other)
                return self.divmod(other)[1]
            def truediv(self, other):
                ans = self.div(other)
                if ans * other == self:
                    return ans
                else:
                    raise basic.NotDivisibleError()
            def __divmod__(self, other):
                try:
                    other = type(self).convert(other)
                except NotImplementedError:
                    return NotImplemented
                else:
                    return self.divmod(other)

            #polynomial stuff
            def degree(self):
                #zero polynomial -> degree -1
                return len(self.rep) - 1
            def __call__(self, obj):
                #function evaluation
                if type(obj) == int:
                    return self(ring.int(obj))
                if self == 0:
                    return ring.int(0) * obj
                return type(obj).sum([c * obj ** p for p, c in enumerate(self.rep)])

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









        
