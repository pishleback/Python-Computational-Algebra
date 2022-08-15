from fractions import Fraction as Frac
from pyalgebra import combinatorics
from pyalgebra import numtheory
from pyalgebra import complex_poly_view
import functools
import itertools




class Poly():
    @staticmethod
    def string(poly):
        if len(poly) == 0:
            return "0"
        ans = ""
        for n in range(len(poly)):
            c = poly[n]
            if c != 0:
                if n == 0:
                    term = str(c)
                elif n == 1:
                    if c == 1:
                        term = "x"
                    elif c == -1:
                        term = "-x"
                    else:
                        term = str(c) + "x"
                else:
                    if c == 1:
                        term = "x^" + str(n)
                    elif c == -1:
                        term = "-x^" + str(n)
                    else:
                        term = str(c) + "x^" + str(n)
                if len(term) != 0 and term[0] != "-":
                    term = "+" + term
                ans += term    
        if len(ans) != 0:
            if ans[0] == "+":
                ans = ans[1:]
        return ans

    @staticmethod
    def validate(poly):
        assert type(poly) in [tuple, list]
        for c in poly:
            assert type(c) in [int, Frac]

    @staticmethod
    def strip(poly):
        while True:
            if len(poly) == 0:
                return poly
            elif poly[-1] == 0:
                poly = poly[:-1]
            else:
                return poly

    @staticmethod
    def add(p, q):
        if len(p) < len(q):
            p, q = q, p
        r = p[:]
        for i, c in enumerate(q):
            r[i] += c
        return Poly.strip(r)

    @staticmethod
    def mult(p, q):
        r = [0 for _ in range(len(p) + len(q))]
        for i, c in enumerate(p):
            for j, d in enumerate(q):
                r[i + j] += c * d
        return Poly.strip(r)

    @staticmethod
    def scale(p, x):
        return Poly.strip([x * c for c in p])

    @staticmethod
    def evaluate(poly, val):
        y = 0
        for coeff in reversed(poly):
            y = y * val + coeff
        return y

    @staticmethod
    def compose(p, q):
        #find p(q(x))
        r = []
        for coeff in reversed(p):
            r = Poly.add(Poly.mult(r, q), [coeff])
        return Poly.strip(r)

    @staticmethod
    def derivative(p):
        return Poly.strip([(i + 1) * p[i + 1] for i in range(len(p) - 1)])

    @staticmethod
    def primitive(p):
        p = Poly.strip(p)
        if len(p) == 0:
            return p
        p = [Frac(x) for x in p]
        m = numtheory.lcm_list([x.denominator for x in p])
        p = [(x * m).numerator for x in p]
        g = numtheory.gcd_list(p)
        p = [x // g for x in p]
        p = Poly.strip(p)
        if p[-1] < 0:
            p = [-x for x in p]
        return p

    @staticmethod
    def divmod(a, b):
        a = Poly.strip(a)
        b = Poly.strip(b)
        q = [0 for _ in range(len(a) - len(b) + 1)]
        r = list(a)
        
        for k in range(len(a) - len(b) + 1):
            m = Frac(r[-k-1], b[-1])
            q[-k-1] = m
            for i in range(len(b)):
                r[len(a) - len(b) - k + i] -= b[i] * m
##            assert tuple(Poly.strip(a)) == tuple(Poly.strip(Poly.add(Poly.mult(q, b), r)))
        return Poly.strip(q), Poly.strip(r)

    @staticmethod
    @lambda gcd : lambda p, q : gcd(tuple(Poly.strip(p)), tuple(Poly.strip(q)))
    @functools.lru_cache()
    def gcd(p, q):
        assert len(p) != 0
        assert len(q) != 0
        while len(q) != 0:
            _, r = Poly.divmod(p, q)
            p, q = q, r
        return Poly.strip(p)

    #primitive gcd, much faster than gcd when we only need the gcd up to multiplication
    @staticmethod
    @lambda gcd : lambda p, q : gcd(tuple(Poly.primitive(p)), tuple(Poly.primitive(q)))
    @functools.lru_cache()
    def pgcd(p, q):
        assert len(p) != 0
        assert len(q) != 0
        while len(q) != 0:
            _, r = Poly.divmod(Poly.scale(p, q[-1] ** (len(p) - len(q) + 1)), q)
            r = Poly.primitive(r)
            p, q = q, r
        return Poly.strip(p)

    @staticmethod
    def sqfree(p):
        p = Poly.primitive(p)
        if len(p) == 1:
            return [1]
        assert len(p) != 0
        p_prime = Poly.derivative(p)
        sqpart = Poly.pgcd(p, p_prime)
        q, r = Poly.divmod(p, sqpart)
        assert len(Poly.strip(r)) == 0
        return Poly.primitive(q)

    @staticmethod
    @lambda f : lambda p, g : f(tuple(p), g)
    @functools.lru_cache()
    def at_fixed_re(poly, gamma):
        #find real and imag polys of poly(gamma + x * i)
        deg = len(poly)
        re, im = [0 for _ in range(deg)], [0 for _ in range(deg)]
        for n in range(deg):
            if poly[n] != 0:
                for k in range(0, n + 1, 2):
                    sign = 1 if k & 2 == 0 else -1
                    re[k] += sign * poly[n] * combinatorics.binomial(n, k) * gamma ** (n - k)
                for k in range(1, n + 1, 2):
                    sign = 1 if k & 2 == 0 else -1
                    im[k] += sign * poly[n] * combinatorics.binomial(n, k) * gamma ** (n - k)
        return Poly.strip(re), Poly.strip(im)

    @staticmethod
    @lambda f : lambda p, g : f(tuple(p), g)
    @functools.lru_cache()
    def at_fixed_im(poly, gamma):
        #find real and imag polys of poly(x + gamma * i)
        deg = len(poly)
        re, im = [0 for _ in range(deg)], [0 for _ in range(deg)]
        for n in range(deg):
            if poly[n] != 0:
                for k in range(0, n + 1, 2):
                    sign = 1 if k & 2 == 0 else -1
                    re[n - k] += sign * poly[n] * combinatorics.binomial(n, k) * gamma ** k
                for k in range(1, n + 1, 2):
                    sign = 1 if k & 2 == 0 else -1
                    im[n - k] += sign * poly[n] * combinatorics.binomial(n, k) * gamma ** k
        return Poly.strip(re), Poly.strip(im)

    @staticmethod
    def sign_variations(poly):
        #https://en.wikipedia.org/wiki/Descartes'_rule_of_signs
        #equals the number of strictly positive real roots modulo 2
        #number of positive real roots is less than this number
        nzc = [x for x in poly if x != 0]
        n = 0
        for i in range(len(nzc) - 1):
            if (nzc[i] < 0) != (nzc[i + 1] < 0):
                n += 1
        return n

    #return (a, b, p) where:
    #p is a polynomial with exactly one squarefree root between a and b
    #poly should be squarefree for this to work correctly
    #finds all real roots in the closed interval [A, B]
    @staticmethod
    def isolate_real_roots(poly, A, B):
        if type(A) == int:
            A = Frac(A, 1)
        if type(B) == int:
            B = Frac(B, 1)
            
        assert len(poly) != 0
        #constant polynomial has no roots
        if len(poly) == 1:
            return ()

        #compute a bound M on the absolute value of any root
        if A is None or B is None:
            M = 2 + Frac(max(abs(poly[i]) for i in range(len(poly) - 1)), poly[-1]) #(Cauchy's bound + 1) https://captainblack.wordpress.com/2009/03/08/cauchys-upper-bound-for-the-roots-of-a-polynomial/
            if A is None:
                A = -M
            if B is None:
                B = M
        else:
            assert A < B

        #there are not roots A < r < B if A >= B
        if A >= B:
            return ()

        #now run Collins and Akritas algorithm %https://en.wikipedia.org/wiki/Real-root_isolation
        def interval(c, k, h):
            d = 2 ** k
            return ((B - A) * Frac(c, d) + A, (B - A) * Frac(c + h, d) + A)

        intervals = []
        p = Poly.compose(poly, [A, B - A])
        if Poly.evaluate(p, 1) == 0:
            intervals.append((B, B))

        #find half open half closed isolating intervals using Collins and Akritas
        L = [(0, 0, p)]
        while len(L) != 0:
            c, k, q = L.pop()
            if q[0] == 0:
                q = q[1:]
                intervals.append(interval(c, k, 0)) #rational root found
            #q'(x) := (x + 1)^nq(1 / (x + 1))
            v = Poly.sign_variations(Poly.compose(list(reversed(q)), [1, 1]))
            if v == 1:                                    
                intervals.append(interval(c, k, 1)) #algebraic root found
            elif v >= 2:
                n = len(q)
                q_small = [c * 2 ** (n - i) for i, c in enumerate(q)]
                L.append((2 * c, k + 1, q_small))
                L.append((2 * c + 1, k + 1, Poly.compose(q_small, [1, 1])))

        #shrink the intervals so that they are fully open
        open_intervals = []
        for a, b in intervals:
            if a == b: #rational root
                open_intervals.append((a, b))
            else: #algebraic root - shrink interval so that neither end point is a rational root
                if Poly.evaluate(poly, a) == 0 or Poly.evaluate(poly, b) == 0:
                    orig_a = a
                    orig_b = b
                    m = Frac(a + b, 2)
                    a = m
                    b = m
                    while True:
                        val_a = Poly.evaluate(poly, a)
                        val_b = Poly.evaluate(poly, b)
                        if val_a != 0 and val_b != 0:
                            if (Poly.evaluate(poly, a) < 0) != (Poly.evaluate(poly, b) < 0):
                                break
                        a = Frac(a + orig_a, 2)
                        b = Frac(b + orig_b, 2)
                open_intervals.append((a, b))
        return tuple(open_intervals)    

    @staticmethod
    def count_complex_roots(poly, a, b, c, d):
        assert a < b
        assert c < d
        poly = Poly.primitive(Poly.sqfree(poly))

        a_vert_re, a_vert_im = Poly.at_fixed_re(poly, a)
        b_vert_re, b_vert_im = Poly.at_fixed_re(poly, b)
        c_horz_re, c_horz_im = Poly.at_fixed_im(poly, c)
        d_horz_re, d_horz_im = Poly.at_fixed_im(poly, d)

        assert (t1r := Poly.evaluate(a_vert_re, c)) == Poly.evaluate(c_horz_re, a)
        assert (t2r := Poly.evaluate(a_vert_re, d)) == Poly.evaluate(d_horz_re, a)
        assert (t3r := Poly.evaluate(b_vert_re, c)) == Poly.evaluate(c_horz_re, b)
        assert (t4r := Poly.evaluate(b_vert_re, d)) == Poly.evaluate(d_horz_re, b)

        assert (t1i := Poly.evaluate(a_vert_im, c)) == Poly.evaluate(c_horz_im, a)
        assert (t2i := Poly.evaluate(a_vert_im, d)) == Poly.evaluate(d_horz_im, a)
        assert (t3i := Poly.evaluate(b_vert_im, c)) == Poly.evaluate(c_horz_im, b)
        assert (t4i := Poly.evaluate(b_vert_im, d)) == Poly.evaluate(d_horz_im, b)
        
        for tr, ti in [[t1r, t1i], [t2r, t2i], [t3r, t3i], [t4r, t4i]]:
            if tr == ti == 0:
                raise BoundaryRoot("Vertex root")
        

        def crossings(re, im, s, t):
            assert len(re) != 0 or len(im) != 0 #otherwise theres a line of infinite zeros which is not possible
            if len(re) == 0:
                roots_im = real_roots(im, s, t)
                if len(roots_im) == 0:
                    if Poly.evaluate(im, Frac(s + t, 2)) > 0:
                        return [1]
                    else:
                        return [3]
                else:
                    raise BoundaryRoot("Edge Root (const)")
            elif len(im) == 0:
                roots_re = real_roots(re, s, t)
                if len(roots_re) == 0:
                    if Poly.evaluate(re, Frac(s + t, 2)) > 0:
                        return [0]
                    else:
                        return [2]
                else:
                    raise BoundaryRoot("Edge Root (const)")
            else:
                roots_re = real_roots(re, s, t)
                roots_im = real_roots(im, s, t)
                for x, y in itertools.product(roots_re, roots_im):
                    if x == y:
                        raise BoundaryRoot("Edge Root")
                    RealRoot.separate(x, y)

                def sign_at(poly, root):
                    at_a = Poly.evaluate(poly, root.a)
                    at_b = Poly.evaluate(poly, root.b)
                    assert at_a != 0 and at_b != 0
                    sign_a = at_a > 0
                    sign_b = at_b > 0
                    if sign_a and sign_b:
                        return True
                    elif not sign_a and not sign_b:
                        return False
                    else:
                        #because the real and imaginary roots have been seperated
                        #for example, the interval around each real root contains no root of the imaginary part
                        #thus the sign of the imaginary part evaluated at the location of the real root is constant
                        #thus the sign is equal to the sign at either end point
                        assert False

                roots = [(x, sign_at(im, x), False) for x in roots_re] + [(x, sign_at(re, x), True) for x in roots_im]

                def to_mod4(sign, reim):
                    return {(True, True) : 0, (True, False) : 1, (False, True) : 2, (False, False) : 3}[(sign, reim)]

                return [to_mod4(tri[1], tri[2]) for tri in sorted(roots, key = lambda tri : tri[0])]

        winding = crossings(a_vert_re, a_vert_im, c, d) + crossings(d_horz_re, d_horz_im, a, b) + list(reversed(crossings(b_vert_re, b_vert_im, c, d))) + list(reversed(crossings(c_horz_re, c_horz_im, a, b)))

        #winding key
        # 0 : crossing positive real axis
        # 1 : crossing positive imag axis
        # 2 : crossing negative real axis
        # 3 : crossing negative imag axis

        #remove double entries and have the start and end be equal

        if len(winding) == 0:
            return 0
        
        norep_winding_loop = [winding[-1]]
        for n in winding:
            if n != norep_winding_loop[-1]:
                norep_winding_loop.append(n)

        quad = 0
        for i in range(len(norep_winding_loop) - 1):
            a = norep_winding_loop[i]
            b = norep_winding_loop[i + 1]
            assert a != b
            if (a - b) % 2 == 0:
                raise Exception("Something went VERY wrong")
            if (a - b) % 4 == 1:
                quad += 1
            else:
                quad -= 1

        turns, rem = divmod(quad, 4)
        assert rem == 0
        assert turns >= 0
        return turns


    #yield islating boxes for the n roots contained in the box (a, b, c, d)
    @staticmethod
    def isolate_complex_roots_boxed(poly, n, a, b, c, d):
        #print(n, a, b, c, d)

        def gen_subs(a, b, c, d):
            ab = Frac(a + b, 2)
            cd = Frac(c + d, 2)
            return [(a, ab, c, cd), (a, ab, cd, d), (ab, b, c, cd), (ab, b, cd, d)]

        def gen_subs(a, b, c, d):
            if b - a > d - c:
                m = Frac(a + b, 2)
                return [(a, m, c, d), (m, b, c, d)]
            else:
                m = Frac(c + d, 2)
                return [(a, b, c, m), (a, b, m, d)]

        subs = gen_subs(a, b, c, d)
        counts = [Poly.count_complex_roots(poly, *box) for box in subs]
        assert sum(counts) == n
        for idx in range(len(subs)):
            sa, sb, sc, sd = subs[idx]
            count = counts[idx]
            if count == 1:
                yield (sa, sb, sc, sd)
            elif count >= 2:
                c = 0
                for sub in Poly.isolate_complex_roots_boxed(poly, count, sa, sb, sc, sd):
                    yield sub
                    c += 1
                assert c == count

    @staticmethod
    def isolate_complex_roots(poly):
        poly = Poly.primitive(Poly.sqfree(poly))

        #find a bound M with all roots in the open box of side lengths 2M
        M = Frac(1, 1)
        while True:
            try:
                n = Poly.count_complex_roots(poly, -M, M, -M, M)
            except BoundaryRoot:
                pass
            else:
                if n == len(poly) - 1:
                    break
            M *= 2

        off = Frac(2, 1)
        M += 1
        while True:
            try:
                boxes = tuple(Poly.isolate_complex_roots_boxed(poly, n, -M + off - 1, M + off - 1, -M + off - 1, M + off - 1))
            except BoundaryRoot as e:
                off = Frac(off, 3)
            else:
                break
        assert len(boxes) == n
        return boxes

    #isolate the roots in the upper half plane
    @staticmethod
    def isolate_uhp_roots(poly):
        poly = Poly.primitive(Poly.sqfree(poly))
        num_real = len(Poly.isolate_real_roots(poly, None, None))
        num_uhp, r = divmod(len(poly) - 1 - num_real, 2)
        assert r == 0
        
        #find a bound M with all roots in the open box in the upper half place which tends to the whole plane as m -> infty
        M = Frac(2, 1)
        while True:
            try:
                n = Poly.count_complex_roots(poly, -M, M, Frac(1, 2 ** M), M)
            except BoundaryRoot:
                pass
            else:
                if n == num_uhp:
                    break
            M *= 2

        off = Frac(2, 1)
        M += 1
        while True:
            try:
                nudge = off - 1
                boxes = tuple(Poly.isolate_complex_roots_boxed(poly, n, -M + nudge, M + nudge, Frac(1, 2 ** M), M + 2 - off))
            except BoundaryRoot as e:
                off = Frac(off, 3)
            else:
                break
        assert len(boxes) == num_uhp
        return boxes
        
                
        
            

#represents both rational and real algebraic numbers
@functools.total_ordering
class RealRoot():
    @staticmethod
    def overlap(x, y):
        return x.b >= y.a and y.b >= x.a
    @classmethod
    def separate(cls, x, y):
        while cls.overlap(x, y):
            x.refine()
            y.refine()
            
    def __init__(self, a, b, p):
        p = Poly.primitive(Poly.sqfree(p))
        #we know that theres exaclty one root in the open interval (a, b)
        #we want to exclude the possibility of a root at an endpoint
        at_a = Poly.evaluate(p, a)
        at_b = Poly.evaluate(p, b)
        if a == b:
            assert at_a == at_b == 0
        else:
            assert at_a != 0
            assert at_b != 0
            assert (at_a < 0) != (at_b < 0)
                    
        if at_a == 0:
            b = a
            at_b = 0
        elif at_b == 0:
            a = b
            at_a = 0             
        
        self.at_a = at_a
        self.a = a
        self.at_b = at_b
        self.b = b
        self.p = p
        
    def __str__(self):
        return str(float(self))
    def __repr__(self):
        return f"Real({str(self)})"

    def __eq__(self, other):
        if (cls := type(self)) == type(other):            
            common = Poly.pgcd(self.p, other.p)
            assert len(common) != 0
            if len(common) == 1:
                return False #no common factor -> not the same root
            else:
                p = Poly.sqfree(common)
                
                if (Poly.evaluate(p, self.a) < 0) == (Poly.evaluate(p, self.b) <= 0):
                    return False
                if (Poly.evaluate(p, other.a) < 0) == (Poly.evaluate(p, other.b) <= 0):
                    return False

                roots = real_roots(p)

                for x, y in itertools.combinations(roots, 2):
                    RealRoot.separate(x, y)
                    
                while True:
                    self_overlap = [root for root in roots if cls.overlap(self, root)]
                    assert len(self_overlap) != 0
                    if len(self_overlap) == 1:
                        break
                    self.refine()
                    
                while True:
                    other_overlap = [root for root in roots if cls.overlap(other, root)]
                    assert len(other_overlap) != 0
                    if len(other_overlap) == 1:
                        break
                    other.refine()

                self_overlap = self_overlap[0]
                other_overlap = other_overlap[0]
                return self_overlap is other_overlap
            
        return False

    #need to check if one real root is less than another in order to get a handle on the order in which axis crossings happen when computing winding numbers
    def __lt__(self, other):
        if (cls := type(self)) == type(other):
            if (ans := (self.a < other.b)) == (self.b < other.a):
                return ans
            else:
                cls.separate(self, other)
                return self < other
        return NotImplemented

    def refine(self):
        if self.a == self.b:
            pass
        m = Frac(self.a + self.b, 2)            
        p_at_m = Poly.evaluate(self.p, m)
        if p_at_m == 0:
            self.at_a = p_at_m
            self.a = m
            self.at_b = p_at_m
            self.b = m
            return
        if (self.at_a < 0) != (p_at_m < 0):
            self.at_b = p_at_m
            self.b = m
        else:
            assert (p_at_m < 0) != (self.at_b < 0)
            self.at_a = p_at_m
            self.a = m

    def __float__(self):
        while self.b - self.a > 10 ** -20:
            self.refine()
        return (float(self.a) + float(self.b)) / 2

class BoundaryRoot(Exception):
    pass


def psqft(p):
    #primitive squarefree tuple version of p
    return tuple(Poly.primitive(Poly.sqfree(p)))

@lambda f : lambda p, A = None, B = None: f(psqft(p), A, B)
@functools.lru_cache()
def isolate_real_roots(p, A, B):
    Poly.validate(p)
    p = Poly.primitive(Poly.sqfree(p))
    if not A is None:
        assert type(A) in {int, Frac}
    if not B is None:
        assert type(B) in {int, Frac}
    if not A is None and not B is None:
        assert A < B
    return Poly.isolate_real_roots(p, A, B)

#count real roots in the closed interval [A, B]
def count_real_roots(p, A = None, B = None):
    Poly.validate(p)
    return len(isolate_real_roots(p, A, B))

#isolated real roots. Allows for refinement and ordering
@lambda f : lambda p, A = None, B = None : f(tuple(Poly.primitive(Poly.sqfree(p))), A, B)
@functools.lru_cache()
def real_roots(poly, A, B):
    return [RealRoot(a, b, poly) for a, b in Poly.isolate_real_roots(poly, A, B)]

def rational_real_root(x):
    assert type(x) == Frac
    return RealRoot(x, x, [-x, 1])
    
#count complex roots strictly within the box
#if a root lies on the boundary then BoundaryRoot is raised
@lambda f : lambda p, a, b, c, d : f(psqft(p), a, b, c, d)
@functools.lru_cache()
def count_complex_roots(p, a, b, c, d):
    Poly.validate(p)
    for x in [a, b, c, d]:
        assert type(x) in {int, Frac}
    assert a < b
    assert c < d
    return Poly.count_complex_roots(p, a, b, c, d)

@lambda f : lambda p: f(psqft(p))
@functools.lru_cache()
def isolate_complex_roots(p):
    Poly.validate(p)
    return Poly.isolate_complex_roots(p)

@lambda f : lambda p: f(psqft(p))
@functools.lru_cache()
def isolate_uhp_roots(p):
    Poly.validate(p)
    return Poly.isolate_uhp_roots(p)

def isolate_lhp_roots(p):
    return tuple([(a, b, -d, -c) for a, b, c, d in Poly.isolate_uhp_roots(p)])

def isolate_imag_roots(p):
    return tuple(list(isolate_uhp_roots(p)) + list(isolate_lhp_roots(p)))



def test():
    def symbol_poly(p):
        import sympy
        x = sympy.symbols("x")
        p = p(x)
        p = sympy.Poly(p, gens = x, domain = "QQ")
        p = [Frac(x.p, x.q) for x in reversed(list(p.all_coeffs()))]
        #p = [Frac(int(x.numerator), int(x.denominator)) for x in reversed(list(p.all_coeffs()))]
        return p
    p = symbol_poly(lambda x : x ** 5 - x - 1)
##    p = symbol_poly(lambda x : 61 - 273 * x + 305 * x ** 2 + x ** 3) #reeeeeeeally close roots on this one
    
    
##    p = symbol_poly(lambda x : x ** 12 - 3 * x ** 5 + 6 * x ** 2 - x + 1)
##    p = symbol_poly(lambda x : x ** 12 - x ** 11 - 2)
##    import random
##    p = symbol_poly(lambda x : sum(random.randint(0, 1) * x ** k for k in range(10)))
##    p = symbol_poly(lambda x : (x + x ** 2) ** 10 - 1)
##    p = symbol_poly(lambda x : x ** 2 - 1)    

    print(p, len(p))
    
    print("REAL")
    for a, b in isolate_real_roots(p):
        print(a, b)

    print("ALL")
    for a, b, c, d in isolate_complex_roots(p):
        print(a, b, c, d)

    print("UHP")
    for a, b, c, d in isolate_uhp_roots(p):
        print(a, b, c, d)

    #boxes for imaginary roots and intervals for real roots
    def gen_boxes(p):
        yield from isolate_imag_roots(p)
        for a, b in isolate_real_roots(p):
            if a == b:
                yield (a-0.01, b+0.01, -0.02, 0.02)
            else:
                yield (a, b, -0.01, 0.01)

    #boxes for all roots
    def gen_boxes(p):
        yield from isolate_complex_roots(p)
            
    complex_poly_view.run(p, list(gen_boxes(p)))

    




































