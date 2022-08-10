import functools
from algebra import base
import itertools

#===SUMMARY===

#==DATA STRUCTURE==
#put numbers into a matrix

#==OVER AN ED==
#=SPAN= (ModuleSpan)
#objects are submodules
#=QUOTIENT= (ModuleQuotient)
#objects are elements mod submodules
#isomorphism type from smith normal form (SNF)

#==OVER A FIELD== (SPECIAL CASE OF ED)
#=SPAN= (VectorSpan)
#objects are subspaces
#=quotient= (VectorQuotient)
#QUOTIENT are vectors mod subspace
#isomorphism type from dimention = how many diagonal entries in the SNF

#===APPLICATIONS===
#homology
#integer systems of equations
#linear systems of equations
#jordan normal form
#character tables
#algebraic numbers

def assoc_opp(f):
    def new_f(cls, values):
        ans = values[0]
        for value in values[1:]:
            ans = f(cls, ans, value)
        return ans
    return new_f


class SizeMismatch(Exception):
    pass

class NoSolution(Exception):
    pass


@functools.cache
def MatrixOver(ring):
    assert issubclass(ring, base.Ring)
    class Matrix(base.MathsSet):
        @classmethod
        def test_axioms(cls, test):
            super().test_axioms(test)
                
        @classmethod
        def typestr(cls):
            return f"Mat({ring})"

        @classmethod
        def sum(cls, rows, cols, mats):
            ans = cls.zero(rows, cols)
            for mat in mats:
                ans += mat
            return ans

        @classmethod
        @assoc_opp
        def join_rows(cls, mat1, mat2):
            assert (cols := mat1.cols) == mat2.cols
            return cls(mat1.rows + mat2.rows, cols, mat1.entries + mat2.entries)
        @classmethod
        @assoc_opp
        def join_cols(cls, mat1, mat2):
            assert (rows := mat1.rows) == mat2.rows
            return cls(rows, mat1.cols + mat2.cols, [mat1.entries[r] + mat2.entries[r] for r in range(rows)])
        @classmethod
        @assoc_opp
        def join_diag(cls, mat1, mat2):
            assert type(mat1) == type(mat2) == cls
            bl_zeros = cls(mat2.rows, mat1.cols, [[0 for c in range(mat1.cols)] for r in range(mat2.rows)])
            tr_zeros = cls(mat1.rows, mat2.cols, [[0 for c in range(mat2.cols)] for r in range(mat1.rows)])
            return cls.join_cols([cls.join_rows([mat1, bl_zeros]), cls.join_rows([tr_zeros, mat2])])

        @classmethod
        def eye(cls, n):
            return cls(n, n, [[(1 if r == c else 0) for c in range(n)] for r in range(n)])
        @classmethod
        def zero(cls, r, c = None):
            if c is None:
                c = r
            return cls(r, c, [[0 for _ in range(c)] for _ in range(r)])

        def __init__(self, rows, cols, entries):
            entries = tuple(tuple(ring.convert(x) for x in row) for row in entries)
            assert len(entries) == rows
            for row in entries:
                assert len(row) == cols

            self.rows = rows
            self.cols = cols
            self.entries = entries

        def __getitem__(self, pos):
            assert type(pos) == tuple and len(pos) == 2
            r, c = pos
            return self.entries[r][c]

        def __str__(self):
            strs = [[str(self[r, c]) for c in range(self.cols)] for r in range(self.rows)]
            lens = [max(len(strs[r][c]) for r in range(self.rows)) for c in range(self.cols)]
            strs = [[strs[r][c] + " " * (lens[c] - len(strs[r][c])) for c in range(self.cols)] for r in range(self.rows)]

            def brac(r):
                if self.rows == 1:
                    return "["
                else:
                    if r == 0:
                        return "/"
                    elif r == self.rows - 1:
                        return "\\"
                    else:
                        return "|"
            def flip(b):
                return {"[" : "]", "/" : "\\", "\\" : "/", "|" : "|"}[b]
                
            return "\n".join(brac(r) + " ".join(strs[r][c] for c in range(self.cols)) + flip(brac(r)) for r in range(self.rows))
        
        def __repr__(self):
            return f"Mat({ring}, [{', '.join('[' + ', '.join(repr(x) for x in row) + ']' for row in self.entries)}])"

        def equal(self, other):
            assert (cls := type(self)) == type(other)
            if (rows := self.rows) == other.rows:
                if (cols := self.cols) == other.cols:
                    for r in range(rows):
                        for c in range(cols):
                            if self[r, c] != other[r, c]:
                                return False
                    return True
            return False
        def add(self, other):
            assert (cls := type(self)) == type(other)
            if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                return cls(rows, cols, [[self[r, c] + other[r, c] for c in range(cols)] for r in range(rows)])
            else:
                raise SizeMismatch("Matricies must have the same dimentions for addition")
        def neg(self):
            return type(self)(self.rows, self.cols, [[-self[r, c] for c in range(self.cols)] for r in range(self.rows)])
        def mul(self, other):
            assert (cls := type(self)) == type(other)
            if (mids := self.cols) == other.rows:
                rows = self.rows
                cols = other.cols
                return cls(rows, cols, [[ring.sum([self[r, m] * other[m, c] for m in range(mids)]) for c in range(cols)] for r in range(rows)])
            else:
                raise SizeMismatch("Matricies dont have compatible dimention for matrix multiplication")
            
        def scalarmul_right(self, other):
            assert isinstance(other, ring)
            return type(self)(self.rows, self.cols, [[self[r, c] * other for c in range(self.cols)] for r in range(self.rows)])
        def scalarmul_left(self, other):
            assert isinstance(other, ring)
            return type(self)(self.rows, self.cols, [[other * self[r, c] for c in range(self.cols)] for r in range(self.rows)])

        def __hash__(self):
            return hash(self.entries)
        def __eq__(self, other):
            if type(self) == type(other):
                return self.equal(other)
            return False
        def __add__(self, other):
            if type(self) == type(other):
                return self.add(other)
            return NotImplemented
        def __neg__(self):
            return self.neg()
        def __sub__(self, other):
            if type(self) == type(other):
                return self.add(other.neg())
            return NotImplemented
        def __mul__(self, other):
            if type(self) == type(other):
                return self.mul(other)
            try:
                other = ring.convert(other)
            except NotImplementedError:
                pass
            else:
                return self.scalarmul_right(other)
            return NotImplemented
        def __rmul__(self, other):
            if type(self) == type(other):
                return other.mul(self)
            try:
                other = ring.convert(other)
            except NotImplementedError:
                pass
            else:
                return self.scalarmul_left(other)
            return NotImplemented
        def __pow__(self, other):
            cls = type(self)
            if (n := self.rows) == self.cols:
                if type(other) == int:
                    if other < 0:
                        return (self ** (-other)).inverse()
                    elif other == 0:
                        return cls.eye(n)
                    elif other == 1:
                        return self
                    elif other == 2:
                        return self * self
                    else:
                        q, r = divmod(other, 2)
                        return (self ** q) ** 2 * self ** r
            return NotImplemented

        def transpose(self):
            return type(self)(self.cols, self.rows, [[self[c, r] for c in range(self.rows)] for r in range(self.cols)])
        def flip_rows(self):
            return type(self)(self.rows, self.cols, [[self[self.rows - r - 1, c] for c in range(self.cols)] for r in range(self.rows)])
        def flip_cols(self):
            return type(self)(self.rows, self.cols, [[self[r, self.cols - c - 1] for c in range(self.cols)] for r in range(self.rows)])
        def row(self, r):
            return type(self)(1, self.cols, [self.entries[r]])
        def col(self, c):
            return self.transpose().row(c).transpose()
        def minor(self, row = None, col = None):
            return type(self)(self.rows if row is None else self.rows - 1, self.cols if col is None else self.cols - 1, [[self[r, c] for c in range(self.cols) if c != col] for r in range(self.rows) if r != row])

        def row_swap(self, r1, r2): #r1 <-> r2
            assert type(r1) == int and 0 <= r1 < self.rows
            assert type(r2) == int and 0 <= r2 < self.rows
            assert r1 != r2
            def perm(r):
                if r == r1:
                    return r2
                elif r == r2:
                    return r1
                else:
                    return r
            return type(self)(self.rows, self.cols, [[self[perm(r), c] for c in range(self.cols)] for r in range(self.rows)])
        def row_mult(self, r1, m): #r1 -> m * r1
            assert type(r1) == int and 0 <= r1 < self.rows
            m = ring.convert(m)
            return type(self)(self.rows, self.cols, [[(m * self[r1, c] if r == r1 else self[r, c]) for c in range(self.cols)] for r in range(self.rows)])
        def row_add(self, r1, r2, m): #r1 -> r1 + m * r2
            assert r1 != r2
            assert type(r1) == int and 0 <= r1 < self.rows
            assert type(r2) == int and 0 <= r2 < self.rows
            m = ring.convert(m)
            return type(self)(self.rows, self.cols, [[(self[r1, c] + m * self[r2, c] if r == r1 else self[r, c]) for c in range(self.cols)] for r in range(self.rows)])

        def col_swap(self, c1, c2): #c1 <-> c2
            return self.transpose().row_swap(c1, c2).transpose()
        def col_mult(self, c1, m): #c1 -> m * c1
            return self.transpose().row_mult(c1, m).transpose()
        def col_add(self, c1, c2, m): #c1 -> c1 + m * c2
            return self.transpose().row_add(c1, c2, m).transpose()

        @functools.lru_cache()
        def hermite_algorithm(self):
            if not issubclass(ring, base.EuclideanDomain):
                raise Exception(f"Cant do hermite algorithm over non-ED {ring} (yet)")
            #put the matrix into hermite normal form / reduced row echelon form
            #this form is unique for a given matrix
            
            #return U, det(U), pivots, H such that U * self == H
            #U is invertible
            #det(U) is det(U) (obviously)
            #H is in hermite normal form
            #pivots is a tuple of the columns of the pivots of H

            #algorithm idea: cancel columns with analouge of euclids algorithm
            #when the matrix is over a field, this is guassian elimination

            n = self.rows
            U = type(self).eye(n)
            H = self
            pivots = []
            det_U = ring.int(1)
            
            r = 0
            c = 0
            while True:
                while True:
                    #if the column is all zero, skip it
                    non_zero = [br for br in range(r, n) if H[br, c] != 0]
                    if len(non_zero) == 0:
                        c += 1
                        break
                    else:
                        #find non-zero element in col c below r with minimal norm
                        min_r = min(non_zero, key = lambda br : H[br, c].norm())
                        if min_r != r:
                            H = H.row_swap(r, min_r)
                            U = U.row_swap(r, min_r)
                            det_U = -det_U
                        
                        if all(H[br, c] == 0 for br in range(r + 1, n)):
                            #if all elements below r are zero, then reduce all elements above r and continue
                            for ar in range(r):
                                if H[ar, c] != 0:
                                    m = -(H[ar, c] // H[r, c])
                                    H = H.row_add(ar, r, m)
                                    U = U.row_add(ar, r, m)
                            #make sure: if pivot elements are units, then they are 1
                            try:
                                m = H[r, c].recip()
                            except base.NotDivisibleError:
                                pass
                            else:
                                H = H.row_mult(r, m)
                                U = U.row_mult(r, m)
                                det_U *= m
                            pivots.append(c)
                            r += 1
                            c += 1
                            break

                        #subtract multiples of the minimal element to reduce all other norms in the column below r (bit like euclids algorithm)
                        for br in range(r + 1, n):
                            m = -(H[br, c] // H[r, c])
                            H = H.row_add(br, r, m)
                            U = U.row_add(br, r, m)
                            
                if r >= H.rows or c >= H.cols:
                    assert U * self == H
                    return U, det_U, tuple(pivots), H
        
        def inverse(self):
            assert (n := self.rows) == self.cols
            U, det_U, pivots, H = self.hermite_algorithm()
            if H == type(self).eye(n):
                return U
            else:
                raise Exception("Matrix is not invertible")

        @functools.lru_cache()
        def smith_algorithm(self):
            if not issubclass(ring, base.EuclideanDomain):
                raise Exception(f"Cant do smith algorithm over non-ED {ring}")
            #find invertible S, T such that S * self * T is in smith normal form

            #idea:
            #for each diagonal point
            # - put the element of smallest norm into that spot
            # - reduce others in the row and column
            # - repeat until all are zero
            
            cls = type(self)
            
            S = cls.eye(self.rows)
            A = self
            T = cls.eye(self.cols)
            k = 0 #which diagonal are we doing rn

            while True:
                #if the matrix is all zero, its already in smith normal form
                non_zero = [(r, c) for r, c in itertools.product(range(k, self.rows), range(k, self.cols)) if A[r, c] != 0]
                if len(non_zero) == 0:
                    return S, A, T
                else:
                    #seach for the non-zero element of minimal norm
                    r, c = min(non_zero, key = lambda pos : A[pos].norm())

                    #place it in the top left of A
                    if r != k:
                        A = A.row_swap(r, k)
                        S = S.row_swap(r, k)
                    if c != k:
                        A = A.col_swap(c, k)
                        T = T.col_swap(c, k)

                    #find non-zero elements in the top row and first column (excluding the top left point)
                    non_zero_row = [r for r in range(k + 1, self.rows) if A[r, k] != 0]
                    non_zero_col = [c for c in range(k + 1, self.cols) if A[k, c] != 0]

                    if len(non_zero_row) == 0 and len(non_zero_col) == 0:
                        #if there are none, then euclidian algorithm is complete

                        #put diagonal elements into standard associate form (e.g. for Z, make them positive. for fields, make them 1)
                        try:
                            m = A[k, k].recip()
                        except base.NotDivisibleError:
                            pass
                        else:
                            A = A.row_mult(k, m)
                            S = S.row_mult(k, m)

                        #now we want all other elements to be multiples of this top left element
                        for r, c in itertools.product(range(k + 1, self.rows), range(k + 1, self.cols)):
                            if A[r, c] % A[k, k] != 0:
                                #if not, then fiddle about a bit and do a single euclid to reduce the norm of this element
                                A = A.row_add(r, k, 1)
                                S = S.row_add(r, k, 1)

                                m = -(A[r, c] // A[r, k])
                                A = A.col_add(c, k, m)
                                T = T.col_add(c, k, m)

                                #now call smith algorithm on this new matrix with smaller norms appearing
                                #will terminate as new top left norm must cannot decrease forever
                                S_prime, A, T_prime = A.smith_algorithm()
                                S = S_prime * S
                                T = T * T_prime
                                return S, A, T

                        #if they are, then move on to a sub matrix or return if we are done
                        k += 1
                        if k >= self.rows or k >= self.cols:
                            return S, A, T
                    else:
                        #if there are non-zero elements, do some euclids to reduce their norm and repeat
                        for r in non_zero_row:
                            m = -(A[r, k] // A[k, k])
                            A = A.row_add(r, k, m)
                            S = S.row_add(r, k, m)
                            
                        for c in non_zero_col:
                            m = -(A[k, c] // A[k, k])
                            A = A.col_add(c, k, m)
                            T = T.col_add(c, k, m)
                        
        def smith_normal_form(self):
            return self.smith_algorithm()[1]
        def smith_diagonal(self):
            snf = self.smith_normal_form()
            return tuple(snf[i, i] for i in range(min(self.rows, self.cols)))

        def rank(self):
            return len(self.hermite_algorithm()[2])
        def trace(self):
            assert (n := self.rows) == self.cols
            return ring.sum([self[i, i] for i in range(n)])
        def det(self):
            assert (n := self.rows) == self.cols
            U, det_U, pivots, H = self.hermite_algorithm()
            #U * self == H
            #det(H) = product of diagonals
            #det_U is given
            # => det(self) = det(H) / det(U)
            return ring.product([H[i, i] for i in range(n)]) / det_U

        def row_span(self):
            return SpanOver(ring)(1, self.cols, [self.row(r) for r in range(self.rows)])
        def col_span(self):
            return SpanOver(ring)(self.rows, 1, [self.col(c) for c in range(self.cols)])
        def row_kernel(self):
            U, _, pivs, H = self.hermite_algorithm()
            return SpanOver(ring)(1, self.rows, [U.row(r) for r in range(len(pivs), self.rows)])
        def col_kernel(self):
            return self.transpose().row_kernel().transpose()

        def col_solve(self, vec):
            #solve self * x = vec for x
            assert self.cols == vec.rows and vec.cols == 1
            U, det_U, pivots, H = self.hermite_algorithm()
            
            mat = type(self).join_cols([self.col(c) for c in range(self.cols)] + [-vec])
            ker_basis = mat.col_kernel().basis
            g, coeffs = ring.xgcd_list([b[self.cols, 0] for b in ker_basis])
            if g.is_unit():
                coeffs = [c / g for c in coeffs]
            else:
                raise NoSolution("No solution")

            sol = Matrix.sum(self.cols + 1, 1, [coeffs[i] * ker_basis[i] for i in range(len(ker_basis))])
            sol = sol.minor(row = self.cols)
            assert self * sol == vec
            return sol
            
    return Matrix
    







@functools.cache
def SpanOver(ring):
    assert issubclass(ring, base.Ring)
    Matrix = MatrixOver(ring)

    def mat_to_row(rows, cols, mat):
        assert mat.rows == rows
        assert mat.cols == cols
        return Matrix(1, cols * rows, [[mat[i % rows, i // rows] for i in range(cols * rows)]])

    def row_to_mat(rows, cols, row_mat):
        assert row_mat.cols == rows * cols
        assert row_mat.rows == 1
        return Matrix(rows, cols, [[row_mat[0, r + c * rows] for c in range(cols)] for r in range(rows)])
    
    class Span():
        def __init__(self, rows, cols, matricies):
            assert type(rows) == int and rows >= 0
            assert type(cols) == int and rows >= 0
            for matrix in matricies:
                assert isinstance(matrix, Matrix)
                assert matrix.rows == rows
                assert matrix.cols == cols
                
            self.rows = rows
            self.cols = cols

            if len(matricies) == 0:
                self.basis = tuple([])
            else:
                _, _, pivs, H = Matrix.join_rows([mat_to_row(self.rows, self.cols, mat) for mat in matricies]).hermite_algorithm()
                self.basis = tuple(row_to_mat(self.rows, self.cols, H.row(i)) for i in range(len(pivs)))
        def __str__(self):
            if len(self.basis) == 0:
                return "{0}"
            else:
                rows = [""] * self.rows
                for idx, mat in enumerate(self.basis):
                    if idx != 0:
                        for r in range(self.rows):
                            if r == self.rows - 1:
                                rows[r] += " , "
                            else:
                                rows[r] += "   "
                    mat_str_rows = str(mat).split("\n")
                    for r in range(self.rows):
                        rows[r] += mat_str_rows[r]
                return "\n".join(rows)
        def __repr__(self):
            return "Span(" + ", ".join(repr(mat) for mat in self.basis) + ")"

        def __add__(self, other):
            if (cls := type(self)) == type(other):
                if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                    return cls(rows, cols, self.basis + other.basis)
            elif type(other) == Matrix:
                return self.as_offsetspan() + other
            return NotImplemented
        def __radd__(self, other):
            if type(other) == Matrix:
                return other + self.as_offsetspan()
            return NotImplemented
        def __and__(self, other):
            if (cls := type(self)) == type(other):
                if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                    metamatrix = Matrix.join_rows([mat_to_row(self.rows, self.cols, mat) for mat in self.basis] + [mat_to_row(other.rows, other.cols, mat) for mat in other.basis])
                    ker = metamatrix.row_kernel()

                    basis_rows = []
                    for coeffs in ker.basis:
                        basis_row = Matrix.zero(1, metamatrix.cols)                    
                        assert coeffs.cols == self.dimention() + other.dimention()
                        for i in range(0, self.dimention()):
                            basis_row += coeffs[0, i] * metamatrix.row(i)
                        basis_rows.append(basis_row)

                    return cls(rows, cols, [row_to_mat(self.rows, self.cols, row) for row in basis_rows])
            return NotImplemented

        def __mul__(self, other):
            if type(other) == Matrix:
                if other.rows == self.cols:
                    return type(self)(self.rows, other.cols, [mat * other for mat in self.basis])
                else:
                    raise Exception("Dimentions don't match for Span * Mat")
            return NotImplemented

        def __rmul__(self, other):
            if type(other) == Matrix:
                if other.cols == self.rows:
                    return type(self)(other.rows, self.cols, [other * mat for mat in self.basis])
                else:
                    raise Exception("Dimentions don't match for Mat * Span")
            return NotImplemented

        def dimention(self):
            return len(self.basis)
        def sample(self):
            if self.dimention == 0:
                raise Exception("No elements to sample from")
            else:
                return self.basis[0]

        def transpose(self):
            return type(self)(self.cols, self.rows, [mat.transpose() for mat in self.basis])

    return Span








if False:
    #old stuff
    @functools.cache
    def AllOver(ring):
        assert issubclass(ring, basic.Ring)
        
        class Matrix(basic.Ring):
            def char_mat(self):
                assert (n := self.rows) == self.cols
                polyring = polynomials.PolyOver(ring)
                Pmat = MatrixOver(polyring)
                pself = Pmat([[polyring.convert(self[r, c]) for c in range(self.cols)] for r in range(self.rows)])
                return pself - polyring.var() * Pmat.eye(n)
            def min_poly(self):
                return self.char_mat().smith_diagonal()[-1]
            def char_poly(self):
                return polynomials.PolyOver(ring).product(self.char_mat().smith_diagonal())

            def row_span(self):
                return SpanOver(ring)(1, self.cols, [self.row(r) for r in range(self.rows)])
            def col_span(self):
                return SpanOver(ring)(self.rows, 1, [self.col(c) for c in range(self.cols)])
            def row_kernel(self):
                U, _, pivs, H = self.hermite_algorithm()
                return SpanOver(ring)(1, self.rows, [U.row(r) for r in range(len(pivs), self.rows)])
            def col_kernel(self):
                return self.transpose().row_kernel().transpose()

            def row_sol_sp(self, a):
                return self.transpose().col_sol_sp(a.transpose()).transpose()
            def col_sol_sp(self, a):
                assert a.cols == 1
                assert (rows := self.rows) == a.rows
                sp = type(self).eye(self.cols).col_span().as_offsetspan()
                for r in range(rows):
                    row = self.row(r)
                    g, coeffs = ring.xgcd_list([row[0, i] for i in range(self.cols)])
                    if g == 0 and a[r, 0] != 0:
                        return OffsetSpanOver(ring).empty(self.cols, 1)
                    elif g == 0:
                        pass
                    elif a[r, 0] % g == 0:
                        mult = a[r, 0] // g
                        sp &= type(self)([[mult * elem] for elem in coeffs]) + row.col_kernel()
                    else:
                        return OffsetSpanOver(ring).empty(self.cols, 1)
                return sp

            def col_unique_sol(self, a):
                sol_sp = self.col_sol_sp(a)
                offset, span = sol_sp.offset_span()
                assert span.dimention() == 0
                return offset

    ##        def fractional_row_sol_sp(self, a):
    ##            raise NotImplementedError()
    ##        def fractional_col_sol_sp(self, a):
    ##            assert issubclass(ring, basic.FieldOfFractions)
    ##            assert type(self) == type(a)
    ##            assert a.rows == self.rows
    ##            assert a.cols == 1
    ##            m1, selfp = self.primitive()
    ##            m2, ap = a.primitive()
    ##            m = ring.ring.lcm(m1, m2)
    ##            selfp = selfp * (m // m1)
    ##            ap = ap * (m // m2)            
    ##            sols_p = selfp.col_sol_sp(ap)
    ##            return m, sols_p
    ##

            def primitive(self):
                if issubclass(ring, basic.FieldOfFractions):
                    elems = [self[r, c] for r, c in itertools.product(range(self.rows), range(self.cols))]
                    m = ring.ring.lcm_list([elem.d for elem in elems])
                    return m, MatrixOver(ring.ring)([[(m * self[r, c]).to_ring() for c in range(self.cols)] for r in range(self.rows)])
                else:
                    raise Exception("Can only decompose into multiple and primitive matrix over a field of fractions")

            def convert_to(self, newring):
                return MatrixOver(newring)([[self[r, c] for c in range(self.cols)] for r in range(self.rows)])








        def mat_to_row(rows, cols, mat):
            assert mat.rows == rows
            assert mat.cols == cols
            return MatrixOver(ring)([[mat[i % rows, i // rows] for i in range(cols * rows)]])

        def row_to_mat(rows, cols, row_mat):
            assert row_mat.cols == rows * cols
            assert row_mat.rows == 1
            return MatrixOver(ring)([[row_mat[0, r + c * rows] for c in range(cols)] for r in range(rows)])

        def projectivize(rows, cols, div, mat):
            row = mat_to_row(rows, cols, mat)
            hom = MatrixOver(ring).join_cols([MatrixOver(ring)([[div]]), row])
            return hom

        def deprojectivize(rows, cols, mat):
            row = mat_to_row(rows, cols, mat)
            return row.minor(col = 0), row[0, 0]




        class OffsetSpan():
            class EmptyError(Exception):
                pass
            
            @classmethod
            def empty(cls, rows, cols):
                return cls(rows, cols, SpanOver(ring)(1, rows * cols + 1, []))
                
            def __init__(self, rows, cols, proj_span):
                assert isinstance(proj_span, SpanOver(ring))
                assert proj_span.rows == 1
                assert proj_span.cols == rows * cols + 1
                self.rows = rows
                self.cols = cols
                self.n = rows * cols + 1
                
                self.proj_span = proj_span

            def __str__(self):
                if self.is_empty():
                    return "EMPTY"
                else:
                    offset, span = self.offset_span()
                    
                    if span.dimention() == 0:
                        return str(offset) + " + <>"

                    offset_lines = str(offset).split("\n")
                    span_lines = str(span).split("\n")

                    return "\n".join([offset_lines[i] + (" + < " if i == len(offset_lines) - 1 else "     ") + span_lines[i] + (" >" if i == len(offset_lines) - 1 else "  ") for i in range(len(offset_lines))])                
                
            def __add__(self, other):
                if (cls := type(self)) == type(other):
                    if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                        raise Exception("im not so sure that this should be implemented")
                        return cls(rows, cols, self.proj_span + other.proj_span)
                elif type(other) == MatrixOver(ring):
                    if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                        return cls(rows, cols, SpanOver(ring)(1, self.n, [mat + mat[0, 0] * projectivize(self.rows, self.cols, 0, other) for mat in self.proj_span.basis]))
                return NotImplemented
            def __radd__(self, other):
                if type(other) == MatrixOver(ring):
                    return self + other
                return NotImplemented
            def __sub__(self, other):
                if type(other) == MatrixOver(ring):
                    return self + (-other)
                return NotImplemented

            def __mul__(self, other):
                if type(other) == MatrixOver(ring):
                    if other.rows == self.cols and self.rows == 1:
                        hom_other = MatrixOver(ring).join_diag([MatrixOver(ring).eye(1), other])
                        return type(self)(other.rows, self.cols, type(self.proj_span)(1, other.rows * self.cols + 1, [hom * hom_other for hom in self.proj_span.basis]))
                    else:
                        raise Exception("Dimentions don't match for ProjSpan * Mat")
                return NotImplemented
            def __rmul__(self, other):
                if type(other) == MatrixOver(ring):
                    if other.cols == self.rows and self.cols == 1:
                        hom_other = MatrixOver(ring).join_diag([MatrixOver(ring).eye(1), other])
                        return type(self)(other.rows, self.cols, type(self.proj_span)(1, other.rows * self.cols + 1, [(hom_other * hom.transpose()).transpose() for hom in self.proj_span.basis]))
                    else:
                        raise Exception("Dimentions don't match for Mat * ProjSpan")
                return NotImplemented

            def __and__(self, other):
                if (cls := type(self)) == type(other):
                    if (rows := self.rows) == other.rows and (cols := self.cols) == other.cols:
                        return cls(rows, cols, self.proj_span & other.proj_span)
                return NotImplemented

            #decompose into an offset and a span
            def offset_span(self):
                m = len(self.proj_span.basis)
                if m == 0:
                    raise self.EmptyError("VERY EMPTY")
                
                metamat = MatrixOver(ring).join_rows(self.proj_span.basis)
                metamat = metamat.hermite_normal_form()
                assert metamat.cols == self.n
                assert metamat.rows == m

                if metamat[0, 0] == 0:
                    raise self.EmptyError("EMPTY")
                else:
                    g = metamat[0, 0]
                    if g == 1:
                        for i in range(1, m):
                            assert metamat[i, 0] == 0

                        offset = row_to_mat(self.rows, self.cols, metamat.row(0).minor(col = 0))
                        span = SpanOver(ring)(self.rows, self.cols, [row_to_mat(self.rows, self.cols, metamat.row(i).minor(col = 0)) for i in range(1, m)])

                        return offset, span
                    else:
                        raise self.EmptyError("KINDA EMPTY")

            def is_empty(self):
                try:
                    self.offset_span()
                except self.EmptyError:
                    return True
                else:
                    return False

            def transpose(self):
                return type(new_self_ex)(self.cols, self.rows, SpanOver(ring)(1, self.n, [MatrixOver(ring).join_cols([mat_to_row(self.cols, self.rows, row_to_mat(self.rows, self.cols, hom.minor(col = self.n - 1)).transpose()), MatrixOver(ring)([[hom[0, self.n - 1]]])]) for hom in self.proj_span.basis]))
            

        return Matrix, Span, OffsetSpan



    def MatrixOver(ring):
        return AllOver(ring)[0]
    def SpanOver(ring):
        return AllOver(ring)[1]
    def OffsetSpanOver(ring):
        return AllOver(ring)[2]







    def QQ_possible_orders(n):
        import pyalgebra
        #yield all possible (finite) orders of nxn matrixies over QQ

        #partition n
        #take all possible phi inverse of the partition
        #take lcm of the phi inverse partition
        #yield that value

        found = set([])
        for part in pyalgebra.combinatorics.int_partitions(n):
            for phi_inv in itertools.product(*[tuple(pyalgebra.numtheory.phi_inverse(k)) for k in part]):
                x = pyalgebra.numtheory.lcm_list(phi_inv)
                if not x in found:
                    found.add(x)
                    yield x

        























