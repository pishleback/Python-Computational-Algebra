import algebra
import pyalgebra


ZZ = algebra.base.ZZ
QQ = algebra.base.ZZ.FractionField
PolyZZ = algebra.polynomials.PolyOver(ZZ)
PolyQQ = algebra.polynomials.PolyOver(QQ)
QQ_bar = algebra.algebraic.Algebraic



def test2():
##    pyalgebra.polynomials.test()
##
##    input()
    
    from fractions import Fraction as Frac

    x = PolyQQ.var()
    poly = x ** 16  - 1

    for f in poly.factor().list():
        print(f)

    for a, k in QQ.AlgebraicClosure.root_powers(poly):
        print(k, a, a.degree())

##    x = PolyZZ.var()
##
####    a = algebra.algebraic._RealRep(x ** 2 - 2 * x - 1, Frac(-1, 1), Frac(0, 1))
####    b = algebra.algebraic._RealRep(x ** 2 - 2 * x - 1, Frac(2, 1), Frac(3, 1))
##
##    t = 0
##    poly = x ** 24 - 1
##    roots = list(algebra.algebraic.Algebraic.poly_roots(poly))
##    for a in roots:
##        print(a, a + 1, (a + 1).min_poly())
    
##    a = algebra.algebraic.RealRep(x ** 5 - x - 1, Frac(1, 1), Frac(2, 1))
##    b = algebra.algebraic.RealRep(x ** 3 - x - 1, Frac(1, 1), Frac(2, 1))
##    c = algebra.algebraic.RealRep(2 * x ** 4 - 6 * x ** 2 + x + 1, Frac(-2, 1), Frac(-1, 1))
##    d = algebra.algebraic.RealRep(2 * x ** 4 - 6 * x ** 2 + x + 1, Frac(-1, 1), Frac(0, 1))
##    e = algebra.algebraic.RealRep(2 * x ** 4 - 6 * x ** 2 + x + 1, Frac(0, 1), Frac(1, 1))
##    f = algebra.algebraic.RealRep(2 * x ** 4 - 6 * x ** 2 + x + 1, Frac(1, 1), Frac(2, 1))
##
##    h = Frac(10, 1)
##    i = Frac(-9, 5)
##
##    print(a, b, a * b)
##    print(c, d, c * d)
##    print(e, f, e * f)
##    print(a, h, a * h)
##
##    print(a, a.recip())
##    print(a, b.recip())

    

def test1():    
    ZZ = algebra.base.ZZ
    QQ = algebra.base.ZZ.FractionField
    G = algebra.base.Gaussian
    i = G(0, 1)

    QQbar = QQ.AlgebraicClosure

    print(QQ.AlgebraicClosure is QQ.AlgebraicClosure)

    M = algebra.matricies.MatrixOver(QQ)

    A = M(4, 4, [[0, -1, 1, 1],
                 [1, 0, 1, 1],
                 [0, 0, 0, -1],
                 [0, 0, 1, 0]])



##    A = M(3, 3, [[1, 0, 0],
##                 [1, 0, 0],
##                 [1, 0, 0]])
##
##    print(A.col_span().compliment())

    

    print(A)
    for x in A.eigen_val_list():
        print(x, x.degree(), x.min_poly())
    print(A.jordan_canonical_form())


def test3():
    ZZ = algebra.base.ZZ
    QQ = algebra.base.ZZ.FractionField
    QQbar = QQ.AlgebraicClosure

    M = algebra.matricies.MatrixOver(QQ)
    SP = algebra.matricies.AffineSubspaceOver(ZZ)
    
    A = M(4, 4, [[2, 0, 0, 0],
                 [0, 2, 0, 0],
                 [0, 0, 2, 0],
                 [0, 0, 0, 2]])
##    V = M(4, 1, [[1],
##                 [1],
##                 [1],
##                 [1]])
##    S = SP(4, 1, V, A.col_list())
##    print(S)
##    
    B = M(4, 4, [[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 1, 0],
                 [0, 0, 0, 1]])
##    W = M(4, 1, [[1],
##                 [1],
##                 [1],
##                 [1]])
##    T = SP(4, 1, W, B.col_list()) 
##    print(T)
    
    print(A.col_span())
    print(B.col_span())
    print(A.col_span() == B.col_span())

    
    

test3()
