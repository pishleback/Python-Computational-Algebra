import algebra
import pyalgebra




def test2():
##    pyalgebra.polynomials.test()
    
    from fractions import Fraction as Frac
    
    ZZ = algebra.base.ZZ
    QQ = algebra.base.ZZ.FractionField
    PolyZZ = algebra.polynomials.PolyOver(ZZ)
    PolyQQ = algebra.polynomials.PolyOver(QQ)

    x = PolyZZ.var()

##    a = algebra.algebraic._RealRep(x ** 2 - 2 * x - 1, Frac(-1, 1), Frac(0, 1))
##    b = algebra.algebraic._RealRep(x ** 2 - 2 * x - 1, Frac(2, 1), Frac(3, 1))

    poly = x ** 12 - 1
    for poly in poly.factor().list():
        q = [int(c) for c in poly.rep]
        #pyalgebra.complex_poly_view.run(q, pyalgebra.polynomials.isolate_complex_roots(q))
        
        for a, b, c, d in pyalgebra.polynomials.isolate_complex_roots([int(c) for c in poly.rep]):
            print(poly)
            print(a, b, c, d)
            root = algebra.algebraic.ComplexRep(poly, a, b, c, d)
            print(root)
    
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


    M = algebra.matricies.MatrixOver(QQ)


    A = M(4, 4, [[1, 1, 1, 1],
                 [0, 1, 1, 1],
                 [0, 0, 2, 1],
                 [0, 0, 0, 1]])



    print(A)
    print(A.char_mat())
    print(A.char_mat().smith_normal_form())
    print(A.min_poly())
    print(A.char_poly())


    input("done")


    P = algebra.polynomials.PolyOver(QQ)
    x = P.var()

    a = 12 * (x ** 5 - x ** 4 - 2 * x ** 3 - 8 * x ** 2 + 6 * x - 1) ** 10 * (x + 3) ** 6
    #a = 12 * (x + 1) ** 3
    print(a.factor())

test2()
