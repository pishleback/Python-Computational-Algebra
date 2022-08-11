import algebra
import pyalgebra




def test2():
    from fractions import Fraction as Frac
    pyalgebra.polynomials.test()
    
    

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
