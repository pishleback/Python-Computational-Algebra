import algebra
import pyalgebra

ZZ = algebra.base.ZZ
QQ = algebra.base.ZZ.FractionField
G = algebra.base.Gaussian
i = G(0, 1)



P = algebra.polynomials.PolyOver(QQ)
x = P.var()

a = 12 * (x ** 5 - x ** 4 - 2 * x ** 3 - 8 * x ** 2 + 6 * x - 1) ** 10 * (x + 3) ** 6
#a = 12 * (x + 1) ** 3
print(a.factor())
