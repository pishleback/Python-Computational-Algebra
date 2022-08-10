import algebra
import pyalgebra

ZZ = algebra.base.ZZ
QQ = algebra.base.QQ
G = algebra.base.Gaussian
i = G(0, 1)


print(G.all_units())

P = algebra.polynomials.PolyOver(ZZ)
x = P.var()

a = 12 * (x ** 5 - x ** 4 - 2 * x ** 3 - 8 * x ** 2 + 6 * x - 1) * (x + 3) ** 6
#a = 12 * (x + 1) ** 3
print(a.factor())
