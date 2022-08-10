import algebra
import pyalgebra

ZZ = algebra.base.ZZ
QQ = algebra.base.QQ
G = algebra.base.Gaussian
i = G(0, 1)

P = algebra.polynomials.PolyOver(G)
x = P.var()

print((x + i + 1) ** 3)
