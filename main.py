import algebra
import pyalgebra

G = algebra.base.Gaussian


n = G(13, 0).factor().list()[0]
Q = G.QuotientRing(n)
x = Q.convert(G(1, 1)) + G(1, 1)
print(x)
