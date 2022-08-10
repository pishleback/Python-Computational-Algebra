import algebra
import pyalgebra

ZZ = algebra.base.ZZ
QQ = algebra.base.QQ
G = algebra.base.Gaussian
i = G(0, 1)

M = algebra.matricies.MatrixOver(G)

X = M(3, 3, [[1, i, 8],
             [-1, i, 9],
             [3, 3*i, 1]])


print(X)
print(X.det())
input()

S, A, T = X.smith_algorithm()

print(X)
print(S)
print(A)
print(T)
