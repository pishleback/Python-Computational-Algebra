import algebra
import pyalgebra

x = algebra.base.ZZ.Factorization(1, {2 : 3})
y = algebra.base.ZZ.Factorization(1, {-2 : 4})

print(x)
print(y)
print(x * y)
