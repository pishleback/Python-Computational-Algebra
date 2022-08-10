import algebra
import pyalgebra

x = algebra.base.ZZ.Factorization(1, {2 : 3, 3 : 1})
y = algebra.base.ZZ.Factorization(1, {-2 : 5, -5 : 1})


print(algebra.base.ZZ.Factorization.lcm([x, y]))

