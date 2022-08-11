import unittest
import algebra
import pyalgebra

class TestBase(unittest.TestCase):
    def test_pyalgebra(self):
        self.assertEqual(pyalgebra.combinatorics.factorial(0), 1)
        self.assertEqual(pyalgebra.combinatorics.factorial(1), 1)
        self.assertEqual(pyalgebra.combinatorics.factorial(2), 2)
        self.assertEqual(pyalgebra.combinatorics.factorial(3), 6)


        self.assertEqual(pyalgebra.combinatorics.binomial(0, -1), 0)
        self.assertEqual(pyalgebra.combinatorics.binomial(0, 0), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(0, 1), 0)

        self.assertEqual(pyalgebra.combinatorics.binomial(1, -1), 0)
        self.assertEqual(pyalgebra.combinatorics.binomial(1, 0), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(1, 1), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(1, 2), 0)

        self.assertEqual(pyalgebra.combinatorics.binomial(2, -1), 0)
        self.assertEqual(pyalgebra.combinatorics.binomial(2, 0), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(2, 1), 2)
        self.assertEqual(pyalgebra.combinatorics.binomial(2, 2), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(2, 3), 0)

        self.assertEqual(pyalgebra.combinatorics.binomial(3, -1), 0)
        self.assertEqual(pyalgebra.combinatorics.binomial(3, 0), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(3, 1), 3)
        self.assertEqual(pyalgebra.combinatorics.binomial(3, 2), 3)
        self.assertEqual(pyalgebra.combinatorics.binomial(3, 3), 1)
        self.assertEqual(pyalgebra.combinatorics.binomial(3, 4), 0)

        

    @unittest.skip
    def test_axioms(self):
        for mset in [algebra.base.ZZ,
                     algebra.base.ZZ.FractionField,
                     algebra.base.Gaussian,
                     algebra.base.PrimQQ,
                     algebra.base.ZZ.QuotientRing(12),
                     algebra.base.ZZ.QuotientRing(13),
                     algebra.polynomials.PolyOver(algebra.base.ZZ),
                     algebra.polynomials.PolyOver(algebra.base.ZZ.FractionField),
                     algebra.polynomials.PolyOver(algebra.base.PrimQQ)]:
            print("Checking", mset, len(mset.test_values()))
            mset.test_axioms(self)
                        
if __name__ == '__main__':
    unittest.main()
