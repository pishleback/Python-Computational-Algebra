import unittest
import algebra

class TestBase(unittest.TestCase):
    def test_axioms(self):
        for mset in [algebra.base.ZZ,
                     algebra.base.Gaussian,
                     algebra.base.QQ,
                     algebra.base.ZZ.QuotientRing(12),
                     algebra.base.ZZ.QuotientRing(13),
                     algebra.polynomials.PolyOver(algebra.base.ZZ),
                     algebra.polynomials.PolyOver(algebra.base.QQ)]:
            mset.test_axioms(self)
                        
if __name__ == '__main__':
    unittest.main()
