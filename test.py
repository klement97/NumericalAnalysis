import unittest
import math
import single_equations as se


class SingleEquationsTestCase(unittest.TestCase):

    def test_bisection(self):
        self.assertAlmostEqual(se.bisection(f=lambda x: x + math.exp(x) - 2,
                                            a=0, b=1),
                               0.442855, 5)

    def test_fixed_point(self):
        self.assertAlmostEqual(se.fixed_point(f=lambda x: math.pi + math.atan(x),
                                              x0=math.pi),
                               4.493410, 5)

    def test_lagrange(self):
        self.assertAlmostEqual(se.lagrange(f=lambda x: x + math.exp(x) - 2,
                                           a=0, b=1),
                               0.442855, 5)

    def test_newton(self):
        self.assertAlmostEqual(se.newton(f=lambda x: x + math.exp(x) - 2,
                                         df=lambda x: 1+math.exp(x),
                                         x0=0),
                               0.442855, 5)

    def test_combined(self):
        self.assertAlmostEqual(se.combined(f=lambda x: x + math.exp(x) - 2,
                                           df=lambda x: 1+math.exp(x),
                                           a=0, b=1),
                               0.442855, 5)

    def test_secant(self):
        self.assertAlmostEqual(se.secant(f=lambda x: x + math.exp(x) - 2, x0=0, x1=1),
                               0.442855, 5)

    def test_dekker_brent(self):
        self.assertAlmostEqual(se.dekker_brent(f=lambda x: x + math.exp(x) - 2, x0=0, x1=1),
                               0.442855, 5)

