import mjlstd
import unittest
import MJLS
from test_cases import eye_one_constants
import numpy.testing as npt


class TestMjlstd(unittest.TestCase):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = eye_one_constants
        args = {'N': self.cf.N,
                'm': self.cf.m,
                'n': self.cf.n,
                'A': self.cf.A,
                'B': self.cf.B,
                'C': self.cf.C,
                'D': self.cf.D,
                'P': self.cf.P,
                'X': self.cf.X,
                'F': self.cf.F}
        self.mjls_obj = MJLS.MJLS(**args)

    def test_get_F(self):
        npt.assert_array_equal(mjlstd.get_F(self.mjls_obj,
                                            self.cf.F_0.copy(),
                                            self.cf.X_0),
                               self.cf.F_0)
