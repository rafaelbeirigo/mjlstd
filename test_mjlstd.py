import mjlstd
import unittest
import MJLS
from test_cases import eye_one_constants
from test_cases import eye_two_constants
import sam_constants
import numpy as np
import numpy.testing as npt
import math


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
        npt.assert_array_almost_equal(mjlstd.get_F(self.mjls_obj,
                                                   self.cf.F_0.copy(),
                                                   self.cf.X_0),
                                      self.cf.F_0)

    def test_get_Upsilon(self):
        Upsilon = math.inf * np.ones_like(self.cf.Upsilon_1)
        for i in range(self.cf.N):
            Upsilon[i] = mjlstd.get_Upsilon(m=self.mjls_obj,
                                            Fs=self.cf.F_0,
                                            Upsilon=self.cf.Upsilon_0[i],
                                            i=i)

        npt.assert_array_almost_equal(Upsilon, self.cf.Upsilon_1)

    def test_get_D(self):
        npt.assert_array_almost_equal(
            mjlstd.get_D(m=self.mjls_obj,
                         Fs=self.cf.F_0,
                         Ys=self.cf.X_0,
                         Upsilon=self.cf.Upsilon_0[self.cf.i1],
                         i=self.cf.i1,
                         j=self.cf.i2),
            self.cf.D_cal_0)


class TestMjlstdEyeTwo(TestMjlstd):
    def setUp(self):
        self.cf = eye_two_constants
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


class TestMjlstdSam(TestMjlstd):
    def setUp(self):
        self.cf = sam_constants
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
