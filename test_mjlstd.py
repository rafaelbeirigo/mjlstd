import mjlstd
import unittest
import MJLS
from test_cases import eye_one_constants
from test_cases import eye_one_constants_F_eye
from test_cases import eye_one_parameters
from test_cases import eye_one_parameters_bigger_K
from test_cases import eye_two_constants
from test_cases import eye_two_constants_F_eye
from test_cases import eye_two_parameters_bigger_K
from test_cases import zero_one_constants
from test_cases import zero_one_parameters
from test_cases import zero_two_constants
from test_cases import zero_two_parameters
from test_cases import eye_one_constants_mjlstd
from test_cases import eye_one_parameters_mjlstd
import sam_constants
import numpy as np
import numpy.testing as npt
import math
import Parameters


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


class TestMjlstdWithParameters(unittest.TestCase):
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

        # The (p)arameters (f)ile
        self.pf = eye_one_parameters
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)

    def test_get_sum_D(self):
        npt.assert_array_almost_equal(mjlstd.get_sum_D(p=self.params_obj,
                                                       m=self.mjls_obj,
                                                       Fs=self.cf.F_get_sum_D,
                                                       Ys=self.cf.X_get_sum_D,
                                                       i=self.pf.i),
                                      self.pf.sum_D)

    def test_get_sum_Ds(self):
        npt.assert_array_almost_equal(
            mjlstd.get_sum_Ds(p=self.params_obj,
                              m=self.mjls_obj,
                              Fs=self.cf.F_get_sum_D,
                              Ys=self.cf.X_get_sum_D),
            self.pf.sum_Ds)

    def test_get_Y(self):
        npt.assert_array_almost_equal(
            mjlstd.get_Y(p=self.params_obj,
                         m=self.mjls_obj,
                         Fs=self.cf.F_get_Y,
                         Ys=self.cf.X_get_Y),
            self.cf.got_Y)


class TestMjlstdWithParametersAllZeroOneD(TestMjlstdWithParameters):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = zero_one_constants
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

        # The (p)arameters (f)ile
        self.pf = zero_one_parameters
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)


class TestMjlstdWithParametersAllZeroTwoD(TestMjlstdWithParameters):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = zero_two_constants
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

        # The (p)arameters (f)ile
        self.pf = zero_two_parameters
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)


class TestMjlstdWithParametersBiggerK(TestMjlstdWithParameters):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = eye_one_constants_F_eye
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

        # The (p)arameters (f)ile
        self.pf = eye_one_parameters_bigger_K
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)


class TestMjlstdWithParametersBiggerKEyeTwo(TestMjlstdWithParameters):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = eye_two_constants_F_eye
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

        # The (p)arameters (f)ile
        self.pf = eye_two_parameters_bigger_K
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)


class TestMjlstdMjlstd(unittest.TestCase):
    def setUp(self):
        # The (c)onstants (f)ile
        self.cf = eye_one_constants_mjlstd
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

        # The (p)arameters (f)ile
        self.pf = eye_one_parameters_mjlstd
        args_p = {'L': self.pf.L,
                  'T': self.pf.T,
                  'K': self.pf.K,
                  'lambda_': self.pf.lambda_,
                  'epsilon': self.pf.epsilon,
                  'c': self.pf.c,
                  'eta': self.pf.eta,
                  'seed': self.pf.seed}
        self.params_obj = Parameters.Parameters(**args_p)

    def test_mjlstd(self):
        npt.assert_array_almost_equal(
            mjlstd.mjlstd(p=self.params_obj,
                          m=self.mjls_obj),
            self.cf.mjlstd_F_Y)
