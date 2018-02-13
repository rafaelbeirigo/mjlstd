from Error import MatricesNumberError, DimensionError
import sam_constants
import eye_one_constants
import MJLS
import unittest
import numpy.testing as npt
import math


class TestMJLS(unittest.TestCase):
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

    def test_N_matrices_per_matrix(self):
        args = {'N': math.inf,
                'm': self.cf.m,
                'n': self.cf.n,
                'A': self.cf.A,
                'B': self.cf.B,
                'C': self.cf.C,
                'D': self.cf.D,
                'P': self.cf.P,
                'X': self.cf.X,
                'F': self.cf.F}
        self.assertRaises(MatricesNumberError, MJLS.MJLS, **args)

    def test_get_and_set_methods(self):
        self.assertIs(self.cf.N, self.mjls_obj.N)
        self.assertEqual(self.cf.m, self.mjls_obj.m)
        self.assertEqual(self.cf.n, self.mjls_obj.n)

        self.assertIs(self.cf.A, self.mjls_obj.A)
        self.assertIs(self.cf.B, self.mjls_obj.B)
        self.assertIs(self.cf.C, self.mjls_obj.C)
        self.assertIs(self.cf.D, self.mjls_obj.D)
        self.assertIs(self.cf.X, self.mjls_obj.X)
        self.assertIs(self.cf.F, self.mjls_obj.F)

    def test_get_bundle_ABCDFX(self):
        ABCD = self.mjls_obj.get_ABCD()
        self.assertIs(self.cf.A, ABCD[0])
        self.assertIs(self.cf.B, ABCD[1])
        self.assertIs(self.cf.C, ABCD[2])
        self.assertIs(self.cf.D, ABCD[3])

    def test_get_individual_ABCDFX(self):
        for i in range(self.cf.N):
            ABCD = self.mjls_obj.get_ABCD(i)
            npt.assert_array_equal(self.cf.A[i], ABCD[0])
            npt.assert_array_equal(self.cf.B[i], ABCD[1])
            npt.assert_array_equal(self.cf.C[i], ABCD[2])
            npt.assert_array_equal(self.cf.D[i], ABCD[3])

    def test_m_is_dimensionally_respected(self):
        args = {'N': self.cf.N,
                'm': math.inf,
                'n': self.cf.n,
                'A': self.cf.A,
                'B': self.cf.B,
                'C': self.cf.C,
                'D': self.cf.D,
                'P': self.cf.P,
                'X': self.cf.X,
                'F': self.cf.F}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)

    def test_n_is_dimensionally_respected(self):
        args = {'N': self.cf.N,
                'm': self.cf.m,
                'n': math.inf,
                'A': self.cf.A,
                'B': self.cf.B,
                'C': self.cf.C,
                'D': self.cf.D,
                'P': self.cf.P,
                'X': self.cf.X,
                'F': self.cf.F}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)

    def test_correct_m_and_n_are_ok(self):
        args = {'N': self.cf.N,
                'm': self.cf.m,      # ok
                'n': self.cf.n,      # ok
                'A': self.cf.A,
                'B': self.cf.B,
                'C': self.cf.C,
                'D': self.cf.D,
                'P': self.cf.P,
                'X': self.cf.X,
                'F': self.cf.F}
        self.assertIsNotNone(MJLS.MJLS(**args))

    def test_N_must_be_positive(self):
        args_N = {'N': -1, 'm': self.cf.m, 'n': self.cf.n,
                  'A': self.cf.A, 'B': self.cf.B, 'C':
                  self.cf.C, 'D': self.cf.D, 'P': self.cf.P,
                  'X': self.cf.X, 'F': self.cf.F}
        self.assertRaises(ValueError, MJLS.MJLS, **args_N)

    def test_m_must_be_positive(self):
        args_m = {'N': self.cf.N, 'm': -1, 'n': self.cf.n,
                  'A': self.cf.A, 'B': self.cf.B, 'C':
                  self.cf.C, 'D': self.cf.D, 'P': self.cf.P,
                  'X': self.cf.X, 'F': self.cf.F}
        self.assertRaises(ValueError, MJLS.MJLS, **args_m)

    def test_n_must_be_positive(self):
        args_n = {'N': self.cf.N, 'm': self.cf.m, 'n': -1,
                  'A': self.cf.A, 'B': self.cf.B, 'C':
                  self.cf.C, 'D': self.cf.D, 'P': self.cf.P,
                  'X': self.cf.X, 'F': self.cf.F}
        self.assertRaises(ValueError, MJLS.MJLS, **args_n)

    def test_P_has_only_non_negative_values(self):
        args_P = {'N': self.cf.N, 'm': self.cf.m, 'n':
                  self.cf.n, 'A': self.cf.A, 'B': self.cf.B,
                  'C': self.cf.C, 'D': self.cf.D, 'P':
                  -self.cf.P, 'X': self.cf.X, 'F':
                  self.cf.F}
        self.assertRaises(ValueError, MJLS.MJLS, **args_P)

    def test_each_row_of_P_sum_to_1(self):
        args_P = {'N': self.cf.N, 'm': self.cf.m, 'n':
                  self.cf.n, 'A': self.cf.A, 'B': self.cf.B,
                  'C': self.cf.C, 'D': self.cf.D, 'P': 0 *
                  self.cf.P, 'X': self.cf.X, 'F':
                  self.cf.F}
        self.assertRaises(ValueError, MJLS.MJLS, **args_P)


class TestMJLSFIsProvided(unittest.TestCase):
    def setUp(self):
        # The (c)onstants (f)ile
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

    def test_lambda_0_makes_it_stable(self):
        self.assertTrue(self.mjls_obj.is_lambda_stable(0))

    def test_lambda_1_makes_it_unstable(self):
        self.assertFalse(self.mjls_obj.is_lambda_stable(1))

    def test_F_stable(self):
        self.assertTrue(self.mjls_obj.is_F_stable())


if __name__ == "__main__":
    unittest.main()
