from Error import MatricesNumberError, DimensionError
import sam_constants as sc
import MJLS
import unittest
import numpy.testing as npt


class TestMJLS(unittest.TestCase):
    def setUp(self):
        args = {'N': sc.N,
                'm': sc.m,
                'n': sc.n,
                'A': sc.A,
                'B': sc.B,
                'C': sc.C,
                'D': sc.D,
                'P': sc.P,
                'X': sc.X_ric,
                'F': sc.F_ric}
        self.mjls_obj = MJLS.MJLS(**args)

    def tearDown(self):
        pass

    def test_N_matrices_per_matrix(self):
        args = {'N': 1,
                'm': sc.m,
                'n': sc.n,
                'A': sc.A,
                'B': sc.B,
                'C': sc.C,
                'D': sc.D,
                'P': sc.P,
                'X': sc.X_ric,
                'F': sc.F_ric}
        self.assertRaises(MatricesNumberError, MJLS.MJLS, **args)

    def test_get_and_set_methods(self):
        self.assertIs(sc.N, self.mjls_obj.N)
        self.assertEqual(sc.m, self.mjls_obj.m)
        self.assertEqual(sc.n, self.mjls_obj.n)

        self.assertIs(sc.A, self.mjls_obj.A)
        self.assertIs(sc.B, self.mjls_obj.B)
        self.assertIs(sc.C, self.mjls_obj.C)
        self.assertIs(sc.D, self.mjls_obj.D)
        self.assertIs(sc.X_ric, self.mjls_obj.X)
        self.assertIs(sc.F_ric, self.mjls_obj.F)

    def test_get_bundle_ABCDFX(self):
        ABCD = self.mjls_obj.get_ABCD()
        self.assertIs(sc.A, ABCD[0])
        self.assertIs(sc.B, ABCD[1])
        self.assertIs(sc.C, ABCD[2])
        self.assertIs(sc.D, ABCD[3])

    def test_get_individual_ABCDFX(self):
        for i in range(sc.N):
            ABCD = self.mjls_obj.get_ABCD(i)
            npt.assert_array_equal(sc.A[i], ABCD[0])
            npt.assert_array_equal(sc.B[i], ABCD[1])
            npt.assert_array_equal(sc.C[i], ABCD[2])
            npt.assert_array_equal(sc.D[i], ABCD[3])

    def test_lambda_0_makes_it_stable(self):
        self.assertTrue(self.mjls_obj.is_stable(0))

    def test_lambda_1_makes_it_unstable(self):
        self.assertFalse(self.mjls_obj.is_stable(1))

    def test_m_is_dimensionally_respected(self):
        args = {'N': sc.N,
                'm': sc.n,      # swapped
                'n': sc.n,
                'A': sc.A,
                'B': sc.B,
                'C': sc.C,
                'D': sc.D,
                'P': sc.P,
                'X': sc.X_ric,
                'F': sc.F_ric}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)

    def test_n_is_dimensionally_respected(self):
        args = {'N': sc.N,
                'm': sc.m,
                'n': sc.m,      # swapped
                'A': sc.A,
                'B': sc.B,
                'C': sc.C,
                'D': sc.D,
                'P': sc.P,
                'X': sc.X_ric,
                'F': sc.F_ric}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)

    def test_correct_m_and_n_are_ok(self):
        args = {'N': sc.N,
                'm': sc.m,      # ok
                'n': sc.n,      # ok
                'A': sc.A,
                'B': sc.B,
                'C': sc.C,
                'D': sc.D,
                'P': sc.P,
                'X': sc.X_ric,
                'F': sc.F_ric}
        self.assertIsNotNone(MJLS.MJLS(**args))

    def test_N_must_be_positive(self):
        args_N = {'N': -1, 'm': sc.m, 'n': sc.n, 'A': sc.A,
                  'B': sc.B, 'C': sc.C, 'D': sc.D, 'P':
                  sc.P, 'X': sc.X_ric, 'F': sc.F_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_N)

    def test_m_must_be_positive(self):
        args_m = {'N': sc.N, 'm': -1, 'n': sc.n, 'A': sc.A,
                  'B': sc.B, 'C': sc.C, 'D': sc.D, 'P': sc.P, 'X':
                  sc.X_ric, 'F': sc.F_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_m)

    def test_n_must_be_positive(self):
        args_n = {'N': sc.N, 'm': sc.m, 'n': -1, 'A': sc.A,
                  'B': sc.B, 'C': sc.C, 'D': sc.D, 'P': sc.P, 'X':
                  sc.X_ric, 'F': sc.F_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_n)

    def test_P_has_only_non_negative_values(self):
        args_P = {'N': sc.N, 'm': sc.m, 'n': sc.n, 'A':
                  sc.A, 'B': sc.B, 'C': sc.C, 'D': sc.D, 'P':
                  -sc.P, 'X': sc.X_ric, 'F': sc.F_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_P)

    def test_each_row_of_P_sum_to_1(self):
        args_P = {'N': sc.N, 'm': sc.m, 'n': sc.n, 'A':
                  sc.A, 'B': sc.B, 'C': sc.C, 'D': sc.D, 'P': 0 *
                  sc.P, 'X': sc.X_ric, 'F': sc.F_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_P)


if __name__ == "__main__":
    unittest.main()
