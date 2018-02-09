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
                'A': sc.As,
                'B': sc.Bs,
                'C': sc.Cs,
                'D': sc.Ds,
                'P': sc.P,
                'X': sc.Xs_ric,
                'F': sc.Fs_ric}
        self.good_boy = MJLS.MJLS(**args)

    def tearDown(self):
        pass

    def test_N_matrices_per_matrix(self):
        args = {'N': 0,
                'm': sc.m,
                'n': sc.n,
                'A': sc.As,
                'B': sc.Bs,
                'C': sc.Cs,
                'D': sc.Ds,
                'P': sc.P,
                'X': sc.Xs_ric,
                'F': sc.Fs_ric}
        self.assertRaises(MatricesNumberError, MJLS.MJLS, **args)

    def test_get_and_set_methods(self):
        self.assertIs(sc.N, self.good_boy.N)
        self.assertEqual(sc.m, self.good_boy.m)
        self.assertEqual(sc.n, self.good_boy.n)

        self.assertIs(sc.As, self.good_boy.A)
        self.assertIs(sc.Bs, self.good_boy.B)
        self.assertIs(sc.Cs, self.good_boy.C)
        self.assertIs(sc.Ds, self.good_boy.D)
        self.assertIs(sc.Xs_ric, self.good_boy.X)
        self.assertIs(sc.Fs_ric, self.good_boy.F)

    def test_get_bundle_ABCDFX(self):
        ABCD = self.good_boy.get_ABCD()
        self.assertIs(sc.As, ABCD[0])
        self.assertIs(sc.Bs, ABCD[1])
        self.assertIs(sc.Cs, ABCD[2])
        self.assertIs(sc.Ds, ABCD[3])

    def test_get_individual_ABCDFX(self):
        for i in range(sc.N):
            ABCD = self.good_boy.get_ABCD(i)
            npt.assert_array_equal(sc.As[i], ABCD[0])
            npt.assert_array_equal(sc.Bs[i], ABCD[1])
            npt.assert_array_equal(sc.Cs[i], ABCD[2])
            npt.assert_array_equal(sc.Ds[i], ABCD[3])

    def test_lambda_0_makes_it_stable(self):
        self.assertTrue(MJLS.is_stable(self.good_boy, 0))

    def test_lambda_1_makes_it_unstable(self):
        self.assertFalse(MJLS.is_stable(self.good_boy, 1))

    def test_m_is_dimensionally_respected(self):
        args = {'N': sc.N,
                'm': sc.n,      # swapped
                'n': sc.n,
                'A': sc.As,
                'B': sc.Bs,
                'C': sc.Cs,
                'D': sc.Ds,
                'P': sc.P,
                'X': sc.Xs_ric,
                'F': sc.Fs_ric}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)

    def test_n_is_dimensionally_respected(self):
        args = {'N': sc.N,
                'm': sc.m,
                'n': sc.m,      # swapped
                'A': sc.As,
                'B': sc.Bs,
                'C': sc.Cs,
                'D': sc.Ds,
                'P': sc.P,
                'X': sc.Xs_ric,
                'F': sc.Fs_ric}
        self.assertRaises(DimensionError, MJLS.MJLS, **args)


if __name__ == "__main__":
    unittest.main()
