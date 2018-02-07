from Error import MatricesNumberError
import sam_constants as sc
import MJLS
import unittest


class TestMJLS(unittest.TestCase):
    def setUp(self):
        args = {'N': sc.N,
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
                'A': sc.As,
                'B': sc.Bs,
                'C': sc.Cs,
                'D': sc.Ds,
                'P': sc.P,
                'X': sc.Xs_ric,
                'F': sc.Fs_ric}
        self.assertRaises(MatricesNumberError, MJLS.MJLS, **args)

    def test_NABCDXF_have_their_set_and_get_working(self):
        self.assertIs(sc.N, self.good_boy.N)
        self.assertIs(sc.As, self.good_boy.A)
        self.assertIs(sc.Bs, self.good_boy.B)
        self.assertIs(sc.Cs, self.good_boy.C)
        self.assertIs(sc.Ds, self.good_boy.D)
        self.assertIs(sc.Xs_ric, self.good_boy.X)
        self.assertIs(sc.Fs_ric, self.good_boy.F)


if __name__ == "__main__":
    unittest.main()
