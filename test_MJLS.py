import sam_constants as sc
import MJLS
import numpy as np
import unittest


class TestMJLS(unittest.TestCase):
    def setUp(self):
        self.good_boy = MJLS.MJLS(sc.As, sc.Bs, sc.Cs,
                                  sc.Ds, sc.P, sc.X_ric, sc.F_ric)
        self.good_boy = MJLS.MJLS(sc.N, sc.As, sc.Bs, sc.Cs,
                                  sc.Ds, sc.P, sc.Xs_ric, sc.Fs_ric)

    def tearDown(self):
        pass

    def test_number_of_dimensions(self):
        pass
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
