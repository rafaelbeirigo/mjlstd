import sam_constants as sc
import MJLS
import numpy as np
import unittest


class TestMJLS(unittest.TestCase):
    def setUp(self):
        self.all_none = MJLS.MJLS(None, None, None, None,
                                  None, None, None)
        self.good_boy = MJLS.MJLS(sc.As, sc.Bs, sc.Cs,
                                  sc.Ds, sc.P, sc.X_ric, sc.F_ric)

    def tearDown(self):
        pass

    def test_number_of_dimensions(self):
        pass

if __name__ == "__main__":
    unittest.main()
