from Error import Not3DError, MatricesNumberError
import numpy as np
import unittest
from Matrices import Matrices


class TestMatrices(unittest.TestCase):
    def setUp(self):
        pass

    def test_M_must_be_3D(self):
        self.assertRaises(Not3DError, Matrices, 1, None)
        self.assertRaises(Not3DError, Matrices, 1, np.zeros((1)))
        self.assertRaises(Not3DError, Matrices, 1, np.zeros((1, 1)))
        self.assertRaises(Not3DError, Matrices, 1, np.zeros((1, 1, 1, 1)))

    def test_there_are_N_matrices_in_M(self):
        self.assertRaises(MatricesNumberError,
                          Matrices, 3, np.zeros((1, 1, 1)))
