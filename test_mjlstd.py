import unittest
import sam_parameters as sp
from mjlstd import Parameters


class TestMjlstd(unittest.TestCase):
    def test_L_must_be_non_negative(self):
        args_L = {'L': -1, 'T': sp.T, 'K': sp.K, 'lambda_':
                  sp.lambda_, 'epsilon': sp.epsilon, 'c': sp.c, 'eta':
                  sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_L)
