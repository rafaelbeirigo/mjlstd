import unittest
import sam_parameters as sp
from mjlstd import Parameters


class TestMjlstd(unittest.TestCase):
    def test_L_must_be_non_negative(self):
        args_L = {'L': -1, 'T': sp.T, 'K': sp.K, 'lambda_':
                  sp.lambda_, 'epsilon': sp.epsilon, 'c': sp.c, 'eta':
                  sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_L)

    def test_T_must_be_non_negative(self):
        args_T = {'L': sp.L, 'T': -1, 'K': sp.K, 'lambda_':
                  sp.lambda_, 'epsilon': sp.epsilon, 'c': sp.c, 'eta':
                  sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_T)

    def test_K_must_be_non_negative(self):
        args_K = {'L': sp.L, 'T': sp.T, 'K': -1, 'lambda_':
                  sp.lambda_, 'epsilon': sp.epsilon, 'c': sp.c, 'eta':
                  sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_K)

    def test_lambda_must_be_non_negative(self):
        args_lambda_ = {'L': sp.L, 'T': sp.T, 'K': sp.K,
                        'lambda_': -1, 'epsilon': sp.epsilon, 'c': sp.c,
                        'eta': sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_lambda_)

    def test_epsilon_must_be_non_negative(self):
        args_epsilon = {'L': sp.L, 'T': sp.T, 'K': sp.K,
                        'lambda_': sp.lambda_, 'epsilon': -1, 'c': sp.c,
                        'eta': sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_epsilon)

    def test_c_must_be_non_negative(self):
        args_c = {'L': sp.L, 'T': sp.T, 'K': sp.K,
                  'lambda_': sp.lambda_, 'epsilon': sp.epsilon, 'c':
                  -1, 'eta': sp.eta, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_c)

    def test_eta_must_be_non_negative(self):
        args_eta = {'L': sp.L, 'T': sp.T, 'K': sp.K,
                    'lambda_': sp.lambda_, 'epsilon': sp.epsilon, 'c':
                    sp.c, 'eta': -1, 'seed': 0}
        self.assertRaises(ValueError, Parameters, **args_eta)

    def test_seed_must_be_non_negative(self):
        args_seed = {'L': sp.L, 'T': sp.T, 'K': sp.K,
                     'lambda_': sp.lambda_, 'epsilon': sp.epsilon, 'c':
                     sp.c, 'eta': sp.eta, 'seed': -1}
        self.assertRaises(ValueError, Parameters, **args_seed)
