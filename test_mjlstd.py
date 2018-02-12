import unittest
import sam_parameters as sp
import sam_constants as sc
import MJLS


class TestMjlstd(unittest.TestCase):
    def test_Parameters_accept_only_reasonable_values(self):
        args_N = {'N': -1, 'm': sc.m, 'n': sc.n, 'A': sc.As,
                  'B': sc.Bs, 'C': sc.Cs, 'D': sc.Ds, 'P':
                  sc.P, 'X': sc.Xs_ric, 'F': sc.Fs_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_N)

        args_m = {'N': sc.N, 'm': -1, 'n': sc.n, 'A': sc.As,
                  'B': sc.Bs, 'C': sc.Cs, 'D': sc.Ds, 'P': sc.P, 'X':
                  sc.Xs_ric, 'F': sc.Fs_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_m)

        args_n = {'N': sc.N, 'm': sc.m, 'n': -1, 'A': sc.As,
                  'B': sc.Bs, 'C': sc.Cs, 'D': sc.Ds, 'P': sc.P, 'X':
                  sc.Xs_ric, 'F': sc.Fs_ric}
        self.assertRaises(ValueError, MJLS.MJLS, **args_n)
