from Error import Not3DError, MatricesNumberError
import numpy as np


class Matrices:
    """Matrices corresponding to each of an MJLS's A, B, C, D, F, and X.
    """
    def __init__(self, N, M):
        """Constructor for the `Matrices' class.

        Args:
            N (:obj:`int`): number of modes of the MJLS.
            M (:obj:`array`): an array of matrices, one matrix corresponding
                to  each one of the `N' modes.
        """
        if np.ndim(M) != 3:
            raise Not3DError('The array of matrices must have 3 dimensions.')
        if M.shape[0] != N:
            raise MatricesNumberError("The number of"
                                      "matrices must be equal" " to `N'.")

        self.N = N
        self.M = M
