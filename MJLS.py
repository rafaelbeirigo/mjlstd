import math
import numpy as np
import scipy.linalg as la


class MJLS:
    """A Markov Jump Linear System."""
    def __init__(self, A, B, C, D, P=None, X=None, F=None):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.P = P
        self.X = X
        self.F = F
        self.N = A.shape[0]

    def get_ABCD(self, i=None):
        """Returns :obj: `A`, `B', `C', and `D' matrices.

        Returns all `A''s, `B''s, `C''s and `D''s, or the
        `i'th one, if `i' is provided.

        Args:
            i (:obj:`int`, optional): index used to retrive
                the matrices.
        """
        if(i is not None):
            return (self.A[i], self.B[i], self.C[i], self.D[i])
        else:
            return (self.A, self.B, self.C, self.D)


def test_stability(m, lambd):
    """Tests parameters for TD(\lambda) stability.

    Based on the article:

    O. L. V. Costa and J. C. C. Aya, "Monte Carlo
    TD(\lambda)-methods for the optimal control of
    discrete-time Markovian jump linear systems",
    Automatica, vol. 38, pp. 217–225, 2002.

    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        lambd (:obj:`float`): the $\lambda$ value to be tested.

    """
    assert m.F is not None, "F must be provided"

    krs = []
    v_max = -math.inf
    for i in range(m.N):
        G = m.A[i] + m.B[i].dot(m.F[i])
        kr = np.kron(G, G)
        krs.append(kr)

        # lambda test
        v = pow(la.norm(kr, ord=2), -1)
        if v > v_max:
            v_max = v

    # Fs test
    I_cal = la.block_diag(*krs)

    if not (max(abs(la.eig(I_cal)[0])) < 1):
        print('test_stability: WARNING: F does not satisfy Lemma 3 '
              '(stabilizability)')

    if not (pow(lambd, 2) < v_max):
        print('test_stability: WARNING: lambda does not satisfy Lemma 2 '
              '(convergence)')
