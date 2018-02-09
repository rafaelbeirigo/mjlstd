import math
import numpy as np
import scipy.linalg as la
from Error import MatricesNumberError, DimensionError


class MJLS:
    """A Markov Jump Linear System."""
    def __init__(self, N, m, n, A, B, C, D, P=None, X=None, F=None):
        self.N = N
        self.m = m
        self.n = n
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.P = P
        self.X = X
        self.F = F

    def __get_matrices_maybe(self, N, M):
        """Returns the `N' matrices of `M' if and only if `N' is exact."""
        if M.shape[0] != N:
            raise MatricesNumberError("Each matrix (with the exception"
                                      " of `P') must contain `N' matrices.")
        return M

    @property
    def N(self):
        return self._N

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, value):
        self._m = value

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

    @N.setter
    def N(self, value):
        self._N = value

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, value):
        if (value.shape[1] != self.n) or (value.shape[2] != self.n):
            raise DimensionError("`A' must have `n' rows and `n' columns.")

        self._A = self.__get_matrices_maybe(self.N, value)

    @property
    def B(self):
        return self._B

    @B.setter
    def B(self, value):
        if (value.shape[1] != self.n) or (value.shape[2] != self.m):
            raise DimensionError("`B' must have `n' rows and `m' columns.")

        self._B = self.__get_matrices_maybe(self.N, value)

    @property
    def C(self):
        return self._C

    @C.setter
    def C(self, value):
        if (value.shape[2] != self.n):
            raise DimensionError("`C' must have `n' columns.")

        self._C = self.__get_matrices_maybe(self.N, value)

    @property
    def D(self):
        return self._D

    @D.setter
    def D(self, value):
        if (value.shape[2] != self.m):
            raise DimensionError("`D' must have `m' columns.")

        self._D = self.__get_matrices_maybe(self.N, value)

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, value):
        self._P = value

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, value):
        if (value.shape[1] != self.n) or (value.shape[2] != self.n):
            raise DimensionError("`X' must have `n' rows and `n' columns.")

        self._X = self.__get_matrices_maybe(self.N, value)

    @property
    def F(self):
        return self._F

    @F.setter
    def F(self, value):
        if (value.shape[1] != self.m) or (value.shape[2] != self.n):
            raise DimensionError("`F' must have `m' rows and `n' columns.")

        self._F = self.__get_matrices_maybe(self.N, value)

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


def is_stable(m, lambda_):
    """Tests parameters for TD(\lambda) stability.

    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        lambda_ (:obj:`float`): the $\lambda$ value to be tested.

    The lambda test is based on Lemma 2 from (1).

    The "Fs" test is based on Equations 3.12b--d from (2).

    1. O. L. V. Costa and J. C. C. Aya, "Monte Carlo
    TD(\lambda)-methods for the optimal control of
    discrete-time Markovian jump linear systems",
    Automatica, vol. 38, pp. 217–225, 2002.

    Implement it according to (1, p. 34, equations 3.12[b--d]).

    2. O. L. V. Costa, M. D. Fragoso, and R. P. Marques,
    Discrete-Time Markov Jump Linear Systems, ser.
    Probability and Its Applications. New York:
    Springer-Verlag, 2005.

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
    C_cal = np.kron(m.P.T, np.eye(m.n * m.n))
    N_cal = la.block_diag(*krs)
    A_cal = C_cal.dot(N_cal)

    if not (max(abs(la.eig(A_cal)[0])) < 1):
        return False

    if not (pow(lambda_, 2) < v_max):
        return False

    return True
