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
