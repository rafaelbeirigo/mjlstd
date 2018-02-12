class Parameters:
    """Parameters for the MJLS TD(\lambda) algorithm."""
    def __init__(self, L, T, K, lambda_, epsilon, c, eta, seed):
        """
        Args:
            L (:obj:`int`): max iterations for calculating Y.
            T (:obj:`int`): max interations for calculating each Y candidate.
            K (:obj:`int`): max iterations for calculating the sum used to
                calculate each Y candidate.
            lambda_ (:obj:`float`): coefficient for the sum used to
                calculate each Y candidate.
            epsilon (:obj:`float`): minimum difference between two numbers to
                consider them equal when testing for convergence.
            c (:obj:`float`): coefficient used to calculate the step size.
            eta (:obj:`float`): exponent used to calculate the step size.
            seed (:obj:`int`): value used to initialize the random number
                generator.
        """
        self.L = L
        self.T = T
        self.K = K
        self.lambda_ = lambda_
        self.epsilon = epsilon
        self.c = c
        self.eta = eta
        self.seed = seed

    def positive_or_error(self, value, name):
        """Raises an error if value is not positive."""

        if value < 0:
            raise ValueError(''.join((name, " must be positive.")))

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, value):
        self.positive_or_error(value, "`L'")
        self._L = value

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value):
        self.positive_or_error(value, "`T'")
        self._T = value

    @property
    def K(self):
        return self._K

    @K.setter
    def K(self, value):
        self.positive_or_error(value, "`K'")
        self._K = value

    @property
    def lambda_(self):
        return self._lambda_

    @lambda_.setter
    def lambda_(self, value):
        self.positive_or_error(value, "`lambda_'")
        self._lambda_ = value

    @property
    def epsilon(self):
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        self.positive_or_error(value, "`epsilon'")
        self._epsilon = value

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, value):
        self.positive_or_error(value, "`c'")
        self._c = value

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, value):
        self.positive_or_error(value, "`eta'")
        self._eta = value

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, value):
        self.positive_or_error(value, "`seed'")
        self._seed = value
