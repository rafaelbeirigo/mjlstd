import numpy as np
from scipy.linalg import inv
from sam_run_episode import get_next_theta
import datetime


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


def log_D(k, got_D):
    f = open('D.log', 'a')
    f.write('{: 6} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e}\n'.format(k, *got_D))
    f.close()


def log(t, *args, filename='log.log'):
    line = ''.join(('{: 6}', ' {: 6.1e}' * len(args), '\n')).format(t, *args)
    f = open(filename, 'a')
    f.write(line)
    f.close()


def get_Y(p, m, Fs, Ys):
    """Calculates Y.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    for t in range(1, p.T + 1):
        gamma = p.c * pow(t, -p.eta)
        incr = gamma * get_sum_Ds(p, m, Fs, Ys)
        Ys += incr

        if abs(incr).max() < p.epsilon:
            return Ys

    return Ys


def get_sum_Ds(p, m, Fs, Ys):
    """Calculates the D sum.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    sum_Ds = 0 * Ys.copy()
    for i in range(m.N):
        sum_Ds[i] = get_sum_D(p, m, Fs, Ys, i)

    return sum_Ds


def get_sum_D(p, m, Fs, Ys, i):
    """Calculates each individual D sum.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    theta = [i, 0]
    sum_D = 0 * Ys[0].copy()
    Upsilon = np.eye(m.A.shape[1])
    for k in range(p.K - 1):
        theta[1] = get_next_theta(theta[0], m.P)
        incr = pow(p.lambda_, k) * get_D(m, Fs, Ys, Upsilon, theta)
        sum_D += incr

        if abs(incr).max() < p.epsilon:
            return sum_D

        Upsilon = get_Upsilon(m, Fs, theta, Upsilon)
        theta[0] = theta[1]

    return sum_D


def get_D(m, Fs, Ys, Upsilon, theta):
    """Calculates each individual D for the sum.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    (A, B, C, D) = m.get_ABCD(theta[1])

    C_, D_ = C.conj().T, D.conj().T

    F = Fs[theta[1]]
    F_ = F.conj().T

    Y1, Y2 = Ys[theta[0]], Ys[theta[1]]

    U = Upsilon
    U_ = U.conj().T

    B_cal = C_.dot(C) + F_.dot(D_.dot(D.dot(F)))
    C_cal = ((A + B.dot(F)).conj().T).dot(Y2.dot(A + B.dot(F))) - Y1
    D_cal = U_.dot((B_cal + C_cal).dot(U))

    return D_cal


def get_Upsilon(m, Fs, theta, Upsilon):
    """Calculate current Upsilon
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    (A, B, _, _) = m.get_ABCD(theta[0])
    F = Fs[theta[0]]

    return ((A + B.dot(F))).dot(Upsilon)


def get_F(m, Fs, Ys):
    """Calculate F.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    for i in range(m.N):
        (A, B, _, D) = m.get_ABCD(i)

        B_, D_ = B.conj().T, D.conj().T

        Y = Ys[i]

        Fs[i] = (-inv(B_.dot(Y).dot(B) + D_.dot(D))).dot(B_).dot(Y).dot(A)

    return Fs


def mjlstd(p, m):
    """Applies the TD(\lambda) method to solve a MJLS.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    if m.X is None:
        m.X = 0 * m.A.copy()
    if m.F is None:
        m.F = np.array([0 * B.T for B in m.B])

    Ys, Fs = m.X.copy(), m.F.copy()

    np.random.seed(p.seed)

    filename = ''.join((datetime.datetime.now().strftime("%Y%m%d%H%M%S"),
                        '.log'))
    for l in range(p.L):
        # Log
        Ys_old = Ys.copy()
        Fs_old = Fs.copy()

        # Calculate updated Ys and Fs
        Ys = get_Y(p, m, Fs, Ys)
        Fs = get_F(m, Fs, Ys)

        # Log
        err_Ys = abs(Ys_old - Ys).max()
        err_Fs = abs(Fs_old - Fs).max()

        err_Ys_par = abs(Ys - m.X).max()
        err_Fs_par = abs(Fs - m.F).max()

        print('err_Ys: {:5.3e}  err_Fs: {:5.3e}  err_Ys_par: {:5.3e}  '
              'err_Fs_par: {:5.3e}'.format(err_Ys, err_Fs, err_Ys_par,
                                           err_Fs_par))

    return (Fs, Ys)
