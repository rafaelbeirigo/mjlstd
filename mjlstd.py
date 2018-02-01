import numpy as np
from scipy.linalg import inv
from sam_run_episode import get_next_theta


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


def log(t, sum_D):
    f = open('sum_D.log', 'a')
    f.write('{: 6} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e}\n'.format(t, *sum_D))
    f.close()


def get_Y(m, Fs, Ys, Theta, T, K, c, eta, lambda_par, epsilon):
    """Calculates Y.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    for t in range(1, T+1):
        Ys_old = Ys.copy()

        sum_D = get_sum_D(m, Fs, Ys, Theta, K, lambda_par, epsilon)

        gamma = c * pow(t, -eta)
        Ys = Ys + gamma * sum_D

        err = abs(Ys_old - Ys).max()
        if err < epsilon:
            break

    return Ys


def get_sum_D(m, Fs, Ys, Theta, K, lambda_par, epsilon):
    """Calculates the D sum.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    sum_D = Ys.copy()
    Upsilons = Ys.copy()
    for i in range(m.N):
        sum_D[i] *= 0
        Upsilons[i] = np.eye(m.A.shape[1])
        for k in range(K-1):
            Theta[i, k+1] = get_next_theta(Theta[i, k], m.P)

            sum_D_old = sum_D[i].copy()
            got_D = get_D(m, Fs, Ys, Upsilons, Theta, i, k)
            sum_D[i] += pow(lambda_par, k) * got_D

            if abs(sum_D_old - sum_D[i]).max() < epsilon:
                break

            Upsilons[i] = get_Upsilon(m, Fs, Theta, Upsilons, i, k)

    return sum_D


def get_D(m, Fs, Ys, Upsilons, Theta, i, k):
    """Calculates each individual D for the sum.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    (A, B, C, D) = m.get_ABCD(Theta[i, k+1])

    C_ = C.conj().T
    D_ = D.conj().T

    F = Fs[Theta[i, k+1]]
    F_ = F.conj().T

    Y1 = Ys[Theta[i, k]]
    Y2 = Ys[Theta[i, k+1]]

    U = Upsilons[i]
    U_ = U.conj().T

    B_cal = C_.dot(C) + F_.dot(D_.dot(D.dot(F)))
    C_cal = ((A + B.dot(F)).conj().T).dot(Y2.dot(A + B.dot(F))) - Y1
    D_cal = U_.dot((B_cal + C_cal).dot(U))

    return D_cal


def get_Upsilon(m, Fs, Theta, Upsilons, i, k):
    """Calculate current Upsilon
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    (A, B, _, _) = m.get_ABCD(Theta[i, k])

    F = Fs[Theta[i, k]]

    return ((A + B.dot(F))).dot(Upsilons[i])


def get_F(m, Fs, Ys):
    """Calculate F.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    for i in range(m.N):
        (A, B, _, D) = m.get_ABCD(i)

        B_ = B.conj().T
        D_ = D.conj().T

        Y = Ys[i]

        Fs[i] = (-inv(B_.dot(Y).dot(B) + D_.dot(D))).dot(B_).dot(Y).dot(A)

    return Fs


def mjlstd(m, lambda_par, L, T, K, epsilon, seed, c, eta):
    """Applies the TD(\lambda) method to solve a MJLS.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
    """
    if m.X is not None:
        Ys = m.X.copy()
    if m.F is not None:
        Fs = m.F.copy()

    Theta = np.zeros((m.N, K), dtype=int)
    Theta[:, 0] = [i for i in range(m.N)]

    np.random.seed(seed)

    for _ in range(L):
        # Log
        Ys_old = Ys.copy()
        Fs_old = Fs.copy()

        # Calculate updated Ys and Fs
        Ys = get_Y(m, Fs, Ys, Theta, T, K, c, eta, lambda_par, epsilon)
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
