import os
import numpy as np
from scipy.linalg import inv
from sam_run_episode import get_next_theta
import datetime


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
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
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
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
    """
    sum_Ds = np.zeros_like(Ys)
    for i in range(m.N):
        sum_Ds[i] = get_sum_D(p, m, Fs, Ys, i)

    return sum_Ds


def get_sum_D(p, m, Fs, Ys, i):
    """Calculates each individual D sum.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
        i: starting mode, i.e., the mode for the first time step in the
            current simulation.
    """
    sum_D = np.zeros_like(Ys[0])
    Upsilon = np.eye(m.m)
    theta = i
    for k in range(p.K + 1):
        next_theta = get_next_theta(theta, m.P)
        incr = pow(p.lambda_, k) * get_D(m, Fs, Ys, Upsilon,
                                         theta, next_theta)
        sum_D += incr

        if abs(incr).max() < p.epsilon:
            return sum_D

        Upsilon = get_Upsilon(m, Fs, Upsilon, theta)
        theta = next_theta

    return sum_D


def get_D(m, Fs, Ys, Upsilon, i, j):
    """Calculates each individual D for the sum.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
        Upsilon: current value of `Upsilon'.
        i: current mode.
        j: next mode.
    """
    (A, B, C, D) = m.get_ABCD(j)

    C_, D_ = C.conj().T, D.conj().T

    F = Fs[j]
    F_ = F.conj().T

    Y1, Y2 = Ys[i], Ys[j]

    U = Upsilon
    U_ = U.conj().T

    B_cal = C_.dot(C) + F_.dot(D_).dot(D).dot(F)
    C_cal = ((A + B.dot(F)).conj().T).dot(Y2).dot(A + B.dot(F)) - Y1
    D_cal = U_.dot(B_cal + C_cal).dot(U)

    return D_cal


def get_Upsilon(m, Fs, Upsilon, i):
    """Calculate current Upsilon
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Upsilon: current value of `Upsilon'.
        i: current mode.
    """
    (A, B, _, _) = m.get_ABCD(i)
    F = Fs[i]

    return ((A + B.dot(F))).dot(Upsilon)


def get_F(m, Fs, Ys):
    """Calculate F.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
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
        m.X = np.zeros_like(m.A)
    if m.F is None:
        m.F = np.array([np.zeros_like(B.T) for B in m.B])

    Ys, Fs = m.X.copy(), m.F.copy()

    np.random.seed(p.seed)

    datetime_string = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join("log", ''.join((datetime_string, ".log")))

    for l in range(p.L):
        Ys = get_Y(p, m, Fs, Ys)
        Fs = get_F(m, Fs, Ys)

        log(l, *(Ys - m.X).flatten(), filename=filename)

    return (Fs, Ys)
