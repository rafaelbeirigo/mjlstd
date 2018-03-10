import numpy as np
from scipy.linalg import inv
from sam_run_episode import get_next_theta
import random


def log_D(k, got_D):
    f = open('D.log', 'a')
    f.write('{: 6} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e}\n'.format(k, *got_D))
    f.close()


def log(t, *args, filename='log.log'):
    line = ''.join(('{: 6}', ' {: 11.6e}' * len(args), '\n')).format(t, *args)
    f = open(filename, 'a')
    f.write(line)
    f.close()


def get_Y(p, m, Fs, Ys, Ys_hist):
    """Calculates Y.
    Args:
        p (:obj:`Parameters`): parameters for the algorithm; for details on
            each parameter, see the :obj:`Parameters` class.
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
        Ys_hist: the `Ys' history (i.e., a sequence with the `Ys' calculated
             so far.)
    """
    gamma = 1.0       # TODO: do we really want it constant?
    lambda_ = 0.1     # TODO: do we really want it constant?

    e = np.zeros(m.N)

    for t in range(1, p.T + 1):
        Upsilon = np.eye(m.n)

        # TODO: parametrize
        theta = random.randint(0, m.N - 1)

        for k in range(1, p.K):
            alpha = 0.1 / k

            next_theta = get_next_theta(theta, m.P)

            delta = get_delta(m, Fs, Ys, Upsilon, gamma, theta, next_theta)

            e[theta] += 1

            for i in range(m.N):
                Ys[i] += alpha * e[i] * delta
                e[i] *= gamma * lambda_

            Upsilon = get_Upsilon(m, Fs, Upsilon, theta)

            theta = next_theta

        Ys_hist.append(Ys.copy())

    return (Ys, Ys_hist)


def get_delta(m, Fs, Ys, Upsilon, gamma, i, j):
    """Calculates the temporal difference `delta'.
    Args:
        m (:obj:`MJLS`): the corresponding Markov Jump Linear System.
        Fs: current approximation of the control gains.
        Ys: current approximation of the CARE solution.
        Upsilon: current value of `Upsilon'.
        gamma: discount rate.
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

    r = C_.dot(C) + F_.dot(D_).dot(D).dot(F)
    V1 = Y1
    V2 = ((A + B.dot(F)).conj().T).dot(Y2).dot(A + B.dot(F))

    delta = U_.dot(r + gamma * V2 - V1).dot(U)

    return delta


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

        DD = D_.dot(D)

        # Insert a "identity-like" row if it is all zeros
        for row in range(DD.shape[0]):
            if not np.any(DD[row]):
                DD[row][row] = 1.0

        Fs[i] = (-inv(B_.dot(Y).dot(B) + DD)).dot(B_).dot(Y).dot(A)

    return Fs


def mjlstd_eligibility(p, m):
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
    Ys_H, Fs_H = [Ys], [Fs]

    np.random.seed(p.seed)

    for l in range(p.L):
        print('mjlstd_eligibility.py: L {:3d} of {:3d} '
              '({:3.0f}%)'.format(l + 1, p.L, 100. * (l + 1)/p.L))

        (Ys, Ys_h) = get_Y(p, m, Fs.copy(), Ys.copy(), [])
        Fs = get_F(m, Fs.copy(), Ys.copy())

        Ys_H.extend(Ys_h)
        Fs_H.append(Fs)

    return (Fs, Ys, Fs_H, Ys_H)
