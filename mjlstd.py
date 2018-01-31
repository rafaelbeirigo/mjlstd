import numpy as np
from scipy.linalg import inv
from sam_run_episode import get_next_theta


def log_D(k, got_D):
    f = open('D.log', 'a')
    f.write('{: 6} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e}\n'.format(k, *got_D))
    f.close()


def log(t, sum_D):
    f = open('sum_D.log', 'a')
    f.write('{: 6} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e} {: 6.1e}\n'.format(t, *sum_D))
    f.close()


def get_Y(As, Bs, Cs, Ds, Fs, P, Ys, Theta, N, T, K, c, eta,
          lambda_par, epsilon):
    """Calculates Y."""

    for t in range(1, T+1):
        Ys_old = Ys.copy()

        sum_D = get_sum_D(As, Bs, Cs, Ds, P, Fs, Ys, Theta, N,
                          K, lambda_par, epsilon)

        gamma = c * pow(t, -eta)
        Ys = Ys + gamma * sum_D

        err = abs(Ys_old - Ys).max()
        if err < epsilon:
            break

    return Ys


def get_sum_D(As, Bs, Cs, Ds, P, Fs, Ys, Theta, N, K, lambda_par, epsilon):
    """Calculates the D sum."""

    sum_D = Ys.copy()
    Upsilons = Ys.copy()
    for i in range(N):
        sum_D[i] *= 0
        Upsilons[i] = np.eye(As.shape[1])
        for k in range(K-1):
            Theta[i, k+1] = get_next_theta(Theta[i, k], P)

            sum_D_old = sum_D[i].copy()
            got_D = get_D(As, Bs, Cs, Ds, Fs, Ys, Upsilons, Theta, i, k)
            sum_D[i] += pow(lambda_par, k) * got_D

            if abs(sum_D_old - sum_D[i]).max() < epsilon:
                break

            Upsilons[i] = get_Upsilon(As, Bs, Fs, Theta, Upsilons, i, k)

    return sum_D


def get_D(As, Bs, Cs, Ds, Fs, Ys, Upsilons, Theta, i, k):
    """Calculates each individual D for the sum."""

    A = As[Theta[i, k+1]]
    B = Bs[Theta[i, k+1]]
    C = Cs[Theta[i, k+1]]
    D = Ds[Theta[i, k+1]]
    F = Fs[Theta[i, k+1]]

    C_ = C.conj().T
    D_ = D.conj().T
    F_ = F.conj().T

    U = Upsilons[i]
    U_ = U.conj().T

    Y1 = Ys[Theta[i, k]]
    Y2 = Ys[Theta[i, k+1]]

    B_cal = C_.dot(C) + F_.dot(D_.dot(D.dot(F)))
    C_cal = ((A + B.dot(F)).conj().T).dot(Y2.dot(A + B.dot(F))) - Y1
    D_cal = U_.dot((B_cal + C_cal).dot(U))

    return D_cal


def get_Upsilon(As, Bs, Fs, Theta, Upsilons, i, k):
    """Calculate current Upsilon"""
    A = As[Theta[i, k]]
    B = Bs[Theta[i, k]]
    F = Fs[Theta[i, k]]

    return ((A + B.dot(F))).dot(Upsilons[i])


def get_F(As, Bs, Ds, Fs, Ys, N):
    """Calculate F."""
    for i in range(N):
        A = As[i]
        B = Bs[i]
        D = Ds[i]
        S = Ys[i]

        B_ = B.conj().T
        D_ = D.conj().T

        Fs[i] = (-inv(B_.dot(S).dot(B) + D_.dot(D))).dot(B_).dot(S).dot(A)

    return Fs


def mjlstd(lambda_par, J, T, K, epsilon, N, P, As, Bs, Cs,
           Ds, Ys_par, Fs_par, seed, c, eta):
    """Applies the TD(\lambda) method to solve a MJLS."""

    Ys = Ys_par.copy()
    Fs = Fs_par.copy()

    Theta = np.zeros((N, K), dtype=int)
    Theta[:, 0] = [i for i in range(N)]

    np.random.seed(seed)

    for j in range(J):
        # Log
        Ys_old = Ys.copy()
        Fs_old = Fs.copy()

        # Calculate updated Ys and Fs
        Ys = get_Y(As, Bs, Cs, Ds, Fs, P, Ys, Theta, N, T,
                   K, c, eta, lambda_par, epsilon)
        Fs = get_F(As, Bs, Ds, Fs, Ys, N)

        # Log
        err_Ys = abs(Ys_old - Ys).max()
        err_Fs = abs(Fs_old - Fs).max()

        err_Ys_par = abs(Ys_par - Ys).max()
        err_Fs_par = abs(Fs - Fs_par).max()

        print('err_Ys: %3.1e\terr_Fs: %3.1e\terr_Ys_par: %e\terr_Fs_par: %e' %
              (err_Ys, err_Fs, err_Ys_par, err_Fs_par))

    return (Fs, Ys)
