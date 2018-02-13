import numpy as np
from numpy import matmul


def get_next_theta(theta, P):
    if P.shape[0] != P.shape[1]:
        raise ValueError("P must be a square matrix.")

    if (P < 0).any():
        raise ValueError("P must contain only non-negative values.")

    row_sum = P.cumsum(1)[:, -1]
    ones = np.ones_like(row_sum)
    if not (np.allclose(row_sum, ones)):
        raise ValueError("Each row of P must sum to 1.")

    d = np.cumsum(P, 1)[theta]
    r = np.random.random()
    return np.nonzero(d > r)[0][0]


def sam_run_episode(T, theta_0, x_0, P, As, Bs, Fs, seed):
    np.random.seed(seed)

    thetas = np.zeros(T, dtype=int)
    xs = np.zeros((T,) + x_0.shape)

    thetas[0] = theta_0
    xs[0] = x_0

    for t in range(1, T):
        thetas[t] = get_next_theta(thetas[t-1], P)

        A = As[thetas[t-1]]
        B = Bs[thetas[t-1]]
        F = Fs[thetas[t-1]]
        xs[t] = matmul((A + matmul(B, F)), xs[t-1])

    return (xs, thetas)
