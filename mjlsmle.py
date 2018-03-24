from file import saverep
from riccati import riccati
from sam_run_episode import get_next_theta
import arm_constants as ac
import arm_parameters as ap
import numpy as np


def normalize(M):
    row_sums = M.sum(axis=1)
    return M / row_sums[:, np.newaxis]


def flat_prior_fix(M):
    '''Fix any row in M that is not stochastic with the flat prior.'''
    rows = M.shape[0]
    for row in range(rows):
        if not np.isclose(M[row, :].sum(), 1.0):
            M[row, :] = np.ones(rows) / rows

    return M


def normalize_flatten(M):
    return flat_prior_fix(normalize(M))


def get_initial_state(N):
    theta = np.random.randint(0, N)
    if theta == 24:
        theta = 24
        pass
    return theta


def get_FX(P, T):
    args = {
        'T': T,
        'N': ac.N,
        'A': ac.A,
        'B': ac.B,
        'C': ac.C,
        'D': ac.D,
        'R': P,
        'epsilon': ap.epsilon,
    }

    return riccati(**args)


for r in range(ap.R):
    print('r: {:4d}'.format(r))
    np.random.seed(r)
    P_count = np.zeros((ac.N, ac.N))
    Fs_H, P_tilde_H = [], []
    for l in range(ap.L):
        print('l: {:4d}'.format(l))
        for t in range(ap.T):
            print('t: {:4d}'.format(t))
            theta = get_initial_state(ac.N)
            for _ in range(ap.K - 1):
                next_theta = get_next_theta(theta, ac.P)
                P_count[theta, next_theta] += 1
                theta = next_theta

        print('Normalizing (and flattening)...')
        P_tilde = normalize_flatten(P_count)

        print('Calculating F...')
        [Fs, _] = get_FX(P_tilde)

        print('Appending F and P to their histories...')
        Fs_H.append(Fs)
        P_tilde_H.append(P_tilde)

    saverep(Fs_H, 'Fs_H_mjlsmle', r)
    saverep(P_tilde_H, 'P_tilde_H_mjlsmle', r)
