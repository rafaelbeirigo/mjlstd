import numpy as np
from sam_run_episode import sam_run_episode


def sam_eval_F(F, R, T, lambda_par, theta_0, x_0, P, A, B, N):
    thetas = np.zeros((R, T))
    xs = np.zeros(((R, T) + x_0.shape))

    for r in range(R):
        xs[r], thetas[r] = sam_run_episode(T, theta_0, x_0, P, A, B, F, r)

    return np.std(xs, 0), np.mean(xs, 0)
