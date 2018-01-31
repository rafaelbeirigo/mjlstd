import math
import numpy as np
import scipy.linalg as la


def test_stability(lambda_par, F, A, B, P, N):
    krs = []
    v_max = -math.inf
    for i in range(N):
        G = A[i] + B[i].dot(F[i])
        kr = np.kron(G, G)
        krs.append(kr)

        # lambda_par test
        v = pow(la.norm(kr, ord=2), -1)
        if v > v_max:
            v_max = v

    # Fs test
    I_cal = la.block_diag(*krs)

    if not (max(abs(la.eig(I_cal)[0])) < 1):
        print('WARNING: F does not satisfy Lemma 3 ' +
              '(stabilizability)')

    if not (pow(lambda_par, 2) < v_max):
        print('WARNING: lambda_par does not satisfy ' +
              ' Lemma 2 (convergence)')