import math
import numpy as np
import scipy.linalg as la


def test_stability(lambda_par, F, A, B, P, N):
    """Tests stability of (A, B) according to (1).

(1) O. L. V. Costa and J. C. C. Aya, "Monte Carlo
TD(\lambda_par)-methods for the optimal control of
discrete-time Markovian jump linear systems, "
Automatica, vol. 38, pp. 217–225, 2002."""

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

    print('max(abs(la.eig(I_cal)[0])): %s' %
          (max(abs(la.eig(I_cal)[0])),))

    if not (max(abs(la.eig(I_cal)[0])) < 1):
        print('WARNING: F does not satisfy Lemma 3 ' +
              '(stabilizability)')

    print('pow(lambda_par, 2): %s' % (pow(lambda_par, 2),))
    print('v_max: %s' % (v_max,))
    if not (pow(lambda_par, 2) < v_max):
        print('WARNING: lambda_par does not satisfy ' +
              ' Lemma 2 (convergence)')
