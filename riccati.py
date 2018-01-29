import numpy as np
from numpy import matmul
from scipy.linalg import inv


def riccati(T, N, A, B, C, D, R):
    """Solve a MJLS by Riccati."""
    X = A * 0
    EX = X.copy()
    F = np.array([np.zeros([D[0].shape[1], A[0].shape[1]])
                  for i in range(N)])

    # Aliases for transpose conjugate
    A_ = np.array([A[i].conj().T for i in range(N)])
    B_ = np.array([B[i].conj().T for i in range(N)])
    C_ = np.array([C[i].conj().T for i in range(N)])
    D_ = np.array([D[i].conj().T for i in range(N)])
    for t in range(T):
        # Estimate X (EX)
        for i in range(N):
            EX[i] *= 0
            for j in range(N):
                EX[i] += R[i, j] * X[j]

        # Calculate F and X based on the estimate of X
        for i in range(N):
            # Aliases (for legibility purposes)
            F0 = matmul(B_[i], matmul(EX[i], B[i]))
            F0 = -inv(F0 + matmul(D_[i], D[i]))
            F[i] = matmul(F0, matmul(B_[i], matmul(EX[i], A[i])))

            X0 = matmul(A_[i], matmul(EX[i], A[i])) + matmul(C_[i], C[i])
            X[i] = X0 + matmul(A_[i], matmul(EX[i], matmul(B[i], F[i])))

    return [F, X]
