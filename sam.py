from riccati import riccati
import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd, get_F
from numpy import zeros_like, arange
import matplotlib.pyplot as plt


def main():
    """Runs the TD(\lambda) algorithm for the Samuelson problem."""
    print('wait for it...')

    seed = 0

    args = {
        'L': sp.L,
        'T': sp.T,
        'K': sp.K,
        'lambda_': sp.lambda_,
        'epsilon': sp.epsilon,
        'c': sp.c,
        'eta': sp.eta,
        'seed': seed,
    }
    p = Parameters(**args)

    factors = []
    factors.extend(arange(0.1, 2.1, 0.1))
    factors.extend(arange(3., 11., 1.))
    factors.extend(arange(20., 110., 10.))
    for factor in factors:
        args = {
            'T': int(1e6),
            'N': sc.N,
            'A': sc.A,
            'B': sc.B,
            'C': sc.C,
            'D': factor * sc.D,
            'R': sc.P,
            'epsilon': sp.epsilon,
        }
        [F_ric, X_ric] = riccati(**args)

        args = {
            'N': sc.N,
            'm': sc.m,
            'n': sc.n,
            'A': sc.A,
            'B': sc.B,
            'C': sc.C,
            'D': sc.D,
            'P': sc.P,
            'X': 0. * sc.X,
            'F': F_ric,
        }
        m = MJLS(**args)

        (Fs, Ys, Ys_H) = mjlstd(p, m)

        Fs_H = []
        for y in Ys_H:
            f = get_F(m, zeros_like(m.F), y).copy()
            Fs_H.append(f)

        Ys_H_error = [(y - X_ric).flatten() for y in Ys_H]
        Fs_H_error = [(f - F_ric).flatten() for f in Fs_H]

        plt.figure()

        plt.subplot(211)
        plt.title('k starting in 0; D * {:05.2f}; c = 0.1'.format(factor))

        plt.plot(Ys_H_error)
        plt.grid(True)
        plt.ylabel("Y(t) - X_ric")

        plt.subplot(212)
        plt.plot(Fs_H_error)
        plt.grid(True)
        plt.ylabel("F(t) - F_ric")
        plt.xlabel("t")

        plt.savefig('k_0_D_{:06.2f}_c_0.1.png'.format(factor),
                    bbox_inches='tight')
