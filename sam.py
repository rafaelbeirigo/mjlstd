from riccati import riccati
import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd, get_F
from numpy import zeros_like, arange
import matplotlib.pyplot as plt
import pickle


def load(filename):
    """Loads persisted data."""
    try:
        with open(filename, 'rb') as f:
            return pickle.load(f)
    except Exception:
        return None


def save(data, filename):
    """Persists data."""
    with open(filename, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def plot_Y_H(m, Ys_H, X_ric, F_ric, factor):
    plt.figure()

    plt.suptitle('Entries of Y at each t,el-step (blue) '
                 'vs true Riccati solution (red)')

    # Pairs (("indexes on X"), ("index on the plot function"))
    plot = [
        ((0, 0, 0), (1)),
        ((0, 0, 1), (2)),
        ((0, 1, 1), (3)),
        ((1, 0, 0), (4)),
        ((1, 0, 1), (5)),
        ((1, 1, 1), (6)),
        ((2, 0, 0), (7)),
        ((2, 0, 1), (8)),
        ((2, 1, 1), (9)),
    ]

    for p in plot:
        # Get the Ys indexes
        i, j, k = p[0][0], p[0][1], p[0][2]

        # Get the actual values to plot
        Ys_plot = [y[i][j][k] for y in Ys_H]
        X_plot = [X_ric[i][j][k] for y in Ys_H]

        # Create the suplot
        plt.subplot(3, 3, p[1])
        plt.plot(Ys_plot, 'blue')
        plt.plot(X_plot, 'red')
        # Configure plot
        plt.ylabel(r'$Y_{}({}, {})$'.format(i +1, j + 1, k + 1))
        plt.xlabel('(t,el)-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0.,w_pad=0.,pad=2)

    plt.savefig('Y_k_0_D_{:06.2f}_c_0.1.png'.format(factor),
                bbox_inches='tight')
    plt.show()
    plt.close()


def plot_F_H(m, Ys_H, X_ric, F_ric, factor):
    F_H = [get_F(m, zeros_like(m.F), y) for y in Ys_H]

    plt.figure()

    plt.suptitle('Entries of F at each t,el-step (blue) '
                 'vs true optimal gain (red)')

    # Pairs (("indexes on X"), ("index on the plot function"))
    plot = [
        ((0, 0, 0), (1)),
        ((0, 0, 1), (2)),
        ((1, 0, 0), (3)),
        ((1, 0, 1), (4)),
        ((2, 0, 0), (5)),
        ((2, 0, 1), (6)),
    ]

    for p in plot:
        # Get the F indexes
        i, j, k = p[0][0], p[0][1], p[0][2]

        # Get the actual values to plot
        F_plot = [f[i][j][k] for f in F_H]
        F_ric_plot = [F_ric[i][j][k] for f in F_H]

        # Create the suplot
        plt.subplot(3, 3, p[1])
        plt.plot(F_plot, 'blue')
        plt.plot(F_ric_plot, 'red')
        # Configure plot
        plt.ylabel(r'$F_{}({}, {})$'.format(i +1, j + 1, k + 1))
        plt.xlabel('el-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0.,w_pad=0.,pad=2)

    plt.savefig('F_k_0_D_{:06.2f}_c_0.1.png'.format(factor),
                bbox_inches='tight')
    plt.show()
    plt.close()


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
        filename = 'k_0_D_{:06.2f}_c_0.1.pickle'.format(factor)
        data = load(filename)
        if data is None:
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

            data = (F_ric, X_ric, m, Fs, Ys, Ys_H)
            save(data, filename)

        (F_ric, X_ric, m, Fs, Ys, Ys_H) = data
        plot_Y_H(m, Ys_H, X_ric, F_ric, factor)
        plot_F_H(m, Ys_H, X_ric, F_ric, factor)


if __name__ == '__main__':
    main()
