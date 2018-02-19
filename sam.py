from riccati import riccati
import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd
from mjlstd_online import mjlstd_online
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


def plot_Y_H(m, Ys_H, Ys_on_H, X_ric, F_ric, factor):
    plt.figure()

    plt.suptitle('Entries of Y at each t,el-step for online (blue) '
                 'and offline (red) vs true Riccati solution (black)')

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
        Ys_on_plot = [y[i][j][k] for y in Ys_on_H]
        X_plot = [X_ric[i][j][k] for y in Ys_H]

        # Create the suplot
        plt.subplot(3, 3, p[1])
        plt.plot(Ys_on_plot, 'blue')
        plt.plot(Ys_plot, 'red')
        plt.plot(X_plot, 'black')
        # Configure plot
        plt.ylabel(r'$Y_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel('(t,el)-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('Y_k_0_D_{:06.2f}_c_0.1_online.png'.format(factor),
                bbox_inches='tight')
    plt.show()
    plt.close()


def plot_F_H(m, Fs_H, Fs_on_H, X_ric, F_ric, factor):
    F_H, F_on_H = Fs_H, Fs_on_H

    plt.figure()

    plt.suptitle('Entries of F at each el-step for online (blue) '
                 'and offline (red) vs true optimal gain (black)')

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
        F_on_plot = [f[i][j][k] for f in F_on_H]
        F_ric_plot = [F_ric[i][j][k] for f in F_H]

        # Create the suplot
        plt.subplot(3, 2, p[1])
        plt.step(range(len(F_on_plot)), F_on_plot, 'blue')
        plt.step(range(len(F_plot)), F_plot, 'red')
        plt.plot(F_ric_plot, 'black')
        # Configure plot
        plt.ylabel(r'$F_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel('el-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('F_k_0_D_{:06.2f}_c_0.1_online.png'.format(factor),
                bbox_inches='tight')
    plt.show()
    plt.close()


def plot_Delta_H(m, Fs_H, Fs_on_H, X_ric, F_ric, factor):
    F_H, F_on_H = Fs_H, Fs_on_H
    Delta_H = [100. * abs((F_ric - f)/F_ric) for f in F_H]
    Delta_on_H = [100. * abs((F_ric - f)/F_ric) for f in F_on_H]

    plt.figure()

    plt.suptitle(r'Entries of $\Delta$ for online (blue) and offline (red) '
                 ' variants at each el-step')

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
        # Get the Delta indexes
        i, j, k = p[0][0], p[0][1], p[0][2]

        # Get the actual values to plot
        Delta_plot = [f[i][j][k] for f in Delta_H]
        Delta_on_plot = [f[i][j][k] for f in Delta_on_H]

        # Create the suplot
        plt.subplot(3, 2, p[1])
        plt.step(range(len(Delta_on_plot)), Delta_on_plot, 'blue')
        plt.step(range(len(Delta_plot)), Delta_plot, 'red')
        # Configure plot
        plt.ylabel(r'$\Delta_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel('el-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('Delta_k_0_D_{:06.2f}_c_0.1_online.png'.format(factor),
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

    factors = [1.]
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

            (Fs, Ys, Fs_H, Ys_H) = mjlstd(p, m)
            (Fs_on, Ys_on, Fs_on_H, Ys_on_H) = mjlstd_online(p, m)

            data = (m, F_ric, X_ric, Fs, Ys, Fs_H, Ys_H, Fs_on, Ys_on,
                    Fs_on_H, Ys_on_H)
            # save(data, filename)

        (m, F_ric, X_ric, Fs, Ys, Fs_H, Ys_H, Fs_on, Ys_on, Fs_on_H,
         Ys_on_H) = data

        plot_Y_H(m, Ys_H, Ys_on_H, X_ric, F_ric, factor)
        plot_F_H(m, Fs_H, Fs_on_H, X_ric, F_ric, factor)
        plot_Delta_H(m, Fs_H, Fs_on_H, X_ric, F_ric, factor)


if __name__ == '__main__':
    main()
