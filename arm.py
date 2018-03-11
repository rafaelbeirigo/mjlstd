from functools import reduce
from riccati import riccati
import arm_constants as sc
import arm_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd
from mjlstd_eligibility import mjlstd_eligibility
import matplotlib
import matplotlib.pyplot as plt
import pickle
import numpy as np
from subprocess import call


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


def get_E_cal_X(m, X):
    E_cal_X = np.zeros_like(X)
    for i in range(m.N):
        for j in range(m.N):
            E_cal_X[i] += m.P[i][j] * X[j]

    return E_cal_X


def plot_Y(m, Y_off, Y_el, X_ric, F_ric):
    plt.figure()

    plt.suptitle(r'Entries of Y at each $t,\ell$-step for '
                 'eligibility traces (blue), '
                 'offline (red), '
                 'vs E_cal(X) (black) ')

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

    E_cal_X = get_E_cal_X(m, X_ric)

    for p in plot:
        # Get the Ys indexes
        i, j, k = p[0][0], p[0][1], p[0][2]

        # Get the actual values to plot
        Y_off_plot = [y[i][j][k] for y in Y_off]
        Y_el_plot = [y[i][j][k] for y in Y_el]
        E_cal_X_plot = [E_cal_X[i][j][k] for y in Y_off]

        # Create the suplot
        plt.subplot(3, 3, p[1])
        plt.plot(Y_el_plot, 'blue')
        plt.plot(Y_off_plot, 'red')
        plt.plot(E_cal_X_plot, 'black')
        # Configure plot
        plt.ylabel(r'$Y_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel('(t,el)-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('Y_k_0_D_c_0.1_eligibility.pdf', bbox_inches='tight')
    plt.close()


def plot_F(m, F_off, F_el, X_ric, F_ric):
    plt.figure()

    plt.suptitle(r'Entries of F at each $\ell$-step for '
                 'eligibility traces (blue), '
                 'offline (red), '
                 'vs true optimal gain (black)')

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
        F_off_plot = [f[i][j][k] for f in F_off]
        F_el_plot = [f[i][j][k] for f in F_el]
        F_ric_plot = [F_ric[i][j][k] for f in F_off]

        # Create the suplot
        plt.subplot(3, 2, p[1])
        plt.step(range(len(F_el_plot)), F_el_plot, 'blue')
        plt.step(range(len(F_off_plot)), F_off_plot, 'red')
        plt.plot(F_ric_plot, 'black')
        # Configure plot
        plt.ylabel(r'$F_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel(r'$\ell$-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('F_k_0_D_c_0.1_eligibility.pdf', bbox_inches='tight')
    plt.close()


def plot_Delta(m, F_off, F_el, X_ric, F_ric):
    Delta_off = [100. * abs((F_ric - f)/F_ric) for f in F_off]
    Delta_el = [100. * abs((F_ric - f)/F_ric) for f in F_el]

    plt.figure()

    plt.suptitle(r'Entries of $\Delta$ for '
                 'eligibility traces (blue), '
                 'offline (red) '
                 'variants at each $\ell$-step')

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
        Delta_off_plot = [f[i][j][k] for f in Delta_off]
        Delta_el_plot = [f[i][j][k] for f in Delta_el]

        # Create the suplot
        plt.subplot(3, 2, p[1])
        plt.step(range(len(Delta_el_plot)), Delta_el_plot, 'blue')
        plt.step(range(len(Delta_off_plot)), Delta_off_plot, 'red')
        # Configure plot
        plt.ylabel(r'$\Delta_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel(r'$\ell$-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('Delta_k_0_D_c_0.1_eligibility.pdf', bbox_inches='tight')
    plt.close()


def plot_Delta_Y(m, Y_off, Y_el, X_ric, F_ric):
    Delta_off = [abs((X_ric - f)/X_ric) for f in Y_off]
    Delta_el = [abs((X_ric - f)/X_ric) for f in Y_el]

    plt.figure()

    plt.suptitle(r'Entries of $\Delta Y$ for '
                 'eligibility traces (blue), '
                 'offline (red) '
                 'variants at each $t,\ell$-step')

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
        Delta_off_plot = [f[i][j][k] for f in Delta_off]
        Delta_el_plot = [f[i][j][k] for f in Delta_el]

        # Create the suplot
        plt.subplot(3, 2, p[1])
        plt.step(range(len(Delta_el_plot)), Delta_el_plot, 'blue')
        plt.step(range(len(Delta_off_plot)), Delta_off_plot, 'red')
        # Configure plot
        plt.ylabel(r'$\Delta_{}({}, {})$'.format(i + 1, j + 1, k + 1))
        plt.xlabel(r'$\ell$-step')
        plt.grid(True)
        plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('Delta_Y_k_0_D_c_0.1_eligibility.pdf', bbox_inches='tight')
    plt.close()


def plot_Delta_Y_sum(m, Y_off_H, Y_el_H, X_ric, F_ric, fontsize=15):
    plt.figure()

    matplotlib.rc('xtick', labelsize=fontsize-2)
    matplotlib.rc('ytick', labelsize=fontsize-2)

    data = [(Y_off_H, 'Offline', 'red'),
            (Y_el_H, 'Online', 'blue')]
    for Y_H, label, color in data:
        F = A(Y_H, X_ric)
        Y_avg, Y_std = np.mean(F, 0), np.std(F, 0)

        x = range(len(Y_avg))
        plt.step(x, Y_avg, label=label, color=color)
        plot_error(x, Y_avg, Y_std)

    plt.legend(loc=1, fontsize=fontsize)

    # Configure plot
    plt.ylabel(r'$\Delta Y$', fontsize=fontsize)
    plt.xlabel(r'Episode', fontsize=fontsize)
    plt.grid(True)
    plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('arm_Y.pdf', bbox_inches='tight')

    plt.close()


def e(x):
    """Calculates the number of elements in x."""
    return reduce((lambda x, y: x * y), x.shape)


def d(x, y):
    """Calculates a normalized difference."""
    return (x - y) / y


def s(x):
    """Calculates the absolute sum of the elements."""
    return sum(sum(sum(abs(x))))


def D(x, y):
    """Calculates the Delta, normalized by the number of coefficients."""
    return s(d(x, y)) / e(x)


def A(x, y):
    """Calculates the Delta for an entire history."""
    Delta = []
    for seq in x:
        delta = []
        for f in seq:
            d = D(f, y)
            delta.append(d)
        Delta.append(delta)

    return np.array(Delta)


def plot_error(x, y, y_err):
    plt.fill_between(x,
                     [x - y for x, y in zip(y, y_err)],
                     [x + y for x, y in zip(y, y_err)],
                     step='pre',
                     color='silver')


def plot_Delta_F_sum(m, F_off_H, F_el_H, F_ric, fontsize=15):
    plt.figure()

    matplotlib.rc('xtick', labelsize=fontsize-2)
    matplotlib.rc('ytick', labelsize=fontsize-2)

    data = [(F_off_H, 'Offline', 'red'),
            (F_el_H, 'Online', 'blue')]
    for F_H, label, color in data:
        F = A(F_H, F_ric)
        F_avg, F_std = np.mean(F, 0), np.std(F, 0)

        x = range(len(F_avg))
        plt.step(x, F_avg, label=label, color=color)
        plot_error(x, F_avg, F_std)

    plt.legend(loc=1, fontsize=fontsize)

    # Configure plot
    plt.ylabel(r'$\Delta F$', fontsize=fontsize)
    plt.xlabel(r'Episode batch', fontsize=fontsize)
    plt.grid(True)
    plt.tight_layout(h_pad=0., w_pad=0., pad=2)

    plt.savefig('arm_F.pdf', bbox_inches='tight')

    plt.close()


def saverep(data, name, r):
    save(data, name + '-{:03d}.pkl'.format(r + 1))


def loadrep(name, r):
    return load(name + '-{:03d}.pkl'.format(r + 1))


def main():
    """Runs the TD(\lambda) algorithm for the Samuelson problem."""
    print('wait for it...')

    data = load('arm.pickle')
    if False and data is not None:
        (m, X_ric, F_ric, Ys_H_, Ys_el_H_, Fs_H_, Fs_el_H_) = data
    else:
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

        args = {
            'T': int(1e6),
            'N': sc.N,
            'A': sc.A,
            'B': sc.B,
            'C': sc.C,
            'D': sc.D,
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
            'X': 0. * X_ric,
            'F': F_ric,
        }
        m = MJLS(**args)

        Fs_H_, Ys_H_, Fs_el_H_, Ys_el_H_ = [], [], [], []
        for r in range(sp.R):
            print('arm.py: Repetition {:3d} of {:3d} '
                  '({:3.0f}%)'.format(r + 1, sp.R, 100. * (r + 1)/sp.R))

            Fs_H = loadrep('Fs_H', r)
            Ys_H = loadrep('Ys_H', r)
            Fs_el_H = loadrep('Fs_el_H', r)
            Ys_el_H = loadrep('Ys_el_H', r)

            if (Fs_H is None or Ys_H is None or Fs_el_H is None or Ys_el_H is None):
                p.seed = r

                (Fs, Ys, Fs_H, Ys_H) = mjlstd(p, m)
                (Fs_el, Ys_el, Fs_el_H, Ys_el_H) = mjlstd_eligibility(p, m)

                saverep(Fs_H, 'Fs_H', r)
                saverep(Ys_H, 'Ys_H', r)
                saverep(Fs_el_H, 'Fs_el_H', r)
                saverep(Ys_el_H, 'Ys_el_H', r)

            Fs_H_.append(Fs_H)
            Ys_H_.append(Ys_H)
            Fs_el_H_.append(Fs_el_H)
            Ys_el_H_.append(Ys_el_H)

    plot_Delta_Y_sum(m, Ys_H_, Ys_el_H_, X_ric, F_ric)


if __name__ == '__main__':
    main()
