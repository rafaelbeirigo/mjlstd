from functools import reduce
from riccati import riccati
import arm_constants as ac
import arm_parameters as ap
from Parameters import Parameters
from MJLS import MJLS
from mjlstdoff import mjlstdoff
from mjlstdon import mjlstdon
from file import saverep, loadrep, load
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


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


def plot_Delta_Y_sum(Y_off_H, Y_el_H, X_ric, F_ric, fontsize=15):
    plt.figure()

    matplotlib.rc('xtick', labelsize=fontsize-2)
    matplotlib.rc('ytick', labelsize=fontsize-2)

    data = [(Y_off_H, 'Offline', 'red'),
            (Y_el_H, 'Online', 'blue')]
    for Y_H, label, color in data:
        Y_avg_orig = load('Y_avg-' + label + '.pkl')
        Y_std_orig = load('Y_std-' + label + '.pkl')

        Y_avg = Y_avg_orig[0::10]/max(Y_avg_orig)
        Y_std = Y_std_orig[0::10]/max(Y_avg_orig)

        x = range(len(Y_avg))
        plt.plot(x, Y_avg, label=label, color=color)
        plot_error(x, Y_avg, Y_std)

    plt.legend(loc=1, fontsize=fontsize)

    plt.xticks((0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000),
               ('0',
                '2500',
                '5000',
                '7500',
                '10000',
                '12500',
                '15000',
                '17500',
                '20000'), fontsize=10)

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
    return (x - y)


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


def main():
    """Runs the TD(\lambda) algorithm for the Samuelson problem."""
    print('wait for it...')

    args = {
        'L': ap.L,
        'T': ap.T,
        'K': ap.K,
        'lambda_': ap.lambda_,
        'epsilon': ap.epsilon,
        'c': ap.c,
        'eta': ap.eta,
        'seed': 0,
    }
    p = Parameters(**args)

    args = {
        'T': int(1e6),
        'N': ac.N,
        'A': ac.A,
        'B': ac.B,
        'C': ac.C,
        'D': ac.D,
        'R': ac.P,
        'epsilon': ap.epsilon,
    }
    [F_ric, X_ric] = riccati(**args)

    args = {
        'N': ac.N,
        'm': ac.m,
        'n': ac.n,
        'A': ac.A,
        'B': ac.B,
        'C': ac.C,
        'D': ac.D,
        'P': ac.P,
        'X': 0. * X_ric,
        'F': F_ric,
    }
    m = MJLS(**args)

    for r in range(ap.R):
        print('arm.py: Repetition {:3d} of {:3d} '
              '({:3.0f}%)'.format(r + 1, ap.R, 100. * (r + 1)/ap.R))

        Ys_H_mjlstdon = loadrep('Ys_H_mjlstdon', r)
        Fs_H_mjlstdon = loadrep('Fs_H_mjlstdon', r)

        if (Ys_H_mjlstdon is None or Fs_H_mjlstdon is None):
            print('Calculating...')
            p.seed = r

            (_, _, Fs_H_mjlstdon, Ys_H_mjlstdon) = mjlstdon(p, m)

            saverep(Ys_H_mjlstdon, 'Ys_H_mjlstdon', r)
            saverep(Fs_H_mjlstdon, 'Fs_H_mjlstdon', r)


if __name__ == '__main__':
    main()
