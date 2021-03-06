import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd, get_F
from numpy import zeros_like
import matplotlib.pyplot as plt

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

args = {
    'N': sc.N,
    'm': sc.m,
    'n': sc.n,
    'A': sc.A,
    'B': sc.B,
    'C': sc.C,
    'D': sc.D,
    'P': sc.P,
    'X': sc.X,
    'F': sc.F,
}
m = MJLS(**args)

(Fs, Ys, Ys_H) = mjlstd(p, m)

Fs_H = [get_F(m, zeros_like(m.F), y).copy() for y in Ys_H]

Ys_H_error = [(y - sc.X).flatten() for y in Ys_H]
Fs_H_error = [(f - sc.F).flatten() for f in Fs_H]

plt.figure(1)

plt.subplot(211)
plt.plot(Ys_H_error)
plt.grid(True)
plt.title("Initializing with \"Riccati\" data")
plt.ylabel("Y(t) - X_ric")
plt.xlabel("t")

plt.subplot(212)
plt.plot(Fs_H_error)
plt.grid(True)
plt.ylabel("F(t) - F_ric")
plt.xlabel("t")

plt.show()
