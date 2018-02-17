from test_cases import acomp_F_trocado as sc
from test_cases import eye_two_parameters_acomp as sp
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
    'F': sc.F,
    'X': zeros_like(sc.X),
}
m = MJLS(**args)

(Fs, Ys, Ys_H) = mjlstd(p, m)

Fs_H = []
for y in Ys_H:
    f = get_F(m, zeros_like(m.F), y).copy()
    Fs_H.append(f)

Ys_H_error = [y.flatten() for y in Ys_H]
Fs_H_error = [f.flatten() for f in Fs_H]

plt.figure(1)

plt.subplot(211)
plt.plot(Ys_H_error)
plt.grid(True)
plt.title("Initializing with \"Riccati\" data")
plt.ylabel("Y(t) - X_ric")

plt.subplot(212)
plt.plot(Fs_H_error)
plt.grid(True)
plt.ylabel("F(t) - F_ric")
plt.xlabel("t")

plt.show()
