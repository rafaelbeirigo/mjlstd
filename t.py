import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd

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

print(Ys)
print(Fs)
