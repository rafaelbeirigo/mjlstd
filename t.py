import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd

print('wait for it...')

seed = 0

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

p = Parameters(sp.L, sp.T, sp.K, sp.lambda_, sp.epsilon,
               sp.c, sp.eta, seed)

test_stability(m, p.lambda_)

(Fs, Ys) = mjlstd(p, m)

print(Ys)
print(Fs)
