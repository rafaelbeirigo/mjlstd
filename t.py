from mjlstd import mjlstd, Parameters
import sam_constants as sc
from test_stability import test_stability
from MJLS import MJLS

print('wait...')

L = 20
T = int(1e6)
K = int(1e6)

lambda_ = 1e-1
c = 1e-6
eta = 1

epsilon = 1e-6

seed = 0

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.X_ric, sc.F_ric)

test_stability(m, lambda_)

p = Parameters(L, T, K, lambda_, epsilon, c, eta, seed)
(Fs, Ys) = mjlstd(p, m)

print(Ys)
print(Fs)
