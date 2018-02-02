from mjlstd import mjlstd, Parameters
import sam_constants as sc
import sam_constants as sp
from test_stability import test_stability
from MJLS import MJLS

print('wait...')

seed = 0

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.X_ric, sc.F_ric)

p = Parameters(sp.L, sp.T, sp.K, sp.lambda_, sp.epsilon,
               sp.c, sp.eta, sp.seed)

test_stability(m, p.lambda_)

(Fs, Ys) = mjlstd(p, m)

print(Ys)
print(Fs)
