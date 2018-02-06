from mjlstd import mjlstd, Parameters
import sam_constants as sc
import sam_parameters as sp
from test_stability import test_stability
from MJLS import MJLS

print('wait...')

seed = 0

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.Xs_ric, sc.Fs_ric)

p = Parameters(sp.L, sp.T, sp.K, sp.lambda_, sp.epsilon,
               sp.c, sp.eta, seed)

test_stability(m, p.lambda_)

(Fs, Ys) = mjlstd(p, m)

print(Ys)
print(Fs)
