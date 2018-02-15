import sam_constants as sc
import sam_parameters as sp
from Parameters import Parameters
from MJLS import MJLS
from mjlstd import mjlstd

print('wait for it...')

seed = 0

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.Xs_ric, sc.Fs_ric)

p = Parameters(sp.L, sp.T, sp.K, sp.lambda_, sp.epsilon,
               sp.c, sp.eta, seed)

test_stability(m, p.lambda_)

(Fs, Ys) = mjlstd(p, m)

print(Ys)
print(Fs)
