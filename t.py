from mjlstd import mjlstd
import sam_constants
import sam_parameters
from test_stability import test_stability
from MJLS import MJLS

print('wait...')

sc = sam_constants
sp = sam_parameters

J = 20
T = int(1e6)
K = int(1e6)

lambd = 1e-1
c = 1e-6
eta = 1

epsilon = 1e-6

seed = 0

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.X_ric, sc.F_ric)

test_stability(m, lambd)

(Fs, Ys) = mjlstd(m, lambd, J, T, K, epsilon, seed, c, eta)

print(Ys)
print(Fs)
