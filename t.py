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

lambda_par = 1e-1
c = 1e-12
eta = 1

epsilon = 1e-6

m = MJLS(sc.As, sc.Bs, sc.Cs, sc.Ds,
         sc.P, sc.F_ric, sc.X_ric)

test_stability(lambda_par, Fs, sc.As, sc.Bs, sc.P, sc.N)

(Fs, Ys) = mjlstd(lambda_par, J, T, K, epsilon, sc.N, 0, c, eta)

print(Ys)
print(Fs)
