from mjlstd import mjlstd
import sam_constants
import sam_parameters
from test_stability import test_stability

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

Fs = sc.F_riccati
Xs = sc.X_riccati

test_stability(lambda_par, Fs, sc.As, sc.Bs, sc.P, sc.N)

(Fs, Ys) = mjlstd(lambda_par, J, T, K, epsilon, sc.N, sc.P, sc.As,
                  sc.Bs, sc.Cs, sc.Ds, Xs, Fs, 0, c, eta)

print(Ys)
print(Fs)
