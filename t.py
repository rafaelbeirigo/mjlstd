from mjlstd import mjlstd
import sam_constants
import sam_parameters
from test_stability import test_stability

print('wait...')

sc = sam_constants
sp = sam_parameters

J = 20
T = int(1e5)
K = int(1e6)
epsilon = 1e-6

Fs = sc.F_sub
Xs = sc.X_riccati

test_stability(sp.lambda_par, Fs, sc.As, sc.Bs, sc.P, sc.N)

(Fs, Ys) = mjlstd(sp.lambda_par, J, T, K, sp.epsilon, sc.N, sc.P, sc.As,
                  sc.Bs, sc.Cs, sc.Ds, Xs, Fs, 0, sp.c, sp.eta)

print(Fs)
print(Ys)
