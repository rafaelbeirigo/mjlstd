import riccati
import sam_constants as sc
import sam_parameters as sp

T = int(1e3)

(F, X) = riccati.riccati(T, sc.N, sc.As, sc.Bs, sc.Cs, sc.Ds, sc.P, sp.epsilon)

print(F)
print(X)
