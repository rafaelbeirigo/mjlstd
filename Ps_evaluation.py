import numpy as np
from stack import stack
import arm_constants as ac
from file import load
from riccati import riccati
from scipy.io import savemat
# from arm_constants import P as P_hardcoded

dir = 'count/'

# P foi obtida de
# cerob/MarkovianSimulator/riccati/1000_iter/loop_loop_iteration_0100.mat
# Git: origin/plot c88519f4d18a719f0fd1401096c22b7600762187
P = load(dir + 'P.pkl')
# P_tilde = load(dir + 'P_tilde.pkl')
Fvec_tilde = load(dir + 'Fvec_tilde.pkl')

# O que acontece se eu
# 1. Calcular F a partir de P no Python
# 2. Stack essa F e salvar como .mat pelo Python
# 3. Abrir esse .mat no Matlab, recuperando a F
# 4. Plotar o controle a partir dessa F
# 5. Comparar com o controle utilizando a F do Matlab (a que
#    está salva junto com aP utilizada nesse experimento)

# 1. Calcular F a partir de P no Python
args = {
    'T': int(1e6),
    'N': ac.N,
    'A': ac.A,
    'B': ac.B,
    'C': ac.C,
    'D': ac.D,
    'R': ac.P,
    'epsilon': 1e-6,
}
[F, _] = riccati(**args)

i = 0
for f in F:
    # Zero the adequate columns
    if i >= 8:
        if i <= 15:
            f[:, 0] *= 0.
            f[:, 3] *= 0.
        else:
            f[:, 1] *= 0.
            f[:, 4] *= 0.
    i += 1

# 2. Stack essa F e salvar como .mat pelo Python
F = stack(F)
F = np.reshape(np.array(F), (len(F), 1))

folder = '/home/rafaelbeirigo/cerob/MarkovianSimulator/Ps_evaluation/P'
filename = 'loop_loop_iteration_0001.mat'
path = '{:s}/{:s}'.format(folder, filename)
savemat(path, mdict={'Fvec_tilde': F})
