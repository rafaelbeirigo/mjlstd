import numpy as np
from stack import stack
import arm_constants as ac
from file import load
from riccati import riccati
from scipy.io import savemat
from arm_constants import P as P_hardcoded

# Pasta com coisas relacionadas à obtenção por
# contagem/maximização de verossimilhança (2017count)
dir = 'count/'

# P.pkl e P_tilde.pkl foram obtidas de
# cerob/MarkovianSimulator/riccati/1000_iter/loop_loop_iteration_0100.mat
# Git: origin/plot c88519f4d18a719f0fd1401096c22b7600762187
# P_hardcoded foi encontrada no código fonte (consigo ser mais específico?)
Ps = (('P', load(dir + 'P.pkl')),
      ('P_tilde', load(dir + 'P_tilde.pkl')),
      ('P_hardcoded', P_hardcoded))

for name, P in Ps:
    # Calcular F a partir de P no Python
    args = {
        'T': int(1e6),
        'N': ac.N,
        'A': ac.A,
        'B': ac.B,
        'C': ac.C,
        'D': ac.D,
        'R': P,
        'epsilon': 1e-6,
    }
    [F, _] = riccati(**args)

    # Zerar as devidas colunas. Isso é necessário por causa
    # da possibilidade de falhas e à estrutura das matrizes
    # do MJLS, que possibilitam divergência caso as devidas
    # colunas não sejam zeradas
    i = 0
    for f in F:
        if i >= 8:
            if i <= 15:
                f[:, 0] *= 0.
                f[:, 3] *= 0.
            else:
                f[:, 1] *= 0.
                f[:, 4] *= 0.
        i += 1

    F = stack(F)
    F = np.reshape(np.array(F), (len(F), 1))

    folder = '/home/rafaelbeirigo/cerob/MarkovianSimulator/Ps_evaluation/' + name
    filename = 'loop_loop_iteration_0001.mat'
    path = '{:s}/{:s}'.format(folder, filename)
    savemat(path, mdict={'Fvec_tilde': F})
