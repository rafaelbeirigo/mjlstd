from file import save
from scipy.io import loadmat


dir = 'count/'

P = loadmat(dir + 'P.mat')['P']
P_tilde = loadmat(dir + 'P_tilde.mat')['P_tilde']
Fvec_tilde = loadmat(dir + 'Fvec_tilde.mat')['Fvec_tilde']

save(P, dir + 'P.pkl')
save(P_tilde, dir + 'P_tilde.pkl')
save(Fvec_tilde, dir + 'Fvec_tilde.pkl')
