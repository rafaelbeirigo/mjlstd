import numpy as np

theta_0 = 0

R = 1

L = 1
T = 1
K = 0

lambda_ = 1.
c = 1.
eta = 1.

epsilon = 1e-6

sum_D = np.zeros((2, 2))

sum_Ds = np.array([np.zeros((2, 2)) for i in range(2)])

seed = 0

i = 0
