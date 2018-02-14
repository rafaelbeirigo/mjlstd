import numpy as np

x_0 = np.array([[0]])
theta_0 = 0

R = 1

L = 1
T = 1
K = 10

lambda_ = 1
c = 1
eta = 1

epsilon = 1e-6

sum_D = 22 * np.eye(2)

sum_Ds = np.array([22. * np.eye(2) for _ in range(2)])

seed = 0

i = 0
