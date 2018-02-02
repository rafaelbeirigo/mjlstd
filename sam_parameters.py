import numpy as np

x_0 = np.array([[1], [1]])
theta_0 = 1

R = 1                           # Number of repetitions of the experiment

L = 20                          # Max iterations for global Y convergence
T = int(1e6)                    # Max iterations for individual Y convergence
K = int(1e6)                    # Max iterations for calculating the sum

lambda_ = 1e-1
c = 1e-6
eta = 1

epsilon = 1e-6
