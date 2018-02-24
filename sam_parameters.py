import numpy as np

x_0 = np.array([[1], [1]])
theta_0 = 1

R = 100                          # Number of repetitions of the experiment

L = 5                            # Max iterations for global Y convergence
T = 100                          # Max iterations for individual Y convergence
K = 10                           # Max iterations for calculating the sum

lambda_ = 1e-1
c = 1e-1
eta = 1e0

epsilon = 1e-3
