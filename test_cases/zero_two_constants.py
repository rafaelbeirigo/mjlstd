import numpy as np

N = 2

# Number of rows in the state array (or vector)
n = 2

# Number of columns in the state array (`1', if it is a vector)
m = 2

# Transition probabilities
P = np.eye(2)

A = np.array([np.zeros((2, 2)) for i in range(2)])
B = np.array([np.zeros((2, 2)) for i in range(2)])
C = np.array([np.zeros((2, 2)) for i in range(2)])
D = np.array([np.zeros((2, 2)) for i in range(2)])

X = None

F = None

# `get_sum_D' tests
X_get_sum_D = np.array([np.zeros((2, 2)) for i in range(2)])
F_get_sum_D = np.array([np.zeros((2, 2)) for i in range(2)])

# `get_Y' tests
X_get_Y = np.array([np.zeros((2, 2)) for i in range(2)])
F_get_Y = np.array([np.zeros((2, 2)) for i in range(2)])
got_Y = np.array([np.zeros((2, 2)) for i in range(2)])

# `get_Y_Ys_hist' tests
Ys_hist_len = 1
