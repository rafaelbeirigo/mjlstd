import numpy as np

N = 2

# Number of rows in the state array (or vector)
n = 2

# Number of columns in the state array (`1', if it is a vector)
m = 2

# Transition probabilities
P = np.eye(2)

A = np.array([0.5 * np.eye(2) for i in range(2)])
B = np.array([0.5 * np.eye(2) for i in range(2)])
C = np.array([np.eye(2) for i in range(2)])
D = np.array([np.eye(2) for i in range(2)])

X = np.array([np.eye(2) for i in range(2)])

F = np.array([np.eye(2) for i in range(2)])

# `get_sum_D' tests
X_get_sum_D = np.array([np.eye(2) for i in range(2)])
F_get_sum_D = np.array([np.eye(2) for i in range(2)])

# `get_Y' tests
X_get_Y = np.array([np.eye(2) for i in range(2)])
F_get_Y = np.array([np.eye(2) for i in range(2)])
got_Y = np.array([23. * np.eye(2) for i in range(2)])
