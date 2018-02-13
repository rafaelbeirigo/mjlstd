import numpy as np

N = 2

# Number of rows in the state array (or vector)
n = 2

# Number of columns in the state array (`1', if it is a vector)
m = 2

# Transition probabilities
P = np.eye(2)

A = np.array([np.eye(2) for i in range(2)])
B = np.array([np.eye(2) for i in range(2)])
C = np.array([np.eye(2) for i in range(2)])
D = np.array([np.eye(2) for i in range(2)])

X = None

F = None

X_0 = np.array([np.eye(2) for i in range(2)])

F_0 = np.array([-0.5 * np.eye(2) for i in range(2)])

Upsilon_0 = np.array([[[1., 0],
                       [0, 1.]],
                      [[1., 0],
                       [0, 1.]],
                      [[1., 0],
                       [0, 1.]]])

Upsilon_1 = np.array([0.5 * np.eye(2) for i in range(2)])
