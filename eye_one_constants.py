import numpy as np

N = 1

# Number of rows in the state array (or vector)
n = 1

# Number of columns in the state array (`1', if it is a vector)
m = 1

# Transition probabilities
P = np.eye(N)

A = np.array([np.eye(n) for i in range(N)])
B = np.array([np.eye(n) for i in range(N)])
C = np.array([np.eye(n) for i in range(N)])
D = np.array([np.eye(n) for i in range(N)])

X = None

F = None
