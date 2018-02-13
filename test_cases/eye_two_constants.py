import numpy as np

N = 2

# Number of rows in the state array (or vector)
n = 2

# Number of columns in the state array (`1', if it is a vector)
m = 2

# Transition probabilities
P = np.eye(N)

A = np.array([np.eye(n) for i in range(N)])
B = np.array([np.eye(n) for i in range(N)])
C = np.array([np.eye(n) for i in range(N)])
D = np.array([np.eye(n) for i in range(N)])

X = None

F = None
