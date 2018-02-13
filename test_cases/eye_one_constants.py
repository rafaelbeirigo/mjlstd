import numpy as np

N = 1

# Number of rows in the state array (or vector)
n = 1

# Number of columns in the state array (`1', if it is a vector)
m = 1

# Transition probabilities
P = np.array([[1.]])

A = np.array([[[1.]]])
B = np.array([[[1.]]])
C = np.array([[[1.]]])
D = np.array([[[1.]]])

X = None

F = None

X_0 = np.array([[[1.]]])

F_0 = np.array([[[-0.5]]])

Upsilon_1 = np.array([[[0.5]]])
