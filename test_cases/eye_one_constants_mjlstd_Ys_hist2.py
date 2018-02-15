import numpy as np

N = 1

# Number of rows in the state array (or vector)
n = 1

# Number of columns in the state array (`1', if it is a vector)
m = 1

# Transition probabilities
P = np.array([[1.]])

A = np.array([[[0.5]]])
B = np.array([[[0.5]]])
C = np.array([[[1.]]])
D = np.array([[[1.]]])

X = np.array([[[1.]]])

F = np.array([[[1.]]])

mjlstd_F_Y = (np.array([[[-0.428571]]]),
              np.array([[[3.]]]))

# `get_Y_Ys_hist' tests
Ys_hist_len = 2
