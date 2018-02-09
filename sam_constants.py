import numpy as np

N = 3

# Number of rows in the state array (or vector)
n = 2

# Number of columns in the state array (`1', if it is a vector)
m = 1

# Transition probabilities
P = np.array([[0.67, 0.17, 0.16],
              [0.30, 0.47, 0.23],
              [0.26, 0.10, 0.64]])

As = np.array([[[0, 1], [-2.5, 3.2]],
               [[0, 1], [-43.7, 45.4]],
               [[0, 1], [5.3, -5.2]]])

Bs = np.array([np.array([[0], [1]]) for i in range(3)])

Cs = np.array([[[1.8974, -2.0028], [0, 0.9268], [0, 0]],
               [[3.1623, -0.9487], [0, 2.6646], [0, 0]],
               [[2.2361, -2.0125], [0, 0.6708], [0, 0]]])

Ds = np.array([[[0], [0], [np.sqrt(2.6)]],
               [[0], [0], [np.sqrt(1.165)]],
               [[0], [0], [np.sqrt(10./9)]]])

# Calculated by Riccati
Xs_ric = np.array([[[19.75510231, -18.57721303],
                    [-18.57721303, 48.283032]],
                   [[2232.50167333, -2263.49846231],
                    [-2263.49846231, 2342.81443288]],
                   [[36.09310699, -40.16022501],
                    [-40.16022501, 71.21455698]]])

# Calculated using Xs_ric above
Fs_ric = np.array([[[2.48538085, -2.27340005]],
                   [[43.65507178, -44.40141795]],
                   [[-5.27993725, 6.05548422]]])

# Some known suboptimal control gains
Fs_sub = np.array([[[2.4, -2.27340005]],
                   [[43.65507178, -44.40141795]],
                   [[-5.27993725, 6.05548422]]])
