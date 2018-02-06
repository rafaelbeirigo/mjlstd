import numpy as np

N = 2

# Transition probabilities
P = np.array([[0.9767, 0.0233],
              [0.0435, 0.9565]])

As = np.array([[0.8353], [0.9646]])

Bs = np.array([[0.0915], [0.0982]])

Cs = np.array([[np.sqrt(0.0355)],
               [np.sqrt(0.0355)]])

Ds = np.array([[1.0], [1.0]])

# # Calculated by Riccati
# Xs_ric = np.array([[0.1295], [0.3603]])

# # Calculated using Xs_ric above
# Fs_ric = np.array([[0.1295], [0.3603]])

# # Some known suboptimal control gains
# Fs_sub = np.array([[[2.4, -2.27340005]],
#                    [[43.65507178, -44.40141795]],
#                    [[-5.27993725, 6.05548422]]])
