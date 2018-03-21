from stack import stack
import pickle
from scipy.io import savemat
import numpy as np
import os


##################################################
# Opens pickle data and saves it as matlab data. #
##################################################

dir = 'count'

# Define names and filenames
a = [['mjlstdon', '{:s}/Fs_el_H-001.pkl'.format(dir)],
     ['mjlstdoff', '{:s}/Fs_H-001.pkl'.format(dir)]]

# Open the data from the files using the filenames
for b in a:
    c = pickle.load(open(b[1], 'rb'))
    i = 1
    for d in c:
        for e in d:
            # Zero the adequate columns
            if i >= 9:
                if i <= 16:
                    e[:, 0] *= 0.
                    e[:, 3] *= 0.
                else:
                    e[:, 1] *= 0.
                    e[:, 4] *= 0.

        f = stack(d)
        f = np.reshape(np.array(f), (len(f), 1))

        # Save the data in matlab format using the names
        # loop_loop_iteration_0001.mat
        dir_ = '/home/rafaelbeirigo/cerob/MarkovianSimulator/{:s}'.format(b[0])
        filename = 'loop_loop_iteration_{:04d}.mat'.format(i)
        path = os.path.join(dir_, filename)
        savemat(path, mdict={'Fvec_tilde': f})

        i = i + 1
