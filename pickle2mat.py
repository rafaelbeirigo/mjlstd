from stack import stack
import pickle
from scipy.io import savemat
import numpy as np
from subprocess import call

##################################################
# Opens pickle data and saves it as matlab data. #
##################################################

# Define names and filenames
a = [['mjlstdon', 'Fs_el_H-001.pkl'],
     ['mjlstdoff', 'Fs_H-001.pkl']]

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

        # Save the data in matlab format using the names
        # loop_loop_iteration_0001.mat
        filename = '{:s}/loop_loop_iteration_{:04d}.mat'.format(b[0], i)
        f = np.array(stack(d))
        savemat(filename, mdict={'Fvec_tilde': f.reshape(len(f), 1)})
        i = i + 1

    call(['cp', '-r', b[0],
          '/home/rafaelbeirigo/cerob/MarkovianSimulator/'])
