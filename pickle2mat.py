from subprocess import call
from stack import stack
import pickle
from scipy.io import savemat
import numpy as np
import os


##################################################
# Opens pickle data and saves it as matlab data. #
##################################################

xp = 'P'
dirin = os.path.join('pkl', xp)
dircerob = '/home/rafaelbeirigo/cerob/MarkovianSimulator/'

# Define names and filenames
a = [['mjlstdon', os.path.join(dirin, 'Fs_el_H-001.pkl')],
     ['mjlstdoff', os.path.join(dirin, 'Fs_H-001.pkl')]]

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
        dirout = os.path.join(dircerob, b[0], xp)
        call(['mkdir', '-p', dirout])
        filename = 'loop_loop_iteration_{:04d}.mat'.format(i)
        path = os.path.join(dirout, filename)
        savemat(path, mdict={'Fvec_tilde': f})

        i = i + 1
