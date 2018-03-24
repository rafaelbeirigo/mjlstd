from subprocess import call
from stack import stack
import pickle
from scipy.io import savemat
import numpy as np
from os.path import join as ospj


##################################################
# Opens pickle data and saves it as matlab data. #
##################################################

xp = 'fair'
dirin = ospj('./pkl', xp)
dircerob = '/home/rafaelbeirigo/cerob/MarkovianSimulator/xp'

# Define names and filenames
algorithms = ['mjlsmle']

# Open the data from the files using the filenames
for algorithm in algorithms:
    pkl_filename = 'Fs_H_' + algorithm + '-001.pkl'
    pkl_path = ospj(dirin, pkl_filename)
    Fs_H = pickle.load(open(pkl_path, 'rb'))
    i = 1
    for Fs in Fs_H:
        j = 1
        for F in Fs:
            # Zero the adequate columns
            if j >= 9:
                if j <= 16:
                    F[:, 0] *= 0.
                    F[:, 3] *= 0.
                else:
                    F[:, 1] *= 0.
                    F[:, 4] *= 0.
            j += 1

        Fs_stack = stack(Fs)
        Fs_stack = np.reshape(np.array(Fs_stack),
                              (len(Fs_stack), 1))

        # Save the data in matlab format
        dirout = ospj(dircerob, algorithm, xp)
        call(['mkdir', '-p', dirout])
        mat_filename = 'loop_loop_iteration_{:04d}.mat'.format(i)
        mat_path = ospj(dirout, mat_filename)
        savemat(mat_path, mdict={'Fvec_tilde': Fs_stack})

        i = i + 1
