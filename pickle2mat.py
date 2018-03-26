from subprocess import call
from stack import stack
import pickle
from scipy.io import savemat
from os.path import join as ospj


##################################################
# Opens pickle data and saves it as matlab data. #
##################################################

xp = 'fair'
dirin = ospj('./pkl', xp)
dircerob = '/home/rafaelbeirigo/cerob/MarkovianSimulator/xp'

# Define names and filenames
algorithms = ['mjlstdon-t']

repetitions = 100

# Open data from the files using the filenames
for repetition_ in range(1, repetitions + 1):
    repetition = '{:03d}'.format(repetition_)
    for algorithm in algorithms:
        pkl_filename = 'Fs_H_' + algorithm + '-' + repetition + '.pkl'
        pkl_path = ospj(dirin, algorithm, pkl_filename)
        Fs_H = pickle.load(open(pkl_path, 'rb'))

        dirout = ospj(dircerob, xp, algorithm, repetition)
        call(['mkdir', '-p', dirout])

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

            # Save the data in matlab format
            mat_filename = 'loop_loop_iteration_{:04d}.mat'.format(i)
            mat_path = ospj(dirout, mat_filename)
            savemat(mat_path, mdict={'Fvec_tilde': Fs_stack})

            i = i + 1
