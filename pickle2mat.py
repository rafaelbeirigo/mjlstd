import pickle
from scipy.io import savemat

# Define names and filenames
a = [['mjlstdon', 'Fs_el_H-001.pkl'],
     ['mjlstdoff', 'Fs_H-001.pkl']]

# Open the data from the files using the filenames
for b in a:
    c = pickle.load(open(b[1], 'rb'))
    i = 1
    for d in c:
        # Save the data in matlab format using the names
        # loop_loop_iteration_0001.mat
        filename = '{:s}/loop_loop_iteration_{:04d}.mat'.format(b[0], i)
        savemat(filename, mdict={'Fvec_tilde': d})
        i = i + 1
