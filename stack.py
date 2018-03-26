import numpy as np


def stack(a):
    '''Stacks a list of matrices in a 1-D list.'''
    c = []                      # shapes (cols and rows)
    d = []                      # data in row format
    count = 0
    for b in a:
        count = count + 1

        c.append(b.shape[0])
        c.append(b.shape[1])

        d = d + b.reshape(-1, order='F').tolist()

    e = [-1] + [count] + c + d

    f = np.reshape(np.array(e), (len(e), 1))

    return f
