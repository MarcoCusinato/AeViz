import numpy as np



def function_average(qt, dim, av_type, dV):
    if dim == 1:
        indices = {'r': 0, 't': None, 'p': None}
    elif dim == 2:
        indices = {'r': 1, 't': 0, 'p': None}
    elif dim == 3:
        indices = {'r': 2, 't': 1, 'p': 0}

    if av_type == 'Omega':
        av = np.sum(qt * dV, axis=tuple(range(dim-1))) / np.sum(dV)
    elif av_type == 'theta':
        print(qt.shape, dV.shape)
        av = np.sum(qt * dV, axis=tuple([i for i in [indices['r'], indices['p'] ] if i is not None])) / np.sum(dV)
        print(av.shape)
    elif av_type == 'phi':
        av = np.sum(qt * dV, axis=tuple([i for i in [indices['r'], indices['t'] ] if i is not None])) / np.sum(dV)
    elif av_type == 'radius':
        av = np.sum(qt * dV, axis=tuple([i for i in [indices['t'], indices['p'] ] if i is not None])) / np.sum(dV)
    elif av_type == 'volume':
        av = np.sum(qt * dV) / np.sum(dV)
    return av