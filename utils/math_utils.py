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

from typing import Literal

def IDL_derivative(x, y, xvariable: Literal['radius', 'theta', 'phi'] = 'radius'):
    """
    Derivatie performed using three point Lagrangian interpolation, as in:
    `https://www.l3harrisgeospatial.com/docs/deriv.html` 
    """
    if not xvariable == 'radius':
        raise TypeError("Not implemented yet, sorry :)")
    assert x.shape == y.shape or x.shape[0] == y[..., :].shape[-1], \
                      "Arrays must have equal last dimension"
    while x.ndim != y.ndim:
        x = x[None, ...]
    assert x.shape[-1] >= 3, "To calculate this derivative you need AT LEAST three points."
    #first point
    x01 = x[..., 0] - x[..., 1]
    x02 = x[..., 0] - x[..., 2]
    x12 = x[..., 1] - x[..., 2]
    derivative = y[..., 0] * (x01 + x02) / \
        (x01 * x02) - y[..., 1] * x02 / (x01 * x12) \
        + y[..., 2] * x01 / (x02 * x12)
    #mid points
    x01 = x[..., : -2] - x[..., 1 : -1]
    x02 = x[..., : -2] - x[..., 2 :]
    x12 = x[..., 1 : -1] - x[..., 2 :]
    derivative = np.concatenate((derivative[..., None], 
                                 y[..., : -2] * x12 / (x01 * x02) + \
                                 y[..., 1 : -1] * (1. / x12 - 1. / x01) - \
                                 y[..., 2 :] * x01 / (x02 * x12)), axis = -1)
    #last point
    x01 = x[..., -3] - x[..., -2]
    x02 = x[..., -3] - x[..., -1]
    x12 = x[..., -2] - x[..., -1]
    
    return np.concatenate((derivative, (-y[...,-3] * x12 / (x01 * x02) +\
                                        y[..., -2] * x02 / (x01 * x12) - \
                                        y[..., -1] * (x02 + x12) / (x02 * x12))[..., None]),
                                        axis = -1)