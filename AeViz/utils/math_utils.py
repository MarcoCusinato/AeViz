import numpy as np
from typing import Literal

def function_average(qt, dim, av_type:Literal['Omega', 'theta', 'phi',
                                              'radius', 'volume'], dV):
    if dim == 1:
        indices = {'r': 0, 't': None, 'p': None}
    elif dim == 2:
        indices = {'r': 1, 't': 0, 'p': None}
    elif dim == 3:
        indices = {'r': 2, 't': 1, 'p': 0}
    mask = np.isnan(qt)
    if av_type == 'Omega':
        if dim < 2:
            return qt
        if dV.ndim != qt.ndim:
            dV = dV[..., None]
        av = np.nansum(qt * dV, axis=tuple(range(dim-1))) / np.sum(dV)
    elif av_type == 'theta':
        if dim < 2:
            return qt
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['r'], 
                                                     indices['p']] 
                                         if i is not None])) / np.sum(dV)
    elif av_type == 'phi':
        if dim < 2:
            return qt
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['r'], 
                                                     indices['t']]
                                         if i is not None])) / np.sum(dV)
    elif av_type == 'radius':
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['t'], 
                                                     indices['p']] 
                                         if i is not None])) / np.sum(dV)
    elif av_type == 'volume':
        av = np.nansum(qt * dV) / np.sum(dV)
    return av

def function_average_radii(qt, dim, dOmega):
    if dim == 1:
        return qt
    else:
        mask = np.isnan(qt)
        return np.nansum(qt * dOmega) / (dOmega[~mask]).sum()
        

def IDL_derivative(x, y, xvariable:Literal['radius', 'theta', 'phi']='radius'):
    """
    Derivatie performed using three point Lagrangian interpolation, as in:
    `https://www.l3harrisgeospatial.com/docs/deriv.html` 
    """
    if xvariable == 'theta':
        y = np.moveaxis(y, -2, -1)
    elif xvariable == 'phi':
        y = np.moveaxis(y, -3, -1)
    
    assert x.shape == y.shape or x.shape[0] == y[..., :].shape[-1], \
                      "Arrays must have equal last dimension"
    while x.ndim != y.ndim:
        x = x[None, ...]
    assert x.shape[-1] >= 3, "To calculate this derivative you need AT LEAST"\
            " three points."
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
    
    derivative = np.concatenate((derivative, (-y[...,-3] * x12 / (x01 * x02) +\
                                        y[..., -2] * x02 / (x01 * x12) - \
                                        y[..., -1] * (x02 + x12) / \
                                            (x02 * x12))[..., None]), 
                            axis = -1)
    
    if xvariable == 'radius':
        return derivative
    elif xvariable == 'theta':
        return np.moveaxis(derivative, -1, -2)
    elif xvariable == 'phi':
        return np.moveaxis(derivative, -1, -3)
    
    return np.concatenate((derivative, (-y[...,-3] * x12 / (x01 * x02) +\
                                        y[..., -2] * x02 / (x01 * x12) - \
                                        y[..., -1] * (x02 + x12) / \
                                            (x02 * x12))[..., None]), 
                            axis = -1)

## The following is based on an Martin Obergaulinger's IDL script

"""
Stramlines calculation by integrating tangents to the vector field.
In this particular case we are interested in calculating the streamlines
of the magnetic field.
"""
def strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane):
    if plane == 'yz':
        dF1 = 0
        dF2 = 0
        dl1 = 1
        dl2 = 2
    elif plane == 'xz':
        dF1 =   b2 * az
        dF2 = - b1 * ax
        dl1 = ly
        dl2 = ly
    elif plane == 'xy':
        dF1 =   b2 * ay
        dF2 = - b1 * ax
        dl1 = lz
        dl2 = lz
    
    sb = b1.shape
    n1 = sb[0]+1 #theta
    n2 = sb[1]+1 #radius
    
    #stream fct
    F = np.zeros((n1, n2))
    dl1inv = 1./dl1
    dl1inv = np.nan_to_num(dl1inv, nan=0)
    dl2inv = 1./dl2
    dl2inv = np.nan_to_num(dl2inv, nan=0)

    #integrate the streamfunction
    for j in range(1, n1):
        jj = min(j, n1-2)
        F[j, 0] = (F[j, 0] * dl1[j-1, 0] + dF2[j-1, 0]) * dl2inv[jj, 0]
    for i in range(1, n2):
        ii = min(i, n2-2)
        F[0, i] = (F[0, i-1] * dl2[0, i-1] + dF1[0, i-1]) * dl2inv[0, ii]
        for j in range(1, n1):
            jj = min(j, n1-2)
            F[j, i] = (F[j-1, i] * dl1[j-1, ii] + dF2[j-1,ii]) * dl1inv[jj, ii]

    F0 = F[0:n1-1, 0:n2-1] + F[1:n1, 0:n2-1] + F [0:n1-1, 1:n2] + F[1:n1, 1:n2]

    if plane == 'yz':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]*lx
    elif plane=='xz':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]*ly
    elif plane == 'xy':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1 ] = F[0:n1-1, 0:n2-1]*lz
    

    FF = F[0:n1-1, 0:n2-1] + F[1:n1, 0:n2-1] + F [0:n1-1, 1:n2] + F[1:n1, 1:n2]
    FF =  FF / 4.
    F0 = F0 / 4.
    return FF

def strfct2D(b, cell, ghost, plane):
    """
    plane: 'yz', 'xz', 'xy'
    cell: object cell
    b magnetic field
    """
    assert plane in ['yz', 'xz', 'xy'], "plane must be one of: \'yz\', "\
        "\'xz\', \'xy\'"
    ax = cell.ax(ghost)
    ay = cell.ay(ghost)
    az = cell.az(ghost)
    lx = cell.lx(ghost)
    ly = cell.ly(ghost)
    lz = cell.lz(ghost)

    if plane == 'yz':
        b1 = b[:,:,1] #theta b
        b2 = b[:,:,2] #phi b
    elif plane == 'xz':
        b1 = b[:,:,0] #rad b
        b2 = b[:,:,2] #phi b
    elif plane == 'xy':
        b1 = b[:,:,0] #rad b
        b2 = b[:,:,1] #theta b

    F =  strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane)  

    return F