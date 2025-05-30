import numpy as np
from typing import Literal
from AeViz.units.aerray import aerray
from AeViz.utils.files.string_utils import merge_strings

def function_average(qt, dim, av_type:Literal['Omega', 'theta', 'phi',
                                              'only_phi', 'radius', 'volume'],
                     dV):
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
        if isinstance(qt, aerray):
            av.set(name=merge_strings(qt.name, '_ang_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_\Omega$'),
                   log=qt.log, limits=qt.limits)
    elif av_type == 'theta':
        if dim < 2:
            return qt
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['r'], 
                                                     indices['p']] 
                                         if i is not None])) / np.sum(dV)
        av.set(name=merge_strings(qt.name, '_th_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_\theta$'),
                   log=qt.log, limits=qt.limits)
    elif av_type == 'phi':
        if dim < 2:
            return qt
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['r'], 
                                                     indices['t']]
                                         if i is not None])) / np.sum(dV)
        av.set(name=merge_strings(qt.name, '_phi_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_\phi$'),
                   log=qt.log, limits=qt.limits)
    elif av_type == 'only_phi':
        if dim < 2:
            return qt
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['p']]
                                         if i is not None])) / np.sum(dV)
        av.set(name=merge_strings(qt.name, '_phi_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_\phi$'),
                   log=qt.log, limits=qt.limits, cmap=qt.cmap)
    elif av_type == 'radius':
        av = np.nansum(qt * dV, axis=tuple([i for i in [indices['t'], 
                                                     indices['p']] 
                                         if i is not None])) / np.sum(dV)
        av.set(name=merge_strings(qt.name, '_r_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_r$'),
                   log=qt.log, limits=qt.limits)
    elif av_type == 'volume':
        av = np.nansum(qt * dV) / np.sum(dV)
        av.set(name=merge_strings(qt.name, '_vol_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_V$'),
                   log=qt.log, limits=qt.limits)
    return av

def function_average_radii(qt, dim, dOmega):
    if dim == 1:
        return qt
    else:
        mask = np.isnan(qt)
        av =  np.nansum(qt * dOmega) / (dOmega[~mask]).sum()
        av.set(name=merge_strings(qt.name, '_ang_avg'), 
                   label=merge_strings(r'$\langle $', qt.label, 
                                        r'$\rangle_\Omega$'),
                   log=qt.log, limits=qt.limits)
        return av
    
def IDL_derivative(x, y, xvariable:Literal['radius', 'theta', 'phi']='radius',
                   axis=None):
    """
    Derivatie performed using three point Lagrangian interpolation, as in:
    `https://www.l3harrisgeospatial.com/docs/deriv.html` 
    """
    if axis is None:
        if xvariable == 'theta':
            y = np.moveaxis(y, -2, -1)
        elif xvariable == 'phi':
            y = np.moveaxis(y, -3, -1)
    else:
        y = np.moveaxis(y, axis, -1)
    
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
    
    if isinstance(derivative, aerray):
        name = merge_strings('d', y.name, '_d', x.name)
        if x.label == r'$t-t_\mathrm{b}$':
            label = merge_strings(r'$\partial_{t}$', y.label)
        else:
            label = merge_strings(r'$\partial_{$', x.label, r'$}$', y.label)
        derivative.set(name=name, label=label, cmap=y.cmap, log=y.log,
                       limits=[derivative.min().value, derivative.max().value])
        
    if axis is None:
        if xvariable == 'radius':
            return derivative
        elif xvariable == 'theta':
            return np.moveaxis(derivative, -1, -2)
        elif xvariable == 'phi':
            return np.moveaxis(derivative, -1, -3)
    else:
        return np.moveaxis(derivative, -1, axis)

def divergence(quantity, x, y, z,
               mode:Literal['spherical', 'cartesian'] = 'spherical'):
    """
    Calculate the divergence of a quantity in spherical or cartesian coordinates
    quantity: np.array of dim 3, z, y, x
              0 - x component
              1 - y component
              2 - z component
    x, x grid data or radius
    y, y grid data or theta
    z, z grid data or phi
    """
    assert quantity.ndim == 4, "The input array MUST be 4-d"
    assert quantity.shape[0] == 3, "The divergence can be computed only on a \
                                    3D array"
    assert quantity.shape[1] == z.shape[0], "Z-data do not match"
    assert quantity.shape[2] == y.shape[0], "Y-data do not match"
    assert quantity.shape[3] == x.shape[0], "X-data do not match"
    
    if mode == 'cartesian':
        return IDL_derivative(x, quantity[0, ...], 'radius') + \
               IDL_derivative(y, quantity[1, ...], 'theta') + \
               IDL_derivative(z, quantity[2, ...], 'phi')
    else:
        dr = 1 / x[None, None, :] ** 2 * IDL_derivative(x, quantity[0, ...] * 
                                                        x[None, None, :] ** 2)
        prefact = 1 / (x[None, None, :] * np.sin(y[None, :, None]))
        dth = prefact * \
            IDL_derivative(y, quantity[1] * np.sin(y[None, :, None]), 'theta')
        dph = prefact * \
            IDL_derivative(z, quantity[3], 'phi')
        return dr + dth + dph

def gradient(quantity, x, y, z,
             mode:Literal['spherical', 'cartesian'] = 'spherical'):
    """
    Calculate the gradient of a function
    Quantity: np array of dim z, y, x
    x, x grid data or radius
    y, y grid data or theta
    z, z grid data or phi
    Return np.array 3, z, y, x
    """
    
    if mode == 'cartesian':
        return np.concatenate((
            IDL_derivative(x, quantity, 'radius')[None, ...],
            IDL_derivative(y, quantity, 'theta')[None, ...],
            IDL_derivative(z, quantity, 'phi')[None, ...]
        ), axis=0)
    else:
        dr = IDL_derivative(x, quantity, 'radius')[None, ...]
        dth = (1 / x[None, None, :] * 
               IDL_derivative(y, quantity, 'theta'))[None, ...]
        dph = (1 / (x[None, None, :] * np.sin(y[None, :, None])) *
               IDL_derivative(z, quantity, 'phi'))[None, ...]
        return np.concatenate((dr, dth, dph), axis=0)

def get_stream_quantities(b1, b2, ax, ay, az, lx, ly, lz, plane):
    if plane != 'xz' and len(b1.shape) < 3:
        raise ValueError("Plane not found")
    if len(b1.shape) == 2:
        return b1, b2, ax, ay, az, lx, ly, lz
    if plane == 'xy':
        index = b1.shape[1] // 2
        b1, b2 = b1[:, index, :], b2[:, index, :]
        ax, ay, az = ax[:, index, :], ay[:, index, :], az[:, index, :]
        lx, ly, lz = lx[:, index, :], ly[:, index, :], lz[:, index, :]
    else:
        if plane == 'xz':
            index1 = b1.shape[0] // 2
            index2 = (index1 + b1.shape[0] // 2) % b1.shape[0]
        elif plane == 'yz':
            index1 = b1.shape[0] // 4
            index2 = (index1 + b1.shape[0] // 2) % b1.shape[0]
        b1 = [b1[index1, :, :], b1[index2, :, :]]
        b2 = [b2[index1, :, :], b2[index2, :, :]]
        ax = [ax[index1, :, :], ax[index2, :, :]]
        ay = [ay[index1, :, :], ay[index2, :, :]]
        az = [az[index1, :, :], az[index2, :, :]]
        lx = [lx[index1, :, :], lx[index2, :, :]]
        ly = [ly[index1, :, :], ly[index2, :, :]]
        lz = [lz[index1, :, :], lz[index2, :, :]]
    return b1, b2, ax, ay, az, lx, ly, lz

## The following is based on an Martin Obergaulinger's IDL script

"""
Stramlines calculation by integrating tangents to the vector field.
In this particular case we are interested in calculating the streamlines
of the magnetic field.
"""
def strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane):
    if isinstance(b1, aerray):
        b1 = b1.value
    if isinstance(b2, aerray):
        b2 = b2.value
    if isinstance(ax, aerray):
        ax = ax.value
    if isinstance(ay, aerray):
        ay = ay.value
    if isinstance(az, aerray):
        az = az.value
    if isinstance(lx, aerray):
        lx = lx.value
    if isinstance(ly, aerray):
        ly = ly.value
    if isinstance(lz, aerray):
        lz = lz.value
    if plane == 'yz':
        dF1 = np.zeros(b1.shape)
        dF2 = np.zeros(b1.shape)
        dl1 = np.ones(b1.shape)
        dl2 = 2 * np.ones(b1.shape)
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
    n1 = sb[0] + 1 #theta
    n2 = sb[1] + 1 #radius
    
    #stream fct
    F = np.zeros((n1, n2))
    dl1inv = 1. / dl1
    dl1inv = np.nan_to_num(dl1inv, nan = 0)
    dl2inv = 1. / dl2
    dl2inv = np.nan_to_num(dl2inv, nan = 0)

    #integrate the streamfunction
    for j in range(1, n1):
        jj = min(j, n1 - 2)
        F[j, 0] = (F[j, 0] * dl1[j - 1, 0] + dF2[j - 1, 0]) * dl2inv[jj, 0]
    for i in range(1, n2):
        ii = min(i, n2-2)
        F[0, i] = (F[0, i - 1] * dl2[0, i - 1] + dF1[0, i - 1]) * dl2inv[0, ii]
        for j in range(1, n1):
            jj = min(j, n1-2)
            F[j, i] = (F[j - 1, i] * dl1[j - 1, ii] + dF2[j - 1,ii]) \
                * dl1inv[jj, ii]

    F0 = F[0 : n1 - 1, 0 : n2 - 1] + F[1 : n1, 0 : n2 - 1] + \
        F [0 : n1 - 1, 1 : n2] + F[1 : n1, 1 : n2]

    if plane == 'yz':
        F0[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1]
        F[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1] * lx
    elif plane == 'xz':
        F0[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1]
        F[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1] * ly
    elif plane == 'xy':
        F0[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1]
        F[0 : n1 - 1, 0 : n2 - 1] = F[0 : n1 - 1, 0 : n2 - 1] * lz
    

    FF = F[0 : n1 - 1, 0 : n2 - 1] + F[1 : n1, 0 : n2 - 1] + \
        F[0 : n1 - 1, 1 : n2] + F[1 : n1, 1 : n2]
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
        b1 = b[..., 1] #theta b
        b2 = b[..., 2] #phi b
    elif plane == 'xz':
        b1 = b[..., 0] #rad b
        b2 = b[..., 2] #phi b
    elif plane == 'xy':
        b1 = b[..., 0] #rad b
        b2 = b[..., 1] #theta b

    b1, b2, ax, ay, az, lx, ly, lz = get_stream_quantities(b1, b2, ax, ay, az,
                                                           lx, ly, lz, plane)

        
    if type(b1) == list:
        F = np.concatenate((np.flip(
            strfunction2D(b1[0], b2[0], ax[0], ay[0], az[0], lx[0], ly[0],
                          lz[0], plane), axis=0),
            strfunction2D(b1[1], b2[1], ax[1], ay[1], az[1], lx[1], ly[1],
                          lz[1], plane)), axis=0)
    else:
        F = strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane)  

    return F