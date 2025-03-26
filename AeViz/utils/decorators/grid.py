from AeViz.units.aerray import aerray
from . import wraps, np, aeseries
from AeViz.grid.grid import grid
from AeViz.utils.math_utils import function_average
import warnings


def profile_grid(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'mesh' not in kwargs:
            kwargs['mesh'] = False
        if not kwargs['mesh']:
            return data
        T, R = np.meshgrid(data.time, data.radius)
        return aeseries(data.data, time=T, radius=R)
    return wrapper

def get_grid(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'plane' not in kwargs:
            return data
        if type(kwargs['plane']) == str:
            kwargs['plane'] = kwargs['plane'].casefold()
        if type(data) == tuple:
            data = list(data)
        if type(data) != list:
            data = [data]
        outdata = []
        if kwargs['plane'] in ['radius', 'r', 'theta', 'th', 'phi', 'ph']:
            for dd in data:
                outdata.append(_get_plane_avgs(args[0], dd, kwargs['plane']))
        elif type(kwargs['plane']) == tuple and len(kwargs['plane']) <= 3:
            for dd in data:
                outdata.append(_get_indices(args[0], dd, kwargs['plane']))
        elif kwargs['plane'] in ['xy', 'yx', 'xz', 'zx', 'zy', 'yz']:
            gr = grid(args[0].dim, args[0].cell.radius(args[0].ghost),
                                       args[0].cell.theta(args[0].ghost),
                                       args[0].cell.phi(args[0].ghost))
            X, Y = gr.cartesian_grid_2D(kwargs['plane'], 64)
            index_theta, index_phi = _get_plane_indices(args[0],
                                                        kwargs['plane'])
            for dd in data:
                outdata.append(aeseries(_plane_cut(args[0].dim, dd, gr, index_theta,
                                     index_phi), **{X.name: X, Y.name: Y})) 
        else:
            raise TypeError(f'plane type {kwargs['plane']} not recognized')
        if len(outdata) == 1:
            return outdata[0]
        else:
            return outdata
    return wrapper

def _get_plane_indices(sim, plane):
    """
    Get the indices associated to the plane. In particular, for 3D
    simulations, in the xz and yz planes, we need to extend the grid
    by taking the opposite slice.
    """
    ## Get the indices associated to the plane
    if sim.dim == 1:
        index_phi = None
        index_theta = None
    elif sim.dim == 2:
        if plane in ['xy', 'yx']:
            index_theta = len(sim.cell.theta(sim.ghost)) // 2
        else:
            index_theta = None
        index_phi = None
    elif plane in ['xy', 'yx']:
        index_phi = None
        index_theta = len(sim.cell.theta(sim.ghost)) // 2
    elif plane in ['xz', 'zx']:
        index_phi = len(sim.cell.phi(sim.ghost)) // 2
        index_theta = None
    elif plane in ['yz', 'zy']:
        index_theta = None
        index_phi = (len(sim.cell.phi(sim.ghost)) - (2 * \
            sim.ghost.ghost - sim.ghost.p_l - sim.ghost.p_r)) // 4 + \
                sim.ghost.ghost - sim.ghost.p_l
    return index_theta, index_phi

def _plane_cut(dim, data, gr, indextheta=None, indexphi=None):
    if dim == 1:
        return gr.map_1D_to_2D(data, 64)
    elif dim == 2:
        if (indexphi == None and indextheta == None):
            return data
        elif indextheta is not None:
            return gr.map_1D_to_2D(data, 64)
    elif dim == 3:
        if indexphi is not None:
            return np.concatenate([np.flip(data[indexphi, :, :], axis=0),
                                data[(indexphi + data.shape[0] // 2) % 
                                        data.shape[0], :, :]], axis=0)
        elif indextheta is not None:
            
            return data[:, indextheta, :]
    else:
        raise ValueError('The simulation dimension is not supported.')
    
def _get_plane_avgs(sim, data, plane):
    if plane in ['radius', 'r']:
        radius = sim.cell.radius(sim.ghost)
        if sim.dim == 1:
            return aeseries(data, radius=radius)
        else:
            return aeseries(function_average(data, sim.dim, 'Omega',
                                             sim.cell.dOmega(sim.ghost)),
                            radius=radius)
    elif plane in ['theta', 'th']:
        theta = sim.cell.theta(sim.ghost)
        if sim.dim == 1:
            raise TypeError("1D simulations do not have a theta plane.")
        elif sim.dim == 2:
            return aeseries(function_average(data, sim.dim, 'theta',
                                             sim.cell.dphi(sim.ghost) * \
                                             sim.cell.dr(sim.ghost)[None, :]),
                            theta=theta)
        elif sim.dim == 3:
            return aeseries(function_average(data, sim.dim, 'theta',
                                             sim.cell.dphi(sim.ghost)[:, None, None] * \
                                             sim.cell.dr(sim.ghost)[None, None, :]),
                            theta=theta)
    elif plane in ['phi', 'ph']:
        if sim.dim < 3:
            raise TypeError("1 and 2D simulations do not have phi dimension.")
        else:
            phi = sim.cell.phi(sim.ghost)
            return aeseries(function_average(data, sim.dim, 'phi',
                                             sim.cell.theta(sim.ghost)[None, :, None] * \
                                             sim.cell.dr(sim.ghost)[None, None, :]),
                            phi=phi)

def _get_indices(sim, data, plane):
    plane = tuple([list(pl) if type(pl) == range else pl for pl in plane])
    assert all([pl is None or type(pl) in [list, int] for pl in plane]), \
        "Indices must either be int, None, list or range."
    radius = sim.cell.radius(sim.ghost)
    theta = sim.cell.theta(sim.ghost)
    phi = sim.cell.phi(sim.ghost)
    if sim.dim == 1:
        if len(plane) > 1:
            warnings.warn("Too many indices, considering only the first.")
            plane = tuple(list(plane[0]))
        if plane is None:
            return aeseries(data, radius=radius)
        elif type(plane) == list:
            return [aeseries(data[i], radius=radius[i]) for i in plane]
        else:
            return aeseries(data[plane], radius[plane])        
    elif sim.dim == 2:
        if len(plane) > 2:
            warnings.warn("Too many indices, considering only the first two.")
            plane = tuple(list(plane[:2]))
        if plane[0] is None and plane[1] is None:
            return [aeseries(data[i, :], radius=radius) for i in
                    range(data.shape[0])]
        elif plane[0] is None:
            if type(plane[1]) == list:
                return [aeseries(data[i, :], radius=radius) for i in plane[1]]
            else:
                return aeseries(data[plane[1], :], radius=radius)
        elif plane[1] is None:
            if type(plane[0]) == list:
                return [aeseries(data[:, i], radius=theta) for i in plane[0]]
            else:
                return aeseries(data[:, plane[0]], radius=theta)
        else:
            raise TypeError("Not supported.")            
    else:
        if len(plane) > 2:
            warnings.warn("Too many indices, considering only the first three.")
            plane = tuple(list(plane[:3]))
        if plane[0] is None and plane[1] is None and plane[2] is None:
            raise TypeError("Three None are not supported.")
        elif plane[0] is None and plane[1] is None:
            if not type(plane[2]) == int:
                raise TypeError("Olnly int index allowed with two None.")
            return [aeseries(data[plane[2], i, :], radius=radius) for i in
                    range(data.shape[1])]
        elif plane[0] is None and plane[2] is None:
            if not type(plane[1]) == int:
                raise TypeError("Olnly int index allowed with two None.")
            return [aeseries(data[i, plane[1], :], radius=radius) for i in
                    range(data.shape[0])]
        elif plane[0] is None and plane[1] is None:
            if not type(plane[2]) == int:
                raise TypeError("Olnly int index allowed with two None.")
            return [aeseries(data[plane[2], i, :], radius=radius) for i in
                    range(data.shape[1])]
        elif plane[1] is None and plane[2] is None:
            if not type(plane[0]) == int:
                raise TypeError("Olnly int index allowed with two None.")
            return [aeseries(data[i, :, plane[0]], theta=theta) for i in
                    range(data.shape[0])]
        elif plane[0] is None:
            if not type(plane[1]) != type(plane[2]):
                raise TypeError("Must have at least a direction.")
            if type(plane[1]) == int:
                if type(plane[2]) == list:
                    return [aeseries(data[i, plane[1], :], radius=radius) 
                            for i in plane[2]]
                else:
                    return aeseries(data[plane[2], plane[1], :], radius=radius)
            else:
                return [aeseries(data[plane[2], i, :], radius=radius) 
                        for i in plane[1]]
        elif plane[1] is None:
            if not type(plane[0]) != type(plane[2]):
                raise TypeError("Must have at least a direction.")
            if type(plane[0]) == int:
                if type(plane[2]) == list:
                    return [aeseries(data[i, :, plane[0]], theta=theta) 
                            for i in plane[2]]
                else:
                    return aeseries(data[plane[2], :, plane[0]], theta=theta)
            else:
                return [aeseries(data[plane[2], :, i], theta=theta) 
                        for i in plane[0]]
        elif plane[2] is None:
            if not type(plane[1]) != type(plane[0]):
                raise TypeError("Must have at least a direction.")
            if type(plane[1]) == int:
                if type(plane[0]) == list:
                    return [aeseries(data[:, plane[1], i], phi=phi) 
                            for i in plane[0]]
                else:
                    return aeseries(data[:, plane[1], plane[0]], phi=phi)
            else:
                return [aeseries(data[:, i, plane[0]], phi=phi) 
                        for i in plane[1]]
        else:
            raise TypeError("Not supported")