from . import wraps, np, aeseries, aerray, u
from AeViz.grid.grid import grid
from AeViz.utils.math_utils import function_average
import warnings
from AeViz.utils.files.string_utils import split_number_and_unit

## PROFILE GRID

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

## 2D PLOTS

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
        if (type(kwargs['plane']) == aerray or \
            isinstance(kwargs['plane'], (int, float)) or \
            (type(kwargs['plane']) == str and \
             kwargs['plane'] not in ['radius', 'r', 'theta', 'th', 'phi', 'ph',
                                     'mass', 'm',
                                     'xy', 'yx', 'xz', 'zx', 'zy', 'yz',
                                     'xz_phi_avg', 'zx_phi_avg', 'yz_phi_avg',
                                     'zy_phi_avg'])) and args[0].dim == 3:
            kwargs['plane'] = (kwargs['plane'], None, None)
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
                outdata.append(aeseries(_plane_cut(args[0].dim, dd, gr,
                                                   index_theta, index_phi),
                                                   **{X.name: X, Y.name: Y}))
        elif kwargs['plane'] in ['mass', 'm']:
            gmass = args[0].mass_grid(args[1])
            for dd in data:
                outdata.append(aeseries(function_average(dd,
                                                         args[0].dim, 'Omega',
                                            args[0].cell.dOmega(args[0].ghost)),
                               mass=gmass)) 
        elif kwargs['plane'] in ['xz_phi_avg', 'zx_phi_avg', 'yz_phi_avg',
                                 'zy_phi_avg']:
            gr = grid(args[0].dim, args[0].cell.radius(args[0].ghost),
                                       args[0].cell.theta(args[0].ghost),
                                       args[0].cell.phi(args[0].ghost))
            X, Y = gr.cartesian_grid_2D(kwargs['plane'], 64)
            for dd in data:
                outdata.append(aeseries(
                    function_average(dd, args[0].dim, 'only_phi', 
                                     args[0].cell.dphi(args[0].ghost)[:, None, None]),
                                     **{X.name: X, Y.name: Y}))
        else:
            raise TypeError(f"plane type {kwargs['plane']} not recognized")
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
    if type(sim).__name__ == 'AeViz':
        sim = sim.loaded_data
    if sim.dim == 1:
        index_phi = None
        index_theta = None
    elif sim.dim == 2:
        if plane.casefold() in ['xy', 'yx']:
            index_theta = len(sim.cell.theta(sim.ghost)) // 2
        else:
            index_theta = None
        index_phi = None
    elif plane.casefold() in ['xy', 'yx']:
        index_phi = None
        index_theta = len(sim.cell.theta(sim.ghost)) // 2
    elif plane.casefold() in ['xz', 'zx']:
        index_phi = len(sim.cell.phi(sim.ghost)) // 2
        index_theta = None
    elif plane.casefold() in ['yz', 'zy']:
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
        elif data.shape == radius.shape:
            return aeseries(data,
                            radius=radius)
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
    assert all([pl is None or type(pl) in [list, int, float] for pl in plane]), \
        "Indices must either be int, None, list, float or range."
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
        if len(plane) > 3:
            warnings.warn("Too many indices, considering only the first three.")
            plane = tuple(list(plane[:3]))
        if isinstance(plane[0], (aerray, float, str)) or \
            (len(plane) == 1 and isinstance(plane[0], int)):
            pl = plane[0]

            if type(pl) == str:
                pl = split_number_and_unit(pl)
            idx = np.argmax(sim.cell.radius(sim.ghost) >= pl)
            return aeseries(data[..., idx].T,
                            phi=sim.cell.phi(sim.ghost),
                            theta=sim.cell.theta(sim.ghost)-(np.pi/2 * u.radian))
        if plane[0] is None and plane[1] is None and plane[2] is None:
            raise TypeError("Three None are not supported.")
        elif plane[0] is None and plane[1] is None:
            if not type(plane[2]) == int:
                raise TypeError("Only int index allowed with two None.")
            return [aeseries(data[plane[2], i, :], radius=radius) for i in
                    range(data.shape[1])]
        elif plane[0] is None and plane[2] is None:
            if not type(plane[1]) == int:
                raise TypeError("Only int index allowed with two None.")
            return [aeseries(data[i, plane[1], :], radius=radius) for i in
                    range(data.shape[0])]
        elif plane[0] is None and plane[1] is None:
            if not type(plane[2]) == int:
                raise TypeError("Only int index allowed with two None.")
            return [aeseries(data[plane[2], i, :], radius=radius) for i in
                    range(data.shape[1])]
        elif plane[1] is None and plane[2] is None:
            if not type(plane[0]) == int:
                raise TypeError("Only int index allowed with two None.")
            return [aeseries(data[i, :, plane[0]], theta=theta) for i in
                    range(data.shape[0])]
        elif plane[0] is None:
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

## RADIUS IN 2D PLOTS

def get_radius(func):
    """
    Returns an aeseries with the 2D slice of the selected radius.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'plot' not in kwargs:
            return data
        if kwargs['rad'] != 'full':
            return data
        if isinstance(kwargs['time'], str):
            time = args[0].time(kwargs['time'])
        else:
            time = kwargs['time']
        radius_index = np.argmax(data[0].time >= time)
        rad = args[0].ghost.remove_ghost_cells_radii(data[0].data[..., radius_index],
                                                     args[0].dim, **data[1])
        if args[0].dim == 2:
            if kwargs['plane'].casefold() == 'xy':
                theta = np.linspace(0, 2 * np.pi, 64, True) * u.radian
                rad = rad[len(rad) // 2]
            else:
                theta = args[0].cell.theta(args[0].ghost)
            rad_x = np.sin(theta) * rad
            rad_y = np.cos(theta) * rad
        elif args[0].dim == 3:
            if kwargs['plane'].casefold() == 'xy':
                rad = rad[:, rad.shape[1] // 2]
                phi = args[0].cell.phi(args[0].ghost)
                rad_x = np.cos(phi) * rad
                rad_y = np.sin(phi) * rad
            elif kwargs['plane'].casefold() == 'xz':
                index_phi_1 = rad.shape[0] // 2
                index_phi_2 = (index_phi_1 + rad.shape[0] // 2) % rad.shape[0]
                theta = args[0].cell.theta(args[0].ghost)
                theta = np.concatenate((theta, theta+np.pi * u.radian))
                rad = np.concatenate((np.flip(rad[index_phi_1, :]),
                                      rad[index_phi_2, :]))
                rad_x = np.sin(theta) * rad
                rad_y = np.cos(theta) * rad
            elif kwargs['plane'].casefold() == 'yz':
                index_phi_1 = rad.shape[0] // 4
                index_phi_2 = (index_phi_1 + rad.shape[0] // 2) % rad.shape[0]
                theta = args[0].cell.theta(args[0].ghost)
                theta = np.concatenate((theta, theta+np.pi * u.radian))
                rad = np.concatenate((np.flip(rad[index_phi_1, :]),
                                      rad[index_phi_2, :]))
                rad_x = np.sin(theta) * rad
                rad_y = np.cos(theta) * rad
            elif kwargs['plane'].casefold() == 'xz_phi_avg':
                theta = args[0].cell.theta(args[0].ghost)
                rad = np.mean(rad, axis=0)
                rad_x = np.sin(theta) * rad
                rad_y = np.cos(theta) * rad
        return aeseries(
                rad_y, X=rad_x
            )
    return wrapper


        