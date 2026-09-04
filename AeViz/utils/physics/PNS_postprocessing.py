from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
import h5py
from AeViz.utils.physics.radii_utils import PNS_radius
from AeViz.utils.math_utils import function_average_radii
from AeViz.units.constants import constants as c
from AeViz.utils.utils import check_existence, units_from_string
from AeViz.utils.files.file_utils import save_hdf

def declare_PNS_dictionary(dim: int,
                           magdim: int) -> dict:
    """
    _summary_

    Parameters
    ----------
    dim : int
        _description_
    magdim : int
        _description_

    Returns
    -------
    dict
        _description_
    """
    PNS_dictionary = {
        'global': {
            'mass': [],
            'rmax': [],
            'rmin': [],
            'ravg': [],
            'er': [],
            'eg': [],
            'ei': [],
            'I' : [],
            'mflux': []
            },
        'local': {
            's': [],
            'ye': [],
            't': [],
            'p': [],
            'vr': [],
            'mflux': [],
            'm': []
        }
    }
    if dim > 1:
        PNS_dictionary['global']['et'] = []
        PNS_dictionary['global']['ep'] = []
        PNS_dictionary['global']['jx'] = []
        PNS_dictionary['global']['jy'] = []
        PNS_dictionary['global']['jz'] = []
        PNS_dictionary['global']['jtot'] = []
        PNS_dictionary['local']['vt'] = []
        PNS_dictionary['local']['vp'] = []
        if magdim > 0:
            PNS_dictionary['global']['ebr'] = []
            PNS_dictionary['global']['ebt'] = []
            PNS_dictionary['global']['ebp'] = []
            PNS_dictionary['local']['br'] = []
            PNS_dictionary['local']['bt'] = []
            PNS_dictionary['local']['bp'] = []
    if dim > 2:
        PNS_dictionary['global']['erot'] = []
        PNS_dictionary['global']['eres'] = []
        PNS_dictionary['local']['vrot']  = []
        PNS_dictionary['local']['omg'] = []
        if magdim > 0:
            PNS_dictionary['global']['ebpol'] = []
            PNS_dictionary['global']['ebtor'] = []
            PNS_dictionary['local']['btor'] = []
            PNS_dictionary['local']['bpol'] = []   
    return PNS_dictionary    

def compute_PNS_postprocessing(PNS_dictionary: dict,
                               simulation: Simulation,
                               file_name: str,
                               gcells: dict,
                               dV: aerray,
                               dOmega: aerray,
                               Lx: aerray | None,
                               Ly: aerray | None,
                               Lz: aerray | None,
                               I: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNS_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    gcells : dict
        _description_
    dV : aerray
        _description_
    dOmega : aerray
        _description_
    Lx : aerray | None
        _description_
    Ly : aerray | None
        _description_
    Lz : aerray | None
        _description_
    I : aerray
        _description_
    """
    ## compute the radius
    PNSr = PNS_radius(simulation, file_name)
    r = simulation.cell.radius(simulation.ghost)
    if simulation.dim == 1:
        PNSmask = (r <= PNSr)
        idx_pns = np.argmax(r >= PNSr)
    else:
        PNSmask = (r <= simulation.ghost.remove_ghost_cells_radii(PNSr,
                                                          simulation.dim,
                                                          **gcells)[..., None])
        while r.ndim <= PNSr.ndim:
            r = r[None, :]
        ixd_pns = np.argmax(r >= PNSr[..., None], axis=-1)
    dmass = (simulation.rho(file_name) * dV).to(u.M_sun)
    __compute_1D_global_local_quantities(PNS_dictionary,
                                         simulation,
                                         file_name,
                                         dV,
                                         dOmega,
                                         PNSmask,
                                         idx_pns,
                                         PNSr,
                                         I,
                                         dmass)
    if simulation.dim > 1:
        __compute_2D_global_local_quantities(PNS_dictionary,
                                            simulation,
                                            file_name,
                                            PNSmask,
                                            idx_pns,
                                            dV,
                                            Lx,
                                            Ly,
                                            Lz,
                                            dmass)
    if simulation.dim == 3:
        __compute_3D_global_local_quantities(PNS_dictionary,
                                            simulation,
                                            file_name,
                                            PNSmask,
                                            idx_pns,
                                            dmass,
                                            dV)

def __compute_1D_global_local_quantities(PNS_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       dV: aerray,
                                       dOmega: aerray,
                                       PNS_mask: np.array[bool],
                                       indx_pns: np.array[int],
                                       PNS_radius: aerray,
                                       I: aerray,
                                       dmass:aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNS_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    dV : aerray
        _description_
    dOmega : aerray
        _description_
    PNS_mask : np.array[bool]
        _description_
    indx_pns : np.array[int]
        _description_
    PNS_radius : aerray
        _description_
    I : aerray
        _description_
    dmass : aerray
        _description_
    """
    mflux = simulation.mass_flux(file_name)
    grav_en = simulation.gravitational_energy(file_name) * dV
    int_en = simulation.internal_energy(file_name) * dV
    vr = simulation.radial_velocity(file_name)
    s = simulation.entropy(file_name)
    ye = simulation.Ye(file_name)
    t = simulation.temperature(file_name)
    p = simulation.gas_pressure(file_name)
    
    ## add the global variables
    PNS_dictionary['global']['mass'].append(np.nansum(dmass[PNS_mask]))
    PNS_dictionary['global']['rmax'].append(np.nanmax(PNS_radius))
    PNS_dictionary['global']['rmin'].append(np.nanmin(PNS_radius))
    PNS_dictionary['global']['ravg'].append(function_average_radii(PNS_radius,
                                                                   simulation.dim,
                                                                   dOmega))
    PNS_dictionary['global']['er'].append(np.nansum(0.5 * dmass[PNS_mask] * 
                                                    vr[PNS_mask] ** 2))
    PNS_dictionary['global']['eg'].append(np.nansum(grav_en[PNS_mask]))
    PNS_dictionary['global']['ei'].append(np.nansum(int_en[PNS_mask]))
    PNS_dictionary['global']['I'].append(np.nansum(I[PNS_mask]))
    
    ## add variables at the surface
    PNS_dictionary['local']['r'].append(PNS_radius)
    if simulation.dim == 1:
        PNS_dictionary['local']['s'].append(s[indx_pns])
        PNS_dictionary['local']['ye'].append(ye[indx_pns])
        PNS_dictionary['local']['t'].append(t[indx_pns])
        PNS_dictionary['local']['p'].append(p[indx_pns])
        PNS_dictionary['local']['vr'].append(vr[indx_pns])
        PNS_dictionary['local']['mflux'].append(mflux[indx_pns])
        PNS_dictionary['local']['m'].append(np.nansum(dmass[indx_pns]))
        PNS_dictionary['global']['mflux'].append(np.nansum(mflux[indx_pns] *
                                                           dOmega))
    elif simulation.dim == 2:
        itheta = np.arange(p.shape[0])
        PNS_dictionary['local']['s'].append(s[itheta, indx_pns])
        PNS_dictionary['local']['ye'].append(ye[itheta, indx_pns])
        PNS_dictionary['local']['t'].append(t[itheta, indx_pns])
        PNS_dictionary['local']['p'].append(p[itheta, indx_pns])
        PNS_dictionary['local']['vr'].append(vr[itheta, indx_pns])
        PNS_dictionary['local']['mflux'].append(mflux[itheta, indx_pns])
        PNS_dictionary['local']['m'].append(np.nancumsum(dmass[indx_pns],
                                                         axis=-1)[itheta,
                                                            indx_pns])
        PNS_dictionary['global']['mflux'].append(np.nansum(mflux[itheta,
                                                                 indx_pns] *
                                                                   dOmega))
    elif simulation.dim == 3:
            itheta = np.arange(p.shape[1])
            iphi = np.arange(p.shape[0])
            PNS_dictionary['local']['s'].append(s[iphi, itheta, indx_pns])
            PNS_dictionary['local']['ye'].append(ye[iphi, itheta, indx_pns])
            PNS_dictionary['local']['t'].append(t[iphi, itheta, indx_pns])
            PNS_dictionary['local']['p'].append(p[iphi, itheta, indx_pns])
            PNS_dictionary['local']['vr'].append(vr[iphi, itheta, indx_pns])
            PNS_dictionary['local']['mflux'].append(mflux[iphi, itheta, 
                                                          indx_pns])
            PNS_dictionary['local']['m'].append(np.nancumsum(dmass[indx_pns],
                                                                axis=-1)[iphi,
                                                                         itheta,
                                                                        indx_pns])
            PNS_dictionary['global']['mflux'].append(np.nansum(mflux[iphi,
                                                                     itheta,
                                                                     indx_pns] *
                                                                       dOmega))
            
def __compute_2D_global_local_quantities(PNS_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       PNS_mask: np.array[bool],
                                       indx_pns: np.array[int],
                                       dV: aerray,
                                       Lx: aerray,
                                       Ly: aerray,
                                       Lz: aerray,
                                       dmass: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNS_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    PNS_mask : np.array[bool]
        _description_
    indx_pns : np.array[int]
        _description_
    dV : aerray
        _description_
    Lx : aerray
        _description_
    Ly : aerray
        _description_
    Lz : aerray
        _description_
    dmass : aerray
        _description_
    """
    vth = simulation.theta_velocity(file_name)
    vph = simulation.phi_velocity(file_name)
    jx = np.nansum(Lx[PNS_mask])
    jy = np.nansum(Ly[PNS_mask])
    jz = np.nansum(Lz[PNS_mask])
    jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
    PNS_dictionary['global']['et'].append(np.nansum(0.5 * dmass[PNS_mask] * 
                                                    vth[PNS_mask]))
    PNS_dictionary['global']['ep'].append(np.nansum(0.5 * dmass[PNS_mask] *
                                                    vph))
    PNS_dictionary['global']['jx'].append(jx)
    PNS_dictionary['global']['jy'].append(jy)
    PNS_dictionary['global']['jz'].append(jz)
    PNS_dictionary['global']['jtot'].append(jtot)
    if simulation.dim == 2:
        itheta = np.arange(vth.shape[0])
        PNS_dictionary['local']['vt'].append(vth[itheta, indx_pns])
        PNS_dictionary['local']['vp'].append(vph[itheta, indx_pns])
    elif simulation.dim == 3:
        iphi = np.arange(vth.shape[0])[:, None]
        itheta = np.arange(vth.shape[1])[None, :]
        PNS_dictionary['local']['vt'].append(vth[iphi, itheta, indx_pns])
        PNS_dictionary['local']['vp'].append(vph[iphi, itheta, indx_pns])
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        if simulation.dim == 2:
            PNS_dictionary['local']['br'].append(br[itheta, indx_pns])
            PNS_dictionary['local']['bt'].append(bt[itheta, indx_pns])
            PNS_dictionary['local']['bp'].append(bp[itheta, indx_pns])
        elif simulation.dim == 3:
            PNS_dictionary['local']['br'].append(br[iphi, itheta, indx_pns])
            PNS_dictionary['local']['bt'].append(bt[iphi, itheta, indx_pns])
            PNS_dictionary['local']['bp'].append(bp[iphi, itheta, indx_pns])
        PNS_dictionary['global']['ebr'].append(np.nansum(0.5 * br[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        PNS_dictionary['global']['ebt'].append(np.nansum(0.5 * bt[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        PNS_dictionary['global']['ebp'].append(np.nansum(0.5 * bp[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        
def __compute_3D_global_local_quantities(PNS_dictionary: dict,
                                         simulation: Simulation,
                                         file_name: str,
                                         PNS_mask: np.array[bool],
                                         indx_pns: np.array[int],
                                         dmass: aerray,
                                         dV: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNS_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    PNS_mask : np.array[bool]
        _description_
    indx_pns : np.array[int]
        _description_
    dmass : aerray
        _description_
    dV : aerray
        _description_
    """
    jtot = PNS_dictionary['global']['jtot'][-1].value
    nx = PNS_dictionary['global']['jx'][-1] / jtot
    ny = PNS_dictionary['global']['jy'][-1] / jtot
    nz = PNS_dictionary['global']['jz'][-1] / jtot
    n = np.array([nx, ny, nz])
    if any([np.isnan(ni) for ni in n]):
        n = np.array([0, 0., 1.])
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    iphi = np.arange(vr.shape[0])[:, None]
    itheta = np.arange(vr.shape[1])[None, :]
    if n[2] >= 0.995:
        vrot = vp
        vturb = vt ** 2 + vr ** 2
        omg = simulation.omega(file_name)
        mask = PNS_mask
    else:
        X = simulation.grid.cartesian_grid()
        X_perp = X - np.tensordot(n, X, axes=(0, 0))[None, ...] * \
            n[:, None, None, None]
        R = np.linalg.norm(X_perp, axis=0) + 1e-12
        mask = ((R > simulation.grid.radius[0].value) & (PNS_mask))
        e_R = X_perp / (R)
        e_rot = np.cross(n[:, None, None, None], e_R, axis=0)
        vel_unit = vr.unit
        vx, vy, vz = simulation.grid.spherical_to_cartesian(
            vr, vt, vp
        )
        v = np.array([vx.value, vy.value, vz.value])
        vrot = np.nansum(v * e_rot, axis=0)
        omg = (vrot / R) * vel_unit / simulation.grid.radius.unit
        vturb = np.nansum((v - vrot * e_rot) ** 2, axis=0) * vel_unit ** 2
        vrot = vrot * vel_unit
    PNS_dictionary['global']['erot'].append(np.nansum(0.5 * (vrot ** 2 * 
                                                             dmass)[mask]))
    PNS_dictionary['global']['eres'].append(np.nansum(0.5 * (vturb *
                                                             dmass)[mask]))
    PNS_dictionary['local']['vrot'].append(vrot[iphi, itheta, indx_pns])
    PNS_dictionary['local']['omg'].append(omg[iphi, itheta, indx_pns])
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name)
        if n[2] >= 0.995:
            bpol = np.sqrt(br ** 2 + bt ** 2)
            btor = bp
        else:
            bunit = br.unit
            bx, by, bz = simulation.grid.spherical_to_cartesian(
            br, bt, bp
            )
            b = np.array([bx.value, by.value, bz.value])
            btor = np.nansum(b * e_rot, axis=0)
            bpol = np.sqrt(np.nansum((b - btor * e_rot) ** 2, axis=0)) * bunit
            btor = btor * bunit
        PNS_dictionary['global']['ebpol'].append(np.nansum(0.5 * (bpol ** 2 * 
                                                                     dV)[mask]) 
                                                 / c.mu0)
        PNS_dictionary['global']['ebtor'].append(np.nansum(0.5 * (btor ** 2 * 
                                                                    dV)[mask]) 
                                                         / c.mu0)
        PNS_dictionary['local']['bpol'].append(bpol[iphi, itheta, indx_pns])
        PNS_dictionary['local']['btor'].append(btor[iphi, itheta, indx_pns])