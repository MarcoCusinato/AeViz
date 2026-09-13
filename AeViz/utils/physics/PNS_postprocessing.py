from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import PNS_radius
from AeViz.utils.math_utils import function_average_radii
from AeViz.units.constants import constants as c

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
            'mflux': [],
            'yeavg': [],
            'tavg': [],
            'savg': [],
            },
        'local': {
            'r': [],
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
        PNS_dictionary['local']['omg'] = []
        if magdim > 0:
            PNS_dictionary['global']['ebr'] = []
            PNS_dictionary['global']['ebt'] = []
            PNS_dictionary['global']['ebp'] = []
            PNS_dictionary['global']['bmax'] = []
            PNS_dictionary['global']['bavg'] = []
            PNS_dictionary['local']['br'] = []
            PNS_dictionary['local']['bt'] = []
            PNS_dictionary['local']['bp'] = []
            PNS_dictionary['local']['beta'] = []
    if dim > 2:
        PNS_dictionary['global']['erot'] = []
        PNS_dictionary['global']['eres'] = []
        PNS_dictionary['global']['cmass_x'] = []
        PNS_dictionary['global']['cmass_y'] = []
        PNS_dictionary['global']['cmass_z'] = []
        PNS_dictionary['local']['vrot']  = []
        if magdim > 0:
            PNS_dictionary['global']['ebpol'] = []
            PNS_dictionary['global']['ebtor'] = []
            PNS_dictionary['local']['btor'] = []
            PNS_dictionary['local']['bpol'] = []
            
    return PNS_dictionary    

def PNS_quantity_metadata() -> dict:
    """
    Metadata for PNS post-processing quantities.
    """

    meta = {
        "global": {
            "mass":   dict(name="PNS_mass",
                           label=r"$M_{\rm PNS}$",
                           limits=[1.0, 2.5],
                           cmap="viridis",
                           log=False),
            "rmax":   dict(name="PNS_radius_max",
                           label=r"$R_{\rm PNS, max}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "rmin":   dict(name="PNS_radius_min",
                           label=r"$R_{\rm PNS, min}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "ravg":   dict(name="PNS_radius_avg",
                           label=r"$\langle R_{\rm PNS}\rangle$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "mflux":  dict(name="PNS_mass_flux",
                           label=r"$\dot{M}_{\rm PNS}$",
                           limits=[1e-4, 5],
                           cmap="cividis",
                           log=True),
            "er":     dict(name="radial_kinetic_energy",
                           label=r"$E_{\rm PNS,kin,r}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "et":     dict(name="theta_kinetic_energy",
                           label=r"$E_{\rm PNS,kin,\theta}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "ep":     dict(name="phi_kinetic_energy",
                           label=r"$E_{\rm PNS,kin,\phi}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "eg":     dict(name="gravitational_energy",
                           label=r"$E_{\rm PNS,grav}$",
                           limits=[-5e53, -1e51],
                           cmap="magma",
                           log=False),
            "ei":     dict(name="internal_energy",
                           label=r"$E_{\rm PNS,int}$",
                           limits=[1e50, 5e53],
                           cmap="hot",
                           log=True),
            "I":      dict(name="moment_of_inertia",
                           label=r"$I_{\rm PNS}$",
                           limits=[1e43, 5e45],
                           cmap="Greens",
                           log=True),
            "yeavg":  dict(name="Ye_average",
                           label=r"$\langle Y_e\rangle_{\rm PNS}$",
                           limits=[0.0, 0.5],
                           cmap="RdYlBu",
                           log=False),
            "tavg":   dict(name="temperature_average",
                           label=r"$\langle T\rangle_{\rm PNS}$",
                           limits=[1, 60],
                           cmap="afmhot",
                           log=False),
            "savg":   dict(name="entropy_average",
                           label=r"$\langle s\rangle_{\rm PNS}$",
                           limits=[0, 30],
                           cmap="Spectral",
                           log=False),
            "jx":   dict(name="Jx",
                         label=r"$J_{x,\rm PNS}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jy":   dict(name="Jy",
                         label=r"$J_{y,\rm PNS}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jz":   dict(name="Jz",
                         label=r"$J_{z,\rm PNS}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jtot": dict(name="Jtot",
                         label=r"$|\mathbf{J}|_{\rm PNS}$",
                         limits=[1e46, 1e50],
                         cmap="viridis",
                         log=True),
            "ebr": dict(name="magnetic_energy_radial",
                        label=r"$E_{\rm mag,r,PNS}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "ebt": dict(name="magnetic_energy_theta",
                        label=r"$E_{\rm mag,\theta,PNS}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "ebp": dict(name="magnetic_energy_phi",
                        label=r"$E_{\rm mag,\phi,PNS}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "bmax": dict(name="maximum_magnetic_field",
                         label=r"$B_{\rm max,PNS}$",
                         limits=[1e10, 1e17],
                         cmap="PuRd",
                         log=True),
            "bavg": dict(name="average_magnetic_field",
                         label=r"$\langle B\rangle_{\rm PNS}$",
                         limits=[1e9, 1e16],
                         cmap="PuRd",
                         log=True),
            "erot": dict(name="rotational_energy",
                         label=r"$E_{\rm rot,PNS}$",
                         limits=[1e46, 1e53],
                         cmap="Oranges",
                         log=True),
            "eres": dict(name="turbulent_energy",
                         label=r"$E_{\rm turb,PNS}$",
                         limits=[1e46, 1e53],
                         cmap="YlOrBr",
                         log=True),
            "cmass_x": dict(name="center_of_mass_x",
                            label=r"$x_{\rm COM,PNS}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "cmass_y": dict(name="center_of_mass_y",
                            label=r"$y_{\rm COM,PNS}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "cmass_z": dict(name="center_of_mass_z",
                            label=r"$z_{\rm COM,PNS}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "ebpol": dict(name="poloidal_magnetic_energy",
                            label=r"$E_{\rm mag,pol,PNS}$",
                            limits=[1e42, 1e51],
                            cmap="PuBuGn",
                            log=True),
            "ebtor": dict(name="toroidal_magnetic_energy",
                            label=r"$E_{\rm mag,tor,PNS}$",
                            limits=[1e42, 1e51],
                            cmap="PuBuGn",
                            log=True),
        },

        "local": {
            "r":      dict(name="PNS_radius_surface",
                           label=r"$R_{\rm PNS}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "s":      dict(name="surface_entropy",
                           label=r"$s_{\rm PNS}$",
                           limits=[0, 40],
                           cmap="Spectral",
                           log=False),
            "ye":     dict(name="surface_Ye",
                           label=r"$Y_{e,\rm PNS}$",
                           limits=[0, 0.6],
                           cmap="RdYlBu",
                           log=False),
            "t":      dict(name="surface_temperature",
                           label=r"$T_{\rm PNS}$",
                           limits=[1, 60],
                           cmap="afmhot",
                           log=False),
            "p":      dict(name="surface_pressure",
                           label=r"$P_{\rm gas,PNS}$",
                           limits=[1e27, 1e35],
                           cmap="magma",
                           log=True),
            "vr":     dict(name="surface_radial_velocity",
                           label=r"$v_{r,\rm PNS}$",
                           limits=[-0.5, 0.5],
                           cmap="coolwarm",
                           log=False),
            "mflux":  dict(name="surface_mass_flux",
                           label=r"$\dot{M}_{\rm PNS}$",
                           limits=[1e-4, 5],
                           cmap="cividis",
                           log=True),
            "m":      dict(name="enclosed_mass",
                           label=r"$M(<R_{\rm PNS})$",
                           limits=[1.0, 2.5],
                           cmap="viridis",
                           log=False),
            "vt":  dict(name="surface_theta_velocity",
                        label=r"$v_{\theta,\rm PNS}$",
                        limits=[-0.5, 0.5],
                        cmap="coolwarm",
                        log=False),
            "vp":  dict(name="surface_phi_velocity",
                        label=r"$v_{\phi,\rm PNS}$",
                        limits=[-0.5, 0.5],
                        cmap="coolwarm",
                        log=False),
            "omg": dict(name="surface_omega",
                        label=r"$\Omega_{\rm PNS}$",
                        limits=[1, 1e4],
                        cmap="twilight",
                        log=True),
            "br": dict(name="surface_Br",
                       label=r"$B_{r,\rm PNS}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "bt": dict(name="surface_Btheta",
                       label=r"$B_{\theta,\rm PNS}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "bp": dict(name="surface_Bphi",
                       label=r"$B_{\phi,\rm PNS}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "beta": dict(name="surface_plasma_beta",
                         label=r"$\beta_{\rm PNS}$",
                         limits=[1e-2, 1e4],
                         cmap="viridis",
                         log=True),
            "vrot": dict(
                        name="surface_rotational_velocity",
                        label=r"$v_{\rm rot,PNS}$",
                        limits=[0, 0.5],
                        cmap="twilight_shifted",
                        log=False,
                        ),
            "bpol": dict(name="surface_Bpol",
                         label=r"$B_{\rm pol,PNS}$",
                         limits=[0, 1e16],
                         cmap="PuRd",
                         log=True),
            
            "btor": dict(name="surface_Btor",
                         label=r"$B_{\rm tor,PNS}$",
                         limits=[0, 1e16],
                         cmap="PuRd",
                         log=True),   
        }
    }

    return meta

def compute_PNS_postprocessing(PNS_dictionary: dict,
                               simulation: Simulation,
                               file_name: str,
                               dV: aerray,
                               dOmega: aerray,
                               Lx: aerray | None,
                               Ly: aerray | None,
                               Lz: aerray | None) -> None:
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
    Lx : aerray | None
        _description_
    Ly : aerray | None
        _description_
    Lz : aerray | None
        _description_
    """
    ## compute the radius
    PNSr = PNS_radius(simulation, file_name)
    r = simulation.cell.radius(simulation.ghost)
    if simulation.dim == 1:
        PNSmask = (r <= PNSr)
        idx_pns = np.argmax(r >= PNSr)
    else:
        while r.ndim <= PNSr.ndim:
            r = r[None, :]
        PNSmask = (r <= PNSr)
        idx_pns = np.argmax(r >= PNSr[..., None], axis=-1)
    dmass = (simulation.rho(file_name) * dV).to(u.M_sun)
    __compute_1D_global_local_quantities(PNS_dictionary,
                                         simulation,
                                         file_name,
                                         dV,
                                         dOmega,
                                         PNSmask,
                                         idx_pns,
                                         PNSr,
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

    mass = np.nansum(dmass[PNS_mask])
    ## add the global variables
    PNS_dictionary['global']['mass'].append(mass)
    PNS_dictionary['global']['rmax'].append(np.nanmax(PNS_radius))
    PNS_dictionary['global']['rmin'].append(np.nanmin(PNS_radius))
    PNS_dictionary['global']['ravg'].append(function_average_radii(PNS_radius,
                                                                   simulation.dim,
                                                                   dOmega))
    PNS_dictionary['global']['er'].append(np.nansum(0.5 * dmass[PNS_mask] * 
                                                    vr[PNS_mask] ** 2))
    PNS_dictionary['global']['eg'].append(np.nansum(grav_en[PNS_mask]))
    PNS_dictionary['global']['ei'].append(np.nansum(int_en[PNS_mask]))
    PNS_dictionary['global']['yeavg'].append(np.nansum(ye[PNS_mask] * dmass[PNS_mask]) /
                                             mass)
    PNS_dictionary['global']['savg'].append(np.nansum(s[PNS_mask] * dmass[PNS_mask]) /
                                                 mass)
    PNS_dictionary['global']['tavg'].append(np.nansum(t[PNS_mask] * dmass[PNS_mask]) /
                                                     mass)

    
    
    ## add variables at the surface
    PNS_dictionary['local']['r'].append(PNS_radius)
    if simulation.dim == 1:
        PNS_dictionary['local']['s'].append(s[indx_pns])
        PNS_dictionary['local']['ye'].append(ye[indx_pns])
        PNS_dictionary['local']['t'].append(t[indx_pns])
        PNS_dictionary['local']['p'].append(p[indx_pns])
        PNS_dictionary['local']['vr'].append(vr[indx_pns])
        PNS_dictionary['local']['mflux'].append(mflux[indx_pns])
        PNS_dictionary['local']['m'].append(mass)
        PNS_dictionary['global']['mflux'].append(np.nansum(mflux[indx_pns] *
                                                           dOmega))
        rperp2 = simulation.cell.radius(simulation.ghost) ** 2
        PNS_dictionary['global']['I'].append(np.nansum((dmass * 
                                                        rperp2)[PNS_mask]))
    elif simulation.dim == 2:
        itheta = np.arange(p.shape[0])
        PNS_dictionary['local']['s'].append(s[itheta, indx_pns])
        PNS_dictionary['local']['ye'].append(ye[itheta, indx_pns])
        PNS_dictionary['local']['t'].append(t[itheta, indx_pns])
        PNS_dictionary['local']['p'].append(p[itheta, indx_pns])
        PNS_dictionary['local']['vr'].append(vr[itheta, indx_pns])
        PNS_dictionary['local']['mflux'].append(mflux[itheta, indx_pns])
        PNS_dictionary['local']['m'].append(np.nancumsum(dmass,
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
            PNS_dictionary['local']['m'].append(np.nancumsum(dmass,
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
                                                    vth[PNS_mask] ** 2))
    PNS_dictionary['global']['ep'].append(np.nansum(0.5 * dmass[PNS_mask] *
                                                    vph[PNS_mask] ** 2))
    PNS_dictionary['global']['jx'].append(jx)
    PNS_dictionary['global']['jy'].append(jy)
    PNS_dictionary['global']['jz'].append(jz)
    PNS_dictionary['global']['jtot'].append(jtot)
    if simulation.dim == 2:
        itheta = np.arange(vth.shape[0])
        rperp2 = (simulation.cell.radius(simulation.ghost)[None, :] *
                          np.sin(simulation.cell.theta(simulation.ghost))[:, None]) ** 2
        omg = simulation.omega(file_name)
        PNS_dictionary['local']['vt'].append(vth[itheta, indx_pns])
        PNS_dictionary['local']['vp'].append(vph[itheta, indx_pns])
        PNS_dictionary['local']['omg'].append(omg[itheta, indx_pns])
        PNS_dictionary['global']['I'].append(np.nansum((dmass * 
                                                        rperp2)[PNS_mask]))
    elif simulation.dim == 3:
        iphi = np.arange(vth.shape[0])[:, None]
        itheta = np.arange(vth.shape[1])[None, :]
        PNS_dictionary['local']['vt'].append(vth[iphi, itheta, indx_pns])
        PNS_dictionary['local']['vp'].append(vph[iphi, itheta, indx_pns])
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        btot = np.sqrt(br ** 2 + bt ** 2 + bp ** 2)
        p = simulation.gas_pressure(file_name)
        if simulation.dim == 2:
            PNS_dictionary['local']['br'].append(br[itheta, indx_pns])
            PNS_dictionary['local']['bt'].append(bt[itheta, indx_pns])
            PNS_dictionary['local']['bp'].append(bp[itheta, indx_pns])
            PNS_dictionary['local']['beta'].append(p[itheta, indx_pns] / 
                                                   (0.5 * btot[itheta, indx_pns] ** 2 /
                                                    c.mu0))
        elif simulation.dim == 3:
            PNS_dictionary['local']['br'].append(br[iphi, itheta, indx_pns])
            PNS_dictionary['local']['bt'].append(bt[iphi, itheta, indx_pns])
            PNS_dictionary['local']['bp'].append(bp[iphi, itheta, indx_pns])
            PNS_dictionary['local']['beta'].append(p[iphi, itheta, indx_pns] / 
                                                  (0.5 * btot[iphi, itheta, indx_pns] ** 2 /
                                                   c.mu0))
        PNS_dictionary['global']['ebr'].append(np.nansum(0.5 * br[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        PNS_dictionary['global']['ebt'].append(np.nansum(0.5 * bt[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        PNS_dictionary['global']['ebp'].append(np.nansum(0.5 * bp[PNS_mask] ** 2 
                                                         * dV[PNS_mask] / c.mu0))
        PNS_dictionary['global']['bmax'].append(np.nanmax(btot[PNS_mask]))
        PNS_dictionary['global']['bavg'].append(np.nansum(btot[PNS_mask] * dmass[PNS_mask]) /
                                                PNS_dictionary['global']['mass'][-1])
        
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
    n = np.array([nx.value, ny.value, nz.value])
    if np.isnan(n).any():
        n = np.array([0, 0., 1.])
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    iphi = np.arange(vr.shape[0])[:, None]
    itheta = np.arange(vr.shape[1])[None, :]
    X = simulation.grid.cartesian_grid()
    x, y, z = X[0] * u.cm, X[1] * u.cm, X[2] * u.cm
    PNS_dictionary['global']['cmass_x'].append(np.nansum(x[PNS_mask] * dmass[PNS_mask]) / 
                                               PNS_dictionary['global']['mass'][-1])
    PNS_dictionary['global']['cmass_y'].append(np.nansum(y[PNS_mask] * dmass[PNS_mask]) / 
                                               PNS_dictionary['global']['mass'][-1])
    PNS_dictionary['global']['cmass_z'].append(np.nansum(z[PNS_mask] * dmass[PNS_mask]) / 
                                               PNS_dictionary['global']['mass'][-1])
    if n[2] >= 0.99:
        vrot = vp
        vturb = vt ** 2 + vr ** 2
        omg = simulation.omega(file_name)
        mask = PNS_mask
        R = (simulation.cell.radius(simulation.ghost)[None, None, :] * 
                 np.sin(simulation.cell.theta(simulation.ghost))[None, :, None])
    else:
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
    PNS_dictionary['global']['I'].append(np.nansum((dmass * 
                                                    R ** 2)[PNS_mask]))
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