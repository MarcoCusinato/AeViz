from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import PNS_core_radius
from AeViz.utils.math_utils import function_average_radii
from AeViz.units.constants import constants as c

def declare_PNScore_dictionary(dim: int,
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
    PNScore_dictionary = {
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
            'rho': [],
            'ye': [],
            't': [],
            'p': [],
            'vr': [],
            'mflux': [],
            'm': []
        }
    }
    if dim > 1:
        PNScore_dictionary['global']['et'] = []
        PNScore_dictionary['global']['ep'] = []
        PNScore_dictionary['global']['jx'] = []
        PNScore_dictionary['global']['jy'] = []
        PNScore_dictionary['global']['jz'] = []
        PNScore_dictionary['global']['jtot'] = []
        PNScore_dictionary['local']['vt'] = []
        PNScore_dictionary['local']['vp'] = []
        PNScore_dictionary['local']['omg'] = []
        if magdim > 0:
            PNScore_dictionary['global']['ebr'] = []
            PNScore_dictionary['global']['ebt'] = []
            PNScore_dictionary['global']['ebp'] = []
            PNScore_dictionary['global']['bmax'] = []
            PNScore_dictionary['global']['bavg'] = []
            PNScore_dictionary['local']['br'] = []
            PNScore_dictionary['local']['bt'] = []
            PNScore_dictionary['local']['bp'] = []
            PNScore_dictionary['local']['beta'] = []
    if dim > 2:
        PNScore_dictionary['global']['erot'] = []
        PNScore_dictionary['global']['eres'] = []
        PNScore_dictionary['global']['cmass_x'] = []
        PNScore_dictionary['global']['cmass_y'] = []
        PNScore_dictionary['global']['cmass_z'] = []
        PNScore_dictionary['local']['vrot']  = []
        if magdim > 0:
            PNScore_dictionary['global']['ebpol'] = []
            PNScore_dictionary['global']['ebtor'] = []
            PNScore_dictionary['local']['btor'] = []
            PNScore_dictionary['local']['bpol'] = []
            
    return PNScore_dictionary    

def PNScore_quantity_metadata() -> dict:
    """
    Metadata for PNS post-processing quantities.
    """

    meta = {
        "global": {
            "mass":   dict(name="core_mass",
                           label=r"$M_{\rm core}$",
                           limits=[1.0, 2.5],
                           cmap="viridis",
                           log=False),
            "rmax":   dict(name="core_radius_max",
                           label=r"$R_{\rm core, max}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "rmin":   dict(name="core_radius_min",
                           label=r"$R_{\rm core, min}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "ravg":   dict(name="core_radius_avg",
                           label=r"$\langle R_{\rm core}\rangle$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "mflux":  dict(name="core_mass_flux",
                           label=r"$\dot{M}_{\rm core}$",
                           limits=[1e-4, 5],
                           cmap="cividis",
                           log=True),
            "er":     dict(name="radial_kinetic_energy",
                           label=r"$E_{\rm core,kin,r}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "et":     dict(name="theta_kinetic_energy",
                           label=r"$E_{\rm core,kin,\theta}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "ep":     dict(name="phi_kinetic_energy",
                           label=r"$E_{\rm core,kin,\phi}$",
                           limits=[1e46, 1e52],
                           cmap="inferno",
                           log=True),
            "eg":     dict(name="gravitational_energy",
                           label=r"$E_{\rm core,grav}$",
                           limits=[-5e53, -1e51],
                           cmap="magma",
                           log=False),
            "ei":     dict(name="internal_energy",
                           label=r"$E_{\rm core,int}$",
                           limits=[1e50, 5e53],
                           cmap="hot",
                           log=True),
            "I":      dict(name="moment_of_inertia",
                           label=r"$I_{\rm core}$",
                           limits=[1e43, 5e45],
                           cmap="Greens",
                           log=True),
            "yeavg":  dict(name="Ye_average",
                           label=r"$\langle Y_e\rangle_{\rm core}$",
                           limits=[0.0, 0.5],
                           cmap="RdYlBu",
                           log=False),
            "tavg":   dict(name="temperature_average",
                           label=r"$\langle T\rangle_{\rm core}$",
                           limits=[1, 60],
                           cmap="afmhot",
                           log=False),
            "savg":   dict(name="entropy_average",
                           label=r"$\langle s\rangle_{\rm core}$",
                           limits=[0, 30],
                           cmap="Spectral",
                           log=False),
            "jx":   dict(name="Jx",
                         label=r"$J_{x,\rm core}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jy":   dict(name="Jy",
                         label=r"$J_{y,\rm core}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jz":   dict(name="Jz",
                         label=r"$J_{z,\rm core}$",
                         limits=[-5e49, 5e49],
                         cmap="coolwarm",
                         log=False),
            "jtot": dict(name="Jtot",
                         label=r"$|\mathbf{J}|_{\rm core}$",
                         limits=[1e46, 1e50],
                         cmap="viridis",
                         log=True),
            "ebr": dict(name="magnetic_energy_radial",
                        label=r"$E_{\rm mag,r,core}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "ebt": dict(name="magnetic_energy_theta",
                        label=r"$E_{\rm mag,\theta,core}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "ebp": dict(name="magnetic_energy_phi",
                        label=r"$E_{\rm mag,\phi,core}$",
                        limits=[1e42, 1e51],
                        cmap="PuBuGn",
                        log=True),
            "bmax": dict(name="maximum_magnetic_field",
                         label=r"$B_{\rm max,core}$",
                         limits=[1e10, 1e17],
                         cmap="PuRd",
                         log=True),
            "bavg": dict(name="average_magnetic_field",
                         label=r"$\langle B\rangle_{\rm core}$",
                         limits=[1e9, 1e16],
                         cmap="PuRd",
                         log=True),
            "erot": dict(name="rotational_energy",
                         label=r"$E_{\rm rot,core}$",
                         limits=[1e46, 1e53],
                         cmap="Oranges",
                         log=True),
            "eres": dict(name="turbulent_energy",
                         label=r"$E_{\rm turb,core}$",
                         limits=[1e46, 1e53],
                         cmap="YlOrBr",
                         log=True),
            "cmass_x": dict(name="center_of_mass_x",
                            label=r"$x_{\rm COM,core}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "cmass_y": dict(name="center_of_mass_y",
                            label=r"$y_{\rm COM,core}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "cmass_z": dict(name="center_of_mass_z",
                            label=r"$z_{\rm COM,core}$",
                            limits=[-20, 20],
                            cmap="coolwarm",
                            log=False),
            "ebpol": dict(name="poloidal_magnetic_energy",
                            label=r"$E_{\rm mag,pol,core}$",
                            limits=[1e42, 1e51],
                            cmap="PuBuGn",
                            log=True),
            "ebtor": dict(name="toroidal_magnetic_energy",
                            label=r"$E_{\rm mag,tor,core}$",
                            limits=[1e42, 1e51],
                            cmap="PuBuGn",
                            log=True),
        },

        "local": {
            "r":      dict(name="core_radius_surface",
                           label=r"$R_{\rm core}$",
                           limits=[5, 60],
                           cmap="plasma",
                           log=False),
            "rho":      dict(name="surface_density",
                           label=r"$\rho_{\rm core}$",
                           limits=[0, 40],
                           cmap="Spectral",
                           log=False),
            "ye":     dict(name="surface_Ye",
                           label=r"$Y_{e,\rm core}$",
                           limits=[0, 0.6],
                           cmap="RdYlBu",
                           log=False),
            "t":      dict(name="surface_temperature",
                           label=r"$T_{\rm core}$",
                           limits=[1, 60],
                           cmap="afmhot",
                           log=False),
            "p":      dict(name="surface_pressure",
                           label=r"$P_{\rm gas,core}$",
                           limits=[1e27, 1e35],
                           cmap="magma",
                           log=True),
            "vr":     dict(name="surface_radial_velocity",
                           label=r"$v_{r,\rm core}$",
                           limits=[-0.5, 0.5],
                           cmap="coolwarm",
                           log=False),
            "mflux":  dict(name="surface_mass_flux",
                           label=r"$\dot{M}_{\rm core}$",
                           limits=[1e-4, 5],
                           cmap="cividis",
                           log=True),
            "m":      dict(name="enclosed_mass",
                           label=r"$M(<R_{\rm core})$",
                           limits=[1.0, 2.5],
                           cmap="viridis",
                           log=False),
            "vt":  dict(name="surface_theta_velocity",
                        label=r"$v_{\theta,\rm core}$",
                        limits=[-0.5, 0.5],
                        cmap="coolwarm",
                        log=False),
            "vp":  dict(name="surface_phi_velocity",
                        label=r"$v_{\phi,\rm core}$",
                        limits=[-0.5, 0.5],
                        cmap="coolwarm",
                        log=False),
            "omg": dict(name="surface_omega",
                        label=r"$\Omega_{\rm core}$",
                        limits=[1, 1e4],
                        cmap="twilight",
                        log=True),
            "br": dict(name="surface_Br",
                       label=r"$B_{r,\rm core}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "bt": dict(name="surface_Btheta",
                       label=r"$B_{\theta,\rm core}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "bp": dict(name="surface_Bphi",
                       label=r"$B_{\phi,\rm core}$",
                       limits=[-1e16, 1e16],
                       cmap="RdBu_r",
                       log=False),
            "beta": dict(name="surface_plasma_beta",
                         label=r"$\beta_{\rm core}$",
                         limits=[1e-2, 1e4],
                         cmap="viridis",
                         log=True),
            "vrot": dict(
                        name="surface_rotational_velocity",
                        label=r"$v_{\rm rot,core}$",
                        limits=[0, 0.5],
                        cmap="twilight_shifted",
                        log=False,
                        ),
            "bpol": dict(name="surface_Bpol",
                         label=r"$B_{\rm pol,core}$",
                         limits=[0, 1e16],
                         cmap="PuRd",
                         log=True),
            
            "btor": dict(name="surface_Btor",
                         label=r"$B_{\rm tor,core}$",
                         limits=[0, 1e16],
                         cmap="PuRd",
                         log=True),   
        }
    }

    return meta

def compute_PNScore_postprocessing(PNScore_dictionary: dict,
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
    PNScore_dictionary : dict
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
    PNScr = PNS_core_radius(simulation, file_name)
    r = simulation.cell.radius(simulation.ghost)
    if simulation.dim == 1:
        PNScmask = (r <= PNScr)
        idx_pns = np.argmax(r >= PNScr)
    else:
        while r.ndim <= PNScr.ndim:
            r = r[None, :]
        PNScmask = (r <= PNScr)
        idx_cpns = np.argmax(r >= PNScr[..., None], axis=-1)
    dmass = (simulation.rho(file_name) * dV).to(u.M_sun)
    __compute_1D_global_local_quantities(PNScore_dictionary,
                                         simulation,
                                         file_name,
                                         dV,
                                         dOmega,
                                         PNScmask,
                                         idx_pns,
                                         PNScr,
                                         dmass)
    if simulation.dim > 1:
        __compute_2D_global_local_quantities(PNScore_dictionary,
                                            simulation,
                                            file_name,
                                            PNScmask,
                                            idx_pns,
                                            dV,
                                            Lx,
                                            Ly,
                                            Lz,
                                            dmass)
    if simulation.dim == 3:
        __compute_3D_global_local_quantities(PNScore_dictionary,
                                            simulation,
                                            file_name,
                                            PNScmask,
                                            idx_pns,
                                            dmass,
                                            dV)

def __compute_1D_global_local_quantities(PNScore_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       dV: aerray,
                                       dOmega: aerray,
                                       PNScore_mask: np.array[bool],
                                       indx_pns: np.array[int],
                                       PNScore_radius: aerray,
                                       dmass:aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNScore_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    dV : aerray
        _description_
    dOmega : aerray
        _description_
    PNScore_mask : np.array[bool]
        _description_
    indx_pns : np.array[int]
        _description_
    PNScore_radius : aerray
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
    rho = simulation.rho(file_name)

    mass = np.nansum(dmass[PNScore_mask])
    ## add the global variables
    PNScore_dictionary['global']['mass'].append(mass)
    PNScore_dictionary['global']['rmax'].append(np.nanmax(PNScore_radius))
    PNScore_dictionary['global']['rmin'].append(np.nanmin(PNScore_radius))
    PNScore_dictionary['global']['ravg'].append(function_average_radii(PNScore_radius,
                                                                   simulation.dim,
                                                                   dOmega))
    PNScore_dictionary['global']['er'].append(np.nansum(0.5 * dmass[PNScore_mask] * 
                                                    vr[PNScore_mask] ** 2))
    PNScore_dictionary['global']['eg'].append(np.nansum(grav_en[PNScore_mask]))
    PNScore_dictionary['global']['ei'].append(np.nansum(int_en[PNScore_mask]))
    PNScore_dictionary['global']['yeavg'].append(np.nansum(ye[PNScore_mask] * dmass[PNScore_mask]) /
                                             mass)
    PNScore_dictionary['global']['savg'].append(np.nansum(s[PNScore_mask] * dmass[PNScore_mask]) /
                                                 mass)
    PNScore_dictionary['global']['tavg'].append(np.nansum(t[PNScore_mask] * dmass[PNScore_mask]) /
                                                     mass)

    
    
    ## add variables at the surface
    PNScore_dictionary['local']['r'].append(PNScore_radius)
    if simulation.dim == 1:
        PNScore_dictionary['local']['rho'].append(rho[indx_pns])
        PNScore_dictionary['local']['ye'].append(ye[indx_pns])
        PNScore_dictionary['local']['t'].append(t[indx_pns])
        PNScore_dictionary['local']['p'].append(p[indx_pns])
        PNScore_dictionary['local']['vr'].append(vr[indx_pns])
        PNScore_dictionary['local']['mflux'].append(mflux[indx_pns])
        PNScore_dictionary['local']['m'].append(mass)
        PNScore_dictionary['global']['mflux'].append(np.nansum(mflux[indx_pns] *
                                                           dOmega))
        rperp2 = simulation.cell.radius(simulation.ghost) ** 2
        PNScore_dictionary['global']['I'].append(np.nansum((dmass * 
                                                        rperp2)[PNScore_mask]))
    elif simulation.dim == 2:
        itheta = np.arange(p.shape[0])
        PNScore_dictionary['local']['s'].append(rho[itheta, indx_pns])
        PNScore_dictionary['local']['ye'].append(ye[itheta, indx_pns])
        PNScore_dictionary['local']['t'].append(t[itheta, indx_pns])
        PNScore_dictionary['local']['p'].append(p[itheta, indx_pns])
        PNScore_dictionary['local']['vr'].append(vr[itheta, indx_pns])
        PNScore_dictionary['local']['mflux'].append(mflux[itheta, indx_pns])
        PNScore_dictionary['local']['m'].append(np.nancumsum(dmass,
                                                         axis=-1)[itheta,
                                                            indx_pns])
        PNScore_dictionary['global']['mflux'].append(np.nansum(mflux[itheta,
                                                                 indx_pns] *
                                                                   dOmega))
    elif simulation.dim == 3:
            itheta = np.arange(p.shape[1])
            iphi = np.arange(p.shape[0])
            PNScore_dictionary['local']['s'].append(rho[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['ye'].append(ye[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['t'].append(t[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['p'].append(p[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['vr'].append(vr[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['mflux'].append(mflux[iphi, itheta, 
                                                          indx_pns])
            PNScore_dictionary['local']['m'].append(np.nancumsum(dmass,
                                                             axis=-1)[iphi,
                                                                      itheta,
                                                                      indx_pns])
            PNScore_dictionary['global']['mflux'].append(np.nansum(mflux[iphi,
                                                                     itheta,
                                                                     indx_pns] *
                                                                       dOmega))
            
def __compute_2D_global_local_quantities(PNScore_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       PNScore_mask: np.array[bool],
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
    PNScore_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    PNScore_mask : np.array[bool]
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
    jx = np.nansum(Lx[PNScore_mask])
    jy = np.nansum(Ly[PNScore_mask])
    jz = np.nansum(Lz[PNScore_mask])
    jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
    PNScore_dictionary['global']['et'].append(np.nansum(0.5 * dmass[PNScore_mask] * 
                                                    vth[PNScore_mask] ** 2))
    PNScore_dictionary['global']['ep'].append(np.nansum(0.5 * dmass[PNScore_mask] *
                                                    vph[PNScore_mask] ** 2))
    PNScore_dictionary['global']['jx'].append(jx)
    PNScore_dictionary['global']['jy'].append(jy)
    PNScore_dictionary['global']['jz'].append(jz)
    PNScore_dictionary['global']['jtot'].append(jtot)
    if simulation.dim == 2:
        itheta = np.arange(vth.shape[0])
        rperp2 = (simulation.cell.radius(simulation.ghost)[None, :] *
                          np.sin(simulation.cell.theta(simulation.ghost))[:, None]) ** 2
        omg = simulation.omega(file_name)
        PNScore_dictionary['local']['vt'].append(vth[itheta, indx_pns])
        PNScore_dictionary['local']['vp'].append(vph[itheta, indx_pns])
        PNScore_dictionary['local']['omg'].append(omg[itheta, indx_pns])
        PNScore_dictionary['global']['I'].append(np.nansum((dmass * 
                                                        rperp2)[PNScore_mask]))
    elif simulation.dim == 3:
        iphi = np.arange(vth.shape[0])[:, None]
        itheta = np.arange(vth.shape[1])[None, :]
        PNScore_dictionary['local']['vt'].append(vth[iphi, itheta, indx_pns])
        PNScore_dictionary['local']['vp'].append(vph[iphi, itheta, indx_pns])
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        btot = np.sqrt(br ** 2 + bt ** 2 + bp ** 2)
        p = simulation.gas_pressure(file_name)
        if simulation.dim == 2:
            PNScore_dictionary['local']['br'].append(br[itheta, indx_pns])
            PNScore_dictionary['local']['bt'].append(bt[itheta, indx_pns])
            PNScore_dictionary['local']['bp'].append(bp[itheta, indx_pns])
            PNScore_dictionary['local']['beta'].append(p[itheta, indx_pns] / 
                                                   (0.5 * btot[itheta, indx_pns] ** 2 /
                                                    c.mu0))
        elif simulation.dim == 3:
            PNScore_dictionary['local']['br'].append(br[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['bt'].append(bt[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['bp'].append(bp[iphi, itheta, indx_pns])
            PNScore_dictionary['local']['beta'].append(p[iphi, itheta, indx_pns] / 
                                                  (0.5 * btot[iphi, itheta, indx_pns] ** 2 /
                                                   c.mu0))
        PNScore_dictionary['global']['ebr'].append(np.nansum(0.5 * br[PNScore_mask] ** 2 
                                                         * dV[PNScore_mask] / c.mu0))
        PNScore_dictionary['global']['ebt'].append(np.nansum(0.5 * bt[PNScore_mask] ** 2 
                                                         * dV[PNScore_mask] / c.mu0))
        PNScore_dictionary['global']['ebp'].append(np.nansum(0.5 * bp[PNScore_mask] ** 2 
                                                         * dV[PNScore_mask] / c.mu0))
        PNScore_dictionary['global']['bmax'].append(np.nanmax(btot[PNScore_mask]))
        PNScore_dictionary['global']['bavg'].append(np.nansum(btot[PNScore_mask] * dmass[PNScore_mask]) /
                                                PNScore_dictionary['global']['mass'][-1])
        
def __compute_3D_global_local_quantities(PNScore_dictionary: dict,
                                         simulation: Simulation,
                                         file_name: str,
                                         PNScore_mask: np.array[bool],
                                         indx_pns: np.array[int],
                                         dmass: aerray,
                                         dV: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    PNScore_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    PNScore_mask : np.array[bool]
        _description_
    indx_pns : np.array[int]
        _description_
    dmass : aerray
        _description_
    dV : aerray
        _description_
    """
    jtot = PNScore_dictionary['global']['jtot'][-1].value
    nx = PNScore_dictionary['global']['jx'][-1] / jtot
    ny = PNScore_dictionary['global']['jy'][-1] / jtot
    nz = PNScore_dictionary['global']['jz'][-1] / jtot
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
    PNScore_dictionary['global']['cmass_x'].append(np.nansum(x[PNScore_mask] * dmass[PNScore_mask]) / 
                                               PNScore_dictionary['global']['mass'][-1])
    PNScore_dictionary['global']['cmass_y'].append(np.nansum(y[PNScore_mask] * dmass[PNScore_mask]) / 
                                               PNScore_dictionary['global']['mass'][-1])
    PNScore_dictionary['global']['cmass_z'].append(np.nansum(z[PNScore_mask] * dmass[PNScore_mask]) / 
                                               PNScore_dictionary['global']['mass'][-1])
    if n[2] >= 0.99:
        vrot = vp
        vturb = vt ** 2 + vr ** 2
        omg = simulation.omega(file_name)
        mask = PNScore_mask
        R = (simulation.cell.radius(simulation.ghost)[None, None, :] * 
                 np.sin(simulation.cell.theta(simulation.ghost))[None, :, None])
    else:
        X_perp = X - np.tensordot(n, X, axes=(0, 0))[None, ...] * \
            n[:, None, None, None]
        R = np.linalg.norm(X_perp, axis=0) + 1e-12
        mask = ((R > simulation.grid.radius[0].value) & (PNScore_mask))
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
    PNScore_dictionary['global']['I'].append(np.nansum((dmass * 
                                                    R ** 2)[PNScore_mask]))
    PNScore_dictionary['global']['erot'].append(np.nansum(0.5 * (vrot ** 2 * 
                                                             dmass)[mask]))
    PNScore_dictionary['global']['eres'].append(np.nansum(0.5 * (vturb *
                                                             dmass)[mask]))
    PNScore_dictionary['local']['vrot'].append(vrot[iphi, itheta, indx_pns])
    PNScore_dictionary['local']['omg'].append(omg[iphi, itheta, indx_pns])
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
        PNScore_dictionary['global']['ebpol'].append(np.nansum(0.5 * (bpol ** 2 * 
                                                                     dV)[mask]) 
                                                 / c.mu0)
        PNScore_dictionary['global']['ebtor'].append(np.nansum(0.5 * (btor ** 2 * 
                                                                    dV)[mask]) 
                                                         / c.mu0)
        PNScore_dictionary['local']['bpol'].append(bpol[iphi, itheta, indx_pns])
        PNScore_dictionary['local']['btor'].append(btor[iphi, itheta, indx_pns])