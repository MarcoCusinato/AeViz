from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import gain_radius
from AeViz.units.constants import constants as c

def declare_gain_dictionary(dim: int,
                            magdim: int) -> dict:
    """
    Declare the dictionary used to store gain-region quantities.

    Parameters
    ----------
    dim : int
        Number of spatial dimensions.
    magdim : int
        Number of magnetic-field dimensions.

    Returns
    -------
    dict
        Nested dictionary containing global and local gain-region
        quantities.
    """
    gain_dictionary = {
        'global': {
            'mass': [],
            'mass_unb': [],

            'rmax': [],
            'rmin': [],
            'ravg': [],

            'tadv': [],
            'tadv_unb': [],

            'er': [],
            'er_unb': [],
            'eg': [],
            'eg_unb': [],
            'ei': [],
            'ei_unb': [],
            'nu_heat': [],
            'nu_heat_unb': [],

            'mdot_net': [],
            'mdot_in': [],
            'mdot_out': [],
        },
        'local': {
            'r': [],
            'mass': [],
            'mass_unb': [],

            'tadv': [],
            'tadv_unb': [],

            'er': [],
            'er_unb': [],
            'eg': [],
            'eg_unb': [],
            'ei': [],
            'ei_unb': [],
            'nu_heat': [],
            'nu_heat_unb': [],
        }
    }

    if dim > 1:
        gain_dictionary['global']['et'] = []
        gain_dictionary['global']['ep'] = []
        gain_dictionary['global']['et_unb'] = []
        gain_dictionary['global']['ep_unb'] = []

        gain_dictionary['local']['et'] = []
        gain_dictionary['local']['ep'] = []
        gain_dictionary['local']['et_unb'] = []
        gain_dictionary['local']['ep_unb'] = []

        gain_dictionary['global']['jx'] = []
        gain_dictionary['global']['jy'] = []
        gain_dictionary['global']['jz'] = []
        gain_dictionary['global']['jtot'] = []

        gain_dictionary['global']['jx_unb'] = []
        gain_dictionary['global']['jy_unb'] = []
        gain_dictionary['global']['jz_unb'] = []
        gain_dictionary['global']['jtot_unb'] = []

        if magdim > 0:
            for key in ('ebr', 'ebt', 'ebp',
                        'ebr_unb', 'ebt_unb', 'ebp_unb'):
                gain_dictionary['global'][key] = []
                gain_dictionary['local'][key] = []

    return gain_dictionary

def gain_quantity_metadata() -> dict:
    """
    Metadata for gain-region post-processing quantities.
    """

    meta = {
        "global": {
            "mass": dict(
                name="gain_mass",
                label=r"$M_{\rm gain}$",
                limits=[1e-5, 1.0],
                cmap="viridis",
                log=True
            ),
            "mass_unb": dict(
                name="gain_unbound_mass",
                label=r"$M_{\rm gain,unb}$",
                limits=[1e-6, 1.0],
                cmap="viridis",
                log=True
            ),

            "rmax": dict(
                name="gain_radius_max",
                label=r"$R_{\rm gain,max}$",
                limits=[1e4, 2e7],
                cmap="plasma",
                log=False
            ),
            "rmin": dict(
                name="gain_radius_min",
                label=r"$R_{\rm gain,min}$",
                limits=[1e4, 2e7],
                cmap="plasma",
                log=False
            ),
            "ravg": dict(
                name="gain_radius_avg",
                label=r"$\langle R_{\rm gain}\rangle$",
                limits=[1e4, 2e7],
                cmap="plasma",
                log=False
            ),

            "tadv": dict(
                name="gain_advection_time",
                label=r"$\tau_{\rm adv,gain}$",
                limits=[1e-4, 1e2],
                cmap="magma",
                log=True
            ),
            "tadv_unb": dict(
                name="gain_unbound_advection_time",
                label=r"$\tau_{\rm adv,gain,unb}$",
                limits=[1e-4, 1e2],
                cmap="magma",
                log=True
            ),

            "er": dict(
                name="gain_radial_kinetic_energy",
                label=r"$E_{\rm kin,r,gain}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),
            "er_unb": dict(
                name="gain_unbound_radial_kinetic_energy",
                label=r"$E_{\rm kin,r,gain,unb}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),

            "eg": dict(
                name="gain_gravitational_energy",
                label=r"$E_{\rm grav,gain}$",
                limits=[-1e52, -1e44],
                cmap="magma",
                log=False
            ),
            "eg_unb": dict(
                name="gain_unbound_gravitational_energy",
                label=r"$E_{\rm grav,gain,unb}$",
                limits=[-1e52, -1e44],
                cmap="magma",
                log=False
            ),

            "ei": dict(
                name="gain_internal_energy",
                label=r"$E_{\rm int,gain}$",
                limits=[1e44, 1e52],
                cmap="hot",
                log=True
            ),
            "ei_unb": dict(
                name="gain_unbound_internal_energy",
                label=r"$E_{\rm int,gain,unb}$",
                limits=[1e44, 1e52],
                cmap="hot",
                log=True
            ),

            "nu_heat": dict(
                name="gain_neutrino_heating_rate",
                label=r"$\dot Q_{\nu,\rm gain}$",
                limits=[-1e52, 1e52],
                cmap="coolwarm",
                log=False
            ),
            "nu_heat_unb": dict(
                name="gain_unbound_neutrino_heating_rate",
                label=r"$\dot Q_{\nu,\rm gain,unb}$",
                limits=[-1e52, 1e52],
                cmap="coolwarm",
                log=False
            ),

            "mdot_net": dict(
                name="gain_mass_flux_net",
                label=r"$\dot M_{\rm net}$",
                limits=[-5.0, 5.0],
                cmap="coolwarm",
                log=False
            ),
            "mdot_in": dict(
                name="gain_mass_flux_inflow",
                label=r"$\dot M_{\rm in}$",
                limits=[1e-5, 5.0],
                cmap="viridis",
                log=True
            ),
            "mdot_out": dict(
                name="gain_mass_flux_outflow",
                label=r"$\dot M_{\rm out}$",
                limits=[1e-5, 5.0],
                cmap="plasma",
                log=True
            ),

            "et": dict(
                name="gain_theta_kinetic_energy",
                label=r"$E_{\rm kin,\theta,gain}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),
            "ep": dict(
                name="gain_phi_kinetic_energy",
                label=r"$E_{\rm kin,\phi,gain}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),
            "et_unb": dict(
                name="gain_unbound_theta_kinetic_energy",
                label=r"$E_{\rm kin,\theta,gain,unb}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),
            "ep_unb": dict(
                name="gain_unbound_phi_kinetic_energy",
                label=r"$E_{\rm kin,\phi,gain,unb}$",
                limits=[1e44, 1e52],
                cmap="inferno",
                log=True
            ),

            "jx": dict(
                name="gain_jx",
                label=r"$J_{x,\rm gain}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jy": dict(
                name="gain_jy",
                label=r"$J_{y,\rm gain}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jz": dict(
                name="gain_jz",
                label=r"$J_{z,\rm gain}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jtot": dict(
                name="gain_jtot",
                label=r"$|\mathbf{J}_{\rm gain}|$",
                limits=[1e46, 1e50],
                cmap="viridis",
                log=True
            ),
            "jx_unb": dict(
                name="gain_jx_unb",
                label=r"$J_{x,\rm gain,unb}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jy_unb": dict(
                name="gain_jy_unb",
                label=r"$J_{y,\rm gain,unb}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jz_unb": dict(
                name="gain_jz_unb",
                label=r"$J_{z,\rm gain,unb}$",
                limits=[-5e49, 5e49],
                cmap="coolwarm",
                log=False
            ),
            "jtot_unb": dict(
                name="gain_jtot_unb",
                label=r"$|\mathbf{J}_{\rm gain,unb}|$",
                limits=[1e46, 1e50],
                cmap="viridis",
                log=True
            ),

            "ebr": dict(
                name="gain_magnetic_energy_radial",
                label=r"$E_{\rm mag,r,gain}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebt": dict(
                name="gain_magnetic_energy_theta",
                label=r"$E_{\rm mag,\theta,gain}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebp": dict(
                name="gain_magnetic_energy_phi",
                label=r"$E_{\rm mag,\phi,gain}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebr_unb": dict(
                name="gain_unbound_magnetic_energy_radial",
                label=r"$E_{\rm mag,r,gain,unb}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebt_unb": dict(
                name="gain_unbound_magnetic_energy_theta",
                label=r"$E_{\rm mag,\theta,gain,unb}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebp_unb": dict(
                name="gain_unbound_magnetic_energy_phi",
                label=r"$E_{\rm mag,\phi,gain,unb}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),

            "ebpol": dict(
                name="gain_poloidal_magnetic_energy",
                label=r"$E_{\rm mag,pol,gain}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebtor": dict(
                name="gain_toroidal_magnetic_energy",
                label=r"$E_{\rm mag,tor,gain}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebpol_unb": dict(
                name="gain_unbound_poloidal_magnetic_energy",
                label=r"$E_{\rm mag,pol,gain,unb}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
            "ebtor_unb": dict(
                name="gain_unbound_toroidal_magnetic_energy",
                label=r"$E_{\rm mag,tor,gain,unb}$",
                limits=[1e40, 1e50],
                cmap="PuBuGn",
                log=True
            ),
        },

        "local": {
            "r": dict(
                name="gain_radius",
                label=r"$r_{\rm gain}$",
                limits=[1e4, 2e7],
                cmap="plasma",
                log=False
            ),
            "mass": dict(
                name="gain_mass_profile",
                label=r"$M_{\rm gain}(\theta,\phi)$",
                limits=[1e-6, 1.0],
                cmap="viridis",
                log=True
            ),
            "mass_unb": dict(
                name="gain_unbound_mass_profile",
                label=r"$M_{\rm gain,unb}(\theta,\phi)$",
                limits=[1e-6, 1.0],
                cmap="viridis",
                log=True
            ),

            "tadv": dict(
                name="gain_advection_time_profile",
                label=r"$\tau_{\rm adv,gain}(\theta,\phi)$",
                limits=[1e-4, 1e2],
                cmap="magma",
                log=True
            ),
            "tadv_unb": dict(
                name="gain_unbound_advection_time_profile",
                label=r"$\tau_{\rm adv,gain,unb}(\theta,\phi)$",
                limits=[1e-4, 1e2],
                cmap="magma",
                log=True
            ),

            "er": dict(
                name="gain_radial_kinetic_energy_profile",
                label=r"$E_{\rm kin,r,gain}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),
            "er_unb": dict(
                name="gain_unbound_radial_kinetic_energy_profile",
                label=r"$E_{\rm kin,r,gain,unb}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),
            "eg": dict(
                name="gain_gravitational_energy_profile",
                label=r"$E_{\rm grav,gain}(\theta,\phi)$",
                limits=[-1e50, -1e42],
                cmap="magma",
                log=False
            ),
            "eg_unb": dict(
                name="gain_unbound_gravitational_energy_profile",
                label=r"$E_{\rm grav,gain,unb}(\theta,\phi)$",
                limits=[-1e50, -1e42],
                cmap="magma",
                log=False
            ),
            "ei": dict(
                name="gain_internal_energy_profile",
                label=r"$E_{\rm int,gain}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="hot",
                log=True
            ),
            "ei_unb": dict(
                name="gain_unbound_internal_energy_profile",
                label=r"$E_{\rm int,gain,unb}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="hot",
                log=True
            ),
            "nu_heat": dict(
                name="gain_neutrino_heating_profile",
                label=r"$\dot Q_{\nu,\rm gain}(\theta,\phi)$",
                limits=[-1e50, 1e50],
                cmap="coolwarm",
                log=False
            ),
            "nu_heat_unb": dict(
                name="gain_unbound_neutrino_heating_profile",
                label=r"$\dot Q_{\nu,\rm gain,unb}(\theta,\phi)$",
                limits=[-1e50, 1e50],
                cmap="coolwarm",
                log=False
            ),

            "et": dict(
                name="gain_theta_kinetic_energy_profile",
                label=r"$E_{\rm kin,\theta,gain}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),
            "ep": dict(
                name="gain_phi_kinetic_energy_profile",
                label=r"$E_{\rm kin,\phi,gain}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),
            "et_unb": dict(
                name="gain_unbound_theta_kinetic_energy_profile",
                label=r"$E_{\rm kin,\theta,gain,unb}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),
            "ep_unb": dict(
                name="gain_unbound_phi_kinetic_energy_profile",
                label=r"$E_{\rm kin,\phi,gain,unb}(\theta,\phi)$",
                limits=[1e42, 1e50],
                cmap="inferno",
                log=True
            ),

            "ebr": dict(
                name="gain_magnetic_energy_radial_profile",
                label=r"$E_{\rm mag,r,gain}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
            "ebt": dict(
                name="gain_magnetic_energy_theta_profile",
                label=r"$E_{\rm mag,\theta,gain}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
            "ebp": dict(
                name="gain_magnetic_energy_phi_profile",
                label=r"$E_{\rm mag,\phi,gain}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
            "ebr_unb": dict(
                name="gain_unbound_magnetic_energy_radial_profile",
                label=r"$E_{\rm mag,r,gain,unb}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
            "ebt_unb": dict(
                name="gain_unbound_magnetic_energy_theta_profile",
                label=r"$E_{\rm mag,\theta,gain,unb}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
            "ebp_unb": dict(
                name="gain_unbound_magnetic_energy_phi_profile",
                label=r"$E_{\rm mag,\phi,gain,unb}(\theta,\phi)$",
                limits=[1e38, 1e48],
                cmap="PuBuGn",
                log=True
            ),
        }
    }

    return meta

def compute_gain_postprocessing(gain_dictionary: dict,
                                simulation: Simulation,
                                file_name: str,
                                dV: aerray,
                                PNS_radius: aerray,
                                shock_radius: aerray,
                                Lx: aerray | None,
                                Ly: aerray | None,
                                Lz: aerray | None,) -> None:
    gainr = gain_radius(simulation,
                        file_name,
                        PNS_radius)
    radius = simulation.grid.radius
    if simulation.dim >1:
        while radius.ndim <= gainr.ndim:
            radius = radius[None, :]
    mhd_ene = simulation.MHD_energy(file_name) + \
            2 * simulation.gravitational_energy(file_name)
    gain_mask = (radius >= gainr) & (radius <= shock_radius)
    unbound_mask = gain_mask & (mhd_ene > 0) & \
        (simulation.radial_velocity(file_name) > 0)
    __compute_1D_global_local_quantities(gain_dictionary,
                                         simulation,
                                         file_name,
                                         dV,
                                         gain_radius,
                                         gain_mask,
                                         unbound_mask)
    if simulation.dim > 1:
        __compute_2D_global_local_quantities(gain_dictionary,
                                             simulation,
                                             file_name,
                                             gain_mask,
                                             unbound_mask,
                                             dV,
                                             Lx,
                                             Ly,
                                             Lz)

def __compute_1D_global_local_quantities(gain_dictionary: dict,
                                         simulation: Simulation,
                                         file_name: str,
                                         dV: aerray,
                                         gain_radius: aerray,
                                         gain_mask: np.array[bool],
                                         unbound_mask: np.array[bool]) -> None:
    rho = simulation.rho(file_name)
    vr = simulation.radial_velocity(file_name)
    grav_en = simulation.gravitational_energy(file_name) * dV
    int_en = simulation.internal_energy(file_name) * dV
    nu_heat = simulation.nu_heat(file_name) * dV
    dr = simulation.cell.dr(simulation.ghost)

    mass = rho * dV
    ekin_r = 0.5 * vr ** 2 * mass
    tadv_gain = np.nansum(dr / np.abs(vr) * gain_mask, axis=-1)
    tadv_unb = np.nansum(dr / np.abs(vr) * unbound_mask, axis=-1)
    tadv_gain = np.where(tadv_gain == 0, np.nan, tadv_gain)
    tadv_unb = np.where(tadv_unb == 0, np.nan, tadv_unb)
    ## add global values
    gain_dictionary['global']['mass'].append(np.nansum(mass[gain_mask]))
    gain_dictionary['global']['mass_unb'].append(np.nansum(mass[unbound_mask]))
    gain_dictionary['global']['rmax'].append(np.nanmax(gain_radius))
    gain_dictionary['global']['rmin'].append(np.nanmin(gain_radius))
    gain_dictionary['global']['ravg'].append(simulation.grid.radii_angular_average(gain_radius))
    gain_dictionary['global']['tadv'].append(simulation.grid.radii_angular_average(tadv_gain))
    gain_dictionary['global']['tadv_unb'].append(simulation.grid.radii_angular_average(tadv_unb))

    gain_dictionary['local']['r'].append(gain_radius)
    gain_dictionary['local']['er'].append(np.nansum(gain_mask * ekin_r * gain_mask, axis=-1))
    gain_dictionary['local']['er_unb'].append(np.nansum(gain_mask * ekin_r * unbound_mask, axis=-1))
    gain_dictionary['local']['eg'].append(np.nansum(gain_mask * grav_en * gain_mask, axis=-1))
    gain_dictionary['local']['eg_unb'].append(np.nansum(gain_mask * grav_en * unbound_mask, axis=-1))
    gain_dictionary['local']['ei'].append(np.nansum(gain_mask * int_en * gain_mask, axis=-1))
    gain_dictionary['local']['ei_unb'].append(np.nansum(gain_mask * int_en * unbound_mask, axis=-1))
    gain_dictionary['local']['nu_heat'].append(np.nansum(gain_mask * nu_heat * gain_mask, axis=-1))
    gain_dictionary['local']['nu_heat_unb'].append(np.nansum(gain_mask * nu_heat * unbound_mask, axis=-1))
    gain_dictionary['local']['tadv_unb'].append(tadv_unb)
    gain_dictionary['local']['tadv'].append(tadv_gain)
    gain_dictionary['local']['mass'].append(np.nansum(mass * gain_mask, axis=-1))
    gain_dictionary['local']['mass_unb'].append(np.nansum(mass * unbound_mask, axis=-1))

    gain_dictionary['global']['er'].append(np.nansum(gain_dictionary['local']['er'][-1]))
    gain_dictionary['global']['er_unb'].append(np.nansum(gain_dictionary['local']['er_unb'][-1]))
    gain_dictionary['global']['eg'].append(np.nansum(gain_dictionary['local']['eg'][-1]))
    gain_dictionary['global']['eg_unb'].append(np.nansum(gain_dictionary['local']['eg_unb'][-1]))
    gain_dictionary['global']['ei'].append(np.nansum(gain_dictionary['local']['ei'][-1]))
    gain_dictionary['global']['ei_unb'].append(np.nansum(gain_dictionary['local']['ei_unb'][-1]))
    gain_dictionary['global']['nu_heat'].append(np.nansum(gain_dictionary['local']['nu_heat'][-1]))
    gain_dictionary['global']['nu_heat_unb'].append(np.nansum(gain_dictionary['local']['nu_heat_unb'][-1]))
    gain_dictionary['global']['ekin'].append(gain_dictionary['global']['er'][-1])
    gain_dictionary['global']['ekin_unb'].append(gain_dictionary['global']['er_unb'][-1])

def __compute_2D_global_local_quantities(gain_dictionary: dict,
                                         simulation: Simulation,
                                         file_name: str,
                                         gain_mask: np.array[bool],
                                         unbound_mask: np.array[bool],
                                         dV: aerray,
                                         Lx: aerray,
                                         Ly: aerray,
                                         Lz: aerray) -> None:
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    rho = simulation.rho(file_name) * dV
    
    et = 0.5 * vt ** 2 * rho
    ep = 0.5 * vp ** 2 * rho

    jx = np.nansum(Lx[gain_mask])
    jy = np.nansum(Ly[gain_mask])
    jz = np.nansum(Lz[gain_mask])
    jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
    jx_unb = np.nansum(Lx[unbound_mask])
    jy_unb = np.nansum(Ly[unbound_mask])
    jz_unb = np.nansum(Lz[unbound_mask])
    jtot_unb = np.sqrt(jx_unb ** 2 + jy_unb ** 2 + jz_unb ** 2)

    gain_dictionary['global']['jx'].append(jx)
    gain_dictionary['global']['jy'].append(jy)
    gain_dictionary['global']['jz'].append(jz)
    gain_dictionary['global']['jtot'].append(jtot)
    gain_dictionary['global']['jx_unb'].append(jx_unb)
    gain_dictionary['global']['jy_unb'].append(jy_unb)
    gain_dictionary['global']['jz_unb'].append(jz_unb)
    gain_dictionary['global']['jtot_unb'].append(jtot_unb)

    gain_dictionary['local']['et'].append(np.nansum(et * gain_mask, axis=-1))
    gain_dictionary['local']['ep'].append(np.nansum(ep * gain_mask, axis=-1))
    gain_dictionary['local']['et_unb'].append(np.nansum(et * unbound_mask, axis=-1))
    gain_dictionary['local']['ep_unb'].append(np.nansum(ep * unbound_mask, axis=-1))

    gain_dictionary['global']['et'].append(np.nansum(gain_dictionary['local']['et'][-1]))
    gain_dictionary['global']['ep'].append(np.nansum(gain_dictionary['local']['ep'][-1]))
    gain_dictionary['global']['et_unb'].append(np.nansum(gain_dictionary['local']['et_unb'][-1]))
    gain_dictionary['global']['ep_unb'].append(np.nansum(gain_dictionary['local']['ep_unb'][-1]))
    gain_dictionary['global']['ekin'][-1] += gain_dictionary['global']['et'][-1] + gain_dictionary['global']['ep'][-1]
    gain_dictionary['global']['ekin_unb'][-1] += gain_dictionary['global']['et_unb'][-1] + gain_dictionary['global']['ep_unb'][-1]
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        ebr = br ** 2 * 0.5 / c.mu0 * dV
        ebt = bt ** 2 * 0.5 / c.mu0 * dV
        ebp = bp ** 2 * 0.5 / c.mu0 * dV
        gain_dictionary['local']['ebr'].append(np.nansum(ebr * gain_mask, axis=-1))
        gain_dictionary['local']['ebt'].append(np.nansum(ebt * gain_mask, axis=-1))
        gain_dictionary['local']['ebp'].append(np.nansum(ebp * gain_mask, axis=-1))
        gain_dictionary['local']['ebr_unb'].append(np.nansum(ebr * unbound_mask, axis=-1))
        gain_dictionary['local']['ebt_unb'].append(np.nansum(ebt * unbound_mask, axis=-1))
        gain_dictionary['local']['ebp_unb'].append(np.nansum(ebp * unbound_mask, axis=-1))

        gain_dictionary['global']['ebr'].append(np.nansum(gain_dictionary['local']['ebr'][-1]))
        gain_dictionary['global']['ebt'].append(np.nansum(gain_dictionary['local']['ebt'][-1]))
        gain_dictionary['global']['ebp'].append(np.nansum(gain_dictionary['local']['ebp'][-1]))
        gain_dictionary['global']['ebr_unb'].append(np.nansum(gain_dictionary['local']['ebr_unb'][-1]))
        gain_dictionary['global']['ebt_unb'].append(np.nansum(gain_dictionary['local']['ebt_unb'][-1]))
        gain_dictionary['global']['ebp_unb'].append(np.nansum(gain_dictionary['local']['ebp_unb'][-1]))

        gain_dictionary['global']['ebtot'].append(gain_dictionary['global']['ebr'][-1] +
                                                  gain_dictionary['global']['ebt'][-1] +
                                                  gain_dictionary['global']['ebp'][-1])
        gain_dictionary['global']['ebtot_unb'].append(gain_dictionary['global']['ebr_unb'][-1] +
                                                      gain_dictionary['global']['ebt_unb'][-1] +
                                                      gain_dictionary['global']['ebp_unb'][-1])