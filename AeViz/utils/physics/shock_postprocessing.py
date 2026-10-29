from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import shock_radius
from AeViz.units.constants import constants as c

def declare_shock_dictionary(dim: int,
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
    shock_dictionary = {
        'global': {
            'mass': [],
            'mass_unb': [],
            'rmax': [],
            'rmin': [],
            'ravg': [],
            'eunb': [],
            'er': [],
            'er_unb': [],
            'eg': [],
            'eg_unb': [],
            'ei': [],
            'ei_unb': [],
            'nu_heat': [],
            'nu_heat_unb': [],
            'I' : [],
            'I_unb': [],
            'ej_ye_045': [],
            'ej_ye_045' :[],
            'ej_rye_045': [],
            'ej_rye_055': [],
            },
        'local': {
            'r': [],
            'm': [],
            'm_unb': [],
            'eunb': [],
            'ei': [],
            'ei_unb': [],
            'er': [],
            'er_unb': [],
            'eg': [],
            'eg_unb': [],
            'nu_heat': [],
            'nu_heat_unb': [],
        }
    }
    if dim > 1:
        shock_dictionary['global']['et'] = []
        shock_dictionary['global']['ep'] = []
        shock_dictionary['global']['et_unb'] = []
        shock_dictionary['global']['ep_unb'] = []
        shock_dictionary['global']['jx'] = []
        shock_dictionary['global']['jy'] = []
        shock_dictionary['global']['jz'] = []
        shock_dictionary['global']['jtot'] = []
        shock_dictionary['global']['jx_unb'] = []
        shock_dictionary['global']['jy_unb'] = []
        shock_dictionary['global']['jz_unb'] = []
        shock_dictionary['global']['jtot_unb'] = []
        shock_dictionary['global']['alpha_mom'] = []
        shock_dictionary['global']['alpha_ej'] = []
        shock_dictionary['local']['et'] = []
        shock_dictionary['local']['ep'] = []
        shock_dictionary['local']['et_unb'] = []
        shock_dictionary['local']['ep_unb'] = []
        if magdim > 0:
            shock_dictionary['local']['ebr'] = []
            shock_dictionary['local']['ebt'] = []
            shock_dictionary['local']['ebp'] = []
            shock_dictionary['local']['ebr_unb'] = []
            shock_dictionary['local']['ebt_unb'] = []
            shock_dictionary['local']['ebp_unb'] = []
            shock_dictionary['global']['ebr'] = []
            shock_dictionary['global']['ebt'] = []
            shock_dictionary['global']['ebp'] = []
            shock_dictionary['global']['ebr_unb'] = []
            shock_dictionary['global']['ebt_unb'] = []
            shock_dictionary['global']['ebp_unb'] = []
    if dim > 2:
        shock_dictionary['global']['erot'] = []
        shock_dictionary['global']['eres'] = []
        shock_dictionary['global']['erot_unb'] = []
        shock_dictionary['global']['eres_unb'] = []
        shock_dictionary['global']['cmass_x'] = []
        shock_dictionary['global']['cmass_y'] = []
        shock_dictionary['global']['cmass_z'] = []
        shock_dictionary['global']['cmass_x_unb'] = []
        shock_dictionary['global']['cmass_y_unb'] = []
        shock_dictionary['global']['cmass_z_unb'] = []
        if magdim > 0:
            shock_dictionary['global']['ebpol'] = []
            shock_dictionary['global']['ebtor'] = []
            shock_dictionary['global']['ebpol_unb'] = []
            shock_dictionary['global']['ebtor_unb'] = []
            
    return shock_dictionary    

def shock_quantity_metadata() -> dict:
    """
    Metadata for shock post-processing quantities.
    """

    meta = {
        "global": {
            "mass": dict(name="shock_mass", label=r"$M_{\rm sh}$",
                         limits=[1e-3, 2.5], cmap="viridis", log=True),
            "mass_unb": dict(name="shock_unbound_mass", label=r"$M_{\rm sh,unb}$",
                             limits=[1e-5, 2.5], cmap="viridis", log=True),

            "rmax": dict(name="shock_radius_max", label=r"$R_{\rm sh,max}$",
                         limits=[5, 1e4], cmap="plasma", log=False),
            "rmin": dict(name="shock_radius_min", label=r"$R_{\rm sh,min}$",
                         limits=[5, 1e4], cmap="plasma", log=False),
            "ravg": dict(name="shock_radius_avg", label=r"$\langle R_{\rm sh}\rangle$",
                         limits=[5, 1e4], cmap="plasma", log=False),

            "eunb": dict(name="shock_unbound_energy", label=r"$E_{\rm sh,unb}$",
                         limits=[1e46, 1e52], cmap="inferno", log=True),

            "er": dict(name="shock_radial_kinetic_energy", label=r"$E_{\rm kin,r,sh}$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "er_unb": dict(name="shock_unbound_radial_kinetic_energy",
                           label=r"$E_{\rm kin,r,sh,unb}$",
                           limits=[1e46, 1e52], cmap="inferno", log=True),

            "eg": dict(name="shock_gravitational_energy", label=r"$E_{\rm grav,sh}$",
                       limits=[-5e53, -1e45], cmap="magma", log=False),
            "eg_unb": dict(name="shock_unbound_gravitational_energy",
                           label=r"$E_{\rm grav,sh,unb}$",
                           limits=[-5e53, -1e45], cmap="magma", log=False),

            "ei": dict(name="shock_internal_energy", label=r"$E_{\rm int,sh}$",
                       limits=[1e46, 5e53], cmap="hot", log=True),
            "ei_unb": dict(name="shock_unbound_internal_energy",
                           label=r"$E_{\rm int,sh,unb}$",
                           limits=[1e46, 5e53], cmap="hot", log=True),

            "nu_heat": dict(name="shock_neutrino_heating_energy",
                            label=r"$E_{\nu,\rm heat,sh}$",
                            limits=[1e46, 1e52], cmap="YlOrRd", log=True),
            "nu_heat_unb": dict(name="shock_unbound_neutrino_heating_energy",
                                label=r"$E_{\nu,\rm heat,sh,unb}$",
                                limits=[1e46, 1e52], cmap="YlOrRd", log=True),

            "I": dict(name="shock_moment_of_inertia", label=r"$I_{\rm sh}$",
                      limits=[1e43, 1e50], cmap="Greens", log=True),
            "I_unb": dict(name="shock_unbound_moment_of_inertia",
                          label=r"$I_{\rm sh,unb}$",
                          limits=[1e43, 1e50], cmap="Greens", log=True),

            "et": dict(name="shock_theta_kinetic_energy",
                       label=r"$E_{\rm kin,\theta,sh}$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "ep": dict(name="shock_phi_kinetic_energy",
                       label=r"$E_{\rm kin,\phi,sh}$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "et_unb": dict(name="shock_unbound_theta_kinetic_energy",
                           label=r"$E_{\rm kin,\theta,sh,unb}$",
                           limits=[1e46, 1e52], cmap="inferno", log=True),
            "ep_unb": dict(name="shock_unbound_phi_kinetic_energy",
                           label=r"$E_{\rm kin,\phi,sh,unb}$",
                           limits=[1e46, 1e52], cmap="inferno", log=True),

            "jx": dict(name="shock_jx", label=r"$J_{x,\rm sh}$",
                       limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jy": dict(name="shock_jy", label=r"$J_{y,\rm sh}$",
                       limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jz": dict(name="shock_jz", label=r"$J_{z,\rm sh}$",
                       limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jtot": dict(name="shock_jtot", label=r"$|\mathbf{J}_{\rm sh}|$",
                         limits=[1e46, 1e50], cmap="viridis", log=True),

            "jx_unb": dict(name="shock_jx_unb", label=r"$J_{x,\rm sh,unb}$",
                           limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jy_unb": dict(name="shock_jy_unb", label=r"$J_{y,\rm sh,unb}$",
                           limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jz_unb": dict(name="shock_jz_unb", label=r"$J_{z,\rm sh,unb}$",
                           limits=[-5e49, 5e49], cmap="coolwarm", log=False),
            "jtot_unb": dict(name="shock_jtot_unb",
                             label=r"$|\mathbf{J}_{\rm sh,unb}|$",
                             limits=[1e46, 1e50], cmap="viridis", log=True),

            "ebr": dict(name="shock_magnetic_energy_radial",
                        label=r"$E_{\rm mag,r,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebt": dict(name="shock_magnetic_energy_theta",
                        label=r"$E_{\rm mag,\theta,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebp": dict(name="shock_magnetic_energy_phi",
                        label=r"$E_{\rm mag,\phi,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebr_unb": dict(name="shock_unbound_magnetic_energy_radial",
                            label=r"$E_{\rm mag,r,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebt_unb": dict(name="shock_unbound_magnetic_energy_theta",
                            label=r"$E_{\rm mag,\theta,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebp_unb": dict(name="shock_unbound_magnetic_energy_phi",
                            label=r"$E_{\rm mag,\phi,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),

            "erot": dict(name="shock_rotational_energy",
                         label=r"$E_{\rm rot,sh}$",
                         limits=[1e46, 1e53], cmap="Oranges", log=True),
            "eres": dict(name="shock_turbulent_energy",
                         label=r"$E_{\rm turb,sh}$",
                         limits=[1e46, 1e53], cmap="YlOrBr", log=True),
            "erot_unb": dict(name="shock_unbound_rotational_energy",
                             label=r"$E_{\rm rot,sh,unb}$",
                             limits=[1e46, 1e53], cmap="Oranges", log=True),
            "eres_unb": dict(name="shock_unbound_turbulent_energy",
                             label=r"$E_{\rm turb,sh,unb}$",
                             limits=[1e46, 1e53], cmap="YlOrBr", log=True),

            "cmass_x": dict(name="shock_center_of_mass_x",
                            label=r"$x_{\rm COM,sh}$",
                            limits=[-1e4, 1e4], cmap="coolwarm", log=False),
            "cmass_y": dict(name="shock_center_of_mass_y",
                            label=r"$y_{\rm COM,sh}$",
                            limits=[-1e4, 1e4], cmap="coolwarm", log=False),
            "cmass_z": dict(name="shock_center_of_mass_z",
                            label=r"$z_{\rm COM,sh}$",
                            limits=[-1e4, 1e4], cmap="coolwarm", log=False),
            "cmass_x_unb": dict(name="shock_unbound_center_of_mass_x",
                                label=r"$x_{\rm COM,sh,unb}$",
                                limits=[-1e4, 1e4], cmap="coolwarm", log=False),
            "cmass_y_unb": dict(name="shock_unbound_center_of_mass_y",
                                label=r"$y_{\rm COM,sh,unb}$",
                                limits=[-1e4, 1e4], cmap="coolwarm", log=False),
            "cmass_z_unb": dict(name="shock_unbound_center_of_mass_z",
                                label=r"$z_{\rm COM,sh,unb}$",
                                limits=[-1e4, 1e4], cmap="coolwarm", log=False),

            "ebpol": dict(name="shock_poloidal_magnetic_energy",
                          label=r"$E_{\rm mag,pol,sh}$",
                          limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebtor": dict(name="shock_toroidal_magnetic_energy",
                          label=r"$E_{\rm mag,tor,sh}$",
                          limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebpol_unb": dict(name="shock_unbound_poloidal_magnetic_energy",
                              label=r"$E_{\rm mag,pol,sh,unb}$",
                              limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebtor_unb": dict(name="shock_unbound_toroidal_magnetic_energy",
                              label=r"$E_{\rm mag,tor,sh,unb}$",
                              limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ej_ye_045": dict(name="ye_ejection_tail_0.45",
                              label=r"$M_{Y_e\leq0.45}$",
                              limits=[1e-5, 1], cmap="PuBuGn", log=True),
            "ej_ye_055": dict(name="ye_ejection_tail_0.45",
                              label=r"$M_{Y_e\geq0.55}$",
                              limits=[1e-5, 1], cmap="PuBuGn", log=True),
            "ej_rye_045": dict(name="centroid_ye_ejection_tail_0.45",
                              label=r"$\langle r\rangle_{Y_e\leq0.45}$",
                              limits=[1e4, 1e8], cmap="PuBuGn", log=True),
            "ej_rye_055": dict(name="centroid_ye_ejection_tail_0.55",
                              label=r"$\langle r\rangle_{Y_e\geq0.55}$",
                              limits=[1e4, 1e8], cmap="PuBuGn", log=True),
                              
        },

        "local": {
            "r": dict(name="shock_radius", label=r"$r_{\rm sh}$",
                      limits=[5, 1e4], cmap="plasma", log=False),
            "m": dict(name="shock_enclosed_mass", label=r"$M_{\rm sh}(r)$",
                      limits=[1e-3, 2.5], cmap="viridis", log=True),
            "m_unb": dict(name="shock_unbound_enclosed_mass",
                          label=r"$M_{\rm sh,unb}(r)$",
                          limits=[1e-5, 2.5], cmap="viridis", log=True),

            "eunb": dict(name="shock_unbound_energy_profile",
                         label=r"$E_{\rm sh,unb}(r)$",
                         limits=[1e46, 1e52], cmap="inferno", log=True),
            "ei": dict(name="shock_internal_energy_profile",
                       label=r"$E_{\rm int,sh}(r)$",
                       limits=[1e46, 5e53], cmap="hot", log=True),
            "ei_unb": dict(name="shock_unbound_internal_energy_profile",
                           label=r"$E_{\rm int,sh,unb}(r)$",
                           limits=[1e46, 5e53], cmap="hot", log=True),
            "er": dict(name="shock_radial_kinetic_energy_profile",
                       label=r"$E_{\rm kin,r,sh}(r)$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "er_unb": dict(name="shock_unbound_radial_kinetic_energy_profile",
                           label=r"$E_{\rm kin,r,sh,unb}(r)$",
                           limits=[1e46, 1e52], cmap="inferno", log=True),
            "eg": dict(name="shock_gravitational_energy_profile",
                       label=r"$E_{\rm grav,sh}(r)$",
                       limits=[-5e53, -1e45], cmap="magma", log=False),
            "eg_unb": dict(name="shock_unbound_gravitational_energy_profile",
                           label=r"$E_{\rm grav,sh,unb}(r)$",
                           limits=[-5e53, -1e45], cmap="magma", log=False),
            "nu_heat": dict(name="shock_neutrino_heating_energy_profile",
                            label=r"$E_{\nu,\rm heat,sh}(r)$",
                            limits=[1e46, 1e52], cmap="YlOrRd", log=True),
            "nu_heat_unb": dict(name="shock_unbound_neutrino_heating_energy_profile",
                                label=r"$E_{\nu,\rm heat,sh,unb}(r)$",
                                limits=[1e46, 1e52], cmap="YlOrRd", log=True),
            "et": dict(name="shock_theta_kinetic_energy_profile",
                       label=r"$E_{\rm kin,\theta,sh}(r)$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "ep": dict(name="shock_phi_kinetic_energy_profile",
                       label=r"$E_{\rm kin,\phi,sh}(r)$",
                       limits=[1e46, 1e52], cmap="inferno", log=True),
            "ebr": dict(name="shock_magnetic_energy_radial_prof",
                        label=r"$E_{\rm mag,r,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebt": dict(name="shock_magnetic_energy_theta_prof",
                        label=r"$E_{\rm mag,\theta,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebp": dict(name="shock_magnetic_energy_phi_prof",
                        label=r"$E_{\rm mag,\phi,sh}$",
                        limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebr_unb": dict(name="shock_unbound_magnetic_energy_radial_prof",
                            label=r"$E_{\rm mag,r,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebt_unb": dict(name="shock_unbound_magnetic_energy_theta_prof",
                            label=r"$E_{\rm mag,\theta,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "ebp_unb": dict(name="shock_unbound_magnetic_energy_phi_prof",
                            label=r"$E_{\rm mag,\phi,sh,unb}$",
                            limits=[1e42, 1e51], cmap="PuBuGn", log=True),
            "alpha_ej": dict(name="momentum_asymmetry",
                              label=r"$\alpha_{\rm mom}$",
                              limits=[-1, 1], cmap="PuBuGn", log=False),
            "alpha_mom": dict(name="mass_asymmetry",
                              label=r"$\alpha_{\rm mom}$",
                              limits=[-1, 1], cmap="PuBuGn", log=False),
        },
    }

    return meta

def compute_shock_postprocessing(shock_dictionary: dict,
                               simulation: Simulation,
                               file_name: str,
                               dV: aerray,
                               dOmega: aerray,
                               Lx: aerray | None,
                               Ly: aerray | None,
                               Lz: aerray | None,
                               pns_radius: aerray,
                               innercore_radius: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    shock_dictionary : dict
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
    shockr = shock_radius(simulation, file_name,
                          pns_radius, innercore_radius)
    r = simulation.cell.radius(simulation.ghost)
    if simulation.dim == 1:
        shockmask = (r <= shockr)
        idx_shock = np.argmax(r >= shockr)
    else:
        while r.ndim <= shockr.ndim:
            r = r[None, :]
        shockmask = (r <= shockr)
        idx_shock = np.argmax(r >= shockr[..., None], axis=-1)
    dmass = (simulation.rho(file_name) * dV).to(u.M_sun)
    mhd_ene = simulation.MHD_energy(file_name) + \
            2 * simulation.gravitational_energy(file_name)
    unbound_mask = (mhd_ene > 0) & (simulation.cell.radius(simulation.ghost) < 1e11) & \
        (simulation.radial_velocity(file_name) > 0)
    __compute_1D_global_local_quantities(shock_dictionary,
                                         simulation,
                                         file_name,
                                         dV,
                                         dOmega,
                                         shockmask,
                                         unbound_mask,
                                         idx_shock,
                                         shockr,
                                         dmass)
    if simulation.dim > 1:
        __compute_2D_global_local_quantities(shock_dictionary,
                                            simulation,
                                            file_name,
                                            shockmask,
                                            unbound_mask,
                                            idx_shock,
                                            dV,
                                            Lx,
                                            Ly,
                                            Lz,
                                            dmass)
    if simulation.dim == 3:
        __compute_3D_global_local_quantities(shock_dictionary,
                                            simulation,
                                            file_name,
                                            shockmask,
                                            unbound_mask,
                                            idx_shock,
                                            dmass,
                                            dV)

def __compute_1D_global_local_quantities(shock_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       dV: aerray,
                                       dOmega: aerray,
                                       shock_mask: np.array[bool],
                                       unbound_mask: np.array[bool],
                                       indx_shock: np.array[int],
                                       shock_radius: aerray,
                                       dmass:aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    shock_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    dV : aerray
        _description_
    dOmega : aerray
        _description_
    shock_mask : np.array[bool]
        _description_
    indx_shock : np.array[int]
        _description_
    shock_radius : aerray
        _description_
    I : aerray
        _description_
    dmass : aerray
        _description_
    """
    grav_en = simulation.gravitational_energy(file_name) * dV
    int_en = simulation.internal_energy(file_name) * dV
    vr = simulation.radial_velocity(file_name)
    mhd_ene = simulation.MHD_energy(file_name)
    nu_heat = simulation.nu_heat(file_name)
    ye = simulation(file_name)

    ## add the global variables
    shock_dictionary['global']['mass'].append(np.nansum(dmass[shock_mask]))
    shock_dictionary['global']['mass_unb'].append(np.nansum(dmass[unbound_mask]))
    shock_dictionary['global']['rmax'].append(np.nanmax(shock_radius))
    shock_dictionary['global']['rmin'].append(np.nanmin(shock_radius))
    shock_dictionary['global']['ravg'].append(simulation.grid.radii_angular_average(shock_radius))
    shock_dictionary['global']['eunb'].append(np.nansum(
                                              mhd_ene[unbound_mask] * 
                                              dV[unbound_mask]).to(u.erg))
    shock_dictionary['global']['er'].append(np.nansum(0.5 * dmass[shock_mask] * 
                                                    vr[shock_mask] ** 2).to(u.erg))
    shock_dictionary['global']['er_unb'].append(np.nansum(0.5 * dmass[unbound_mask] * 
                                                        vr[unbound_mask] ** 2).to(u.erg))
    shock_dictionary['global']['eg'].append(np.nansum(grav_en[shock_mask]))
    shock_dictionary['global']['eg_unb'].append(np.nansum(grav_en[unbound_mask]))
    shock_dictionary['global']['ei'].append(np.nansum(int_en[shock_mask]))
    shock_dictionary['global']['ei_unb'].append(np.nansum(int_en[unbound_mask]))
    shock_dictionary['global']['nu_heat'].append(np.nansum(nu_heat[shock_mask] *
                                                           dV[shock_mask]))
    shock_dictionary['global']['nu_heat_unb'].append(np.nansum(nu_heat[unbound_mask] *
                                                               dV[unbound_mask]))
    ## ejecta
    shock_dictionary['global']['ej_ye_045'].append(np.nansum(dmass[unbound_mask & (ye <= 0.45)]))
    shock_dictionary['global']['ej_ye_055'].append(np.nansum(dmass[unbound_mask & (ye >= 0.55)]))
    centroid = dmass * simulation.cell.radius(simulation.ghost)
    shock_dictionary['global']['ej_rye_045'].append(np.nansum(centroid[unbound_mask & (ye <= 0.45)]) / 
                                                    shock_dictionary['global']['ej_ye_045'][-1])
    shock_dictionary['global']['ej_rye_055'].append(np.nansum(centroid[unbound_mask & (ye >= 0.55)]) / 
                                                    shock_dictionary['global']['ej_ye_055'][-1])
    ## add variables at the surface
    shock_dictionary['local']['r'].append(shock_radius)
    if simulation.dim == 1:
        shock_dictionary['local']['m'].append(shock_dictionary['global']['mass'][-1])
        shock_dictionary['local']['m_unb'].append(shock_dictionary['global']['mass_unb'][-1])
        shock_dictionary['local']['eunb'].append(shock_dictionary['global']['eunb'][-1])
        shock_dictionary['local']['ei'].append(shock_dictionary['global']['ei'][-1])
        shock_dictionary['local']['ei_unb'].append(shock_dictionary['global']['ei_unb'][-1])
        shock_dictionary['local']['er'].append(shock_dictionary['global']['er'][-1])
        shock_dictionary['local']['er_unb'].append(shock_dictionary['global']['er_unb'][-1])
        shock_dictionary['local']['eg'].append(shock_dictionary['global']['eg'][-1])
        shock_dictionary['local']['eg_unb'].append(shock_dictionary['global']['eg_unb'][-1])
        shock_dictionary['local']['nu_heat'].append(shock_dictionary['global']['nu_heat'][-1])
        shock_dictionary['local']['nu_heat_unb'].append(shock_dictionary['global']['nu_heat_unb'][-1])
        rperp2 = simulation.cell.radius(simulation.ghost) ** 2
        shock_dictionary['global']['I'].append(np.nansum((dmass.to(u.g) * 
                                                        rperp2)[shock_mask]))
        shock_dictionary['global']['I_unb'].append(np.nansum((dmass.to(u.g) * 
                                                                rperp2)[unbound_mask]))
    elif simulation.dim == 2:
        itheta = np.arange(vr.shape[0])
        shock_dictionary['local']['m'].append(np.nancumsum(dmass,
                                                         axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['m_unb'].append(np.nansum(np.where(unbound_mask,
                                                                        dmass,
                                                                        0), axis=-1))
        shock_dictionary['local']['eunb'].append(np.nansum(np.where(unbound_mask,
                                                                        mhd_ene * dV,
                                                                        0), axis=-1))
        shock_dictionary['local']['ei'].append(np.nancumsum(int_en,
                                                         axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['ei_unb'].append(np.nansum(np.where(unbound_mask,
                                                                        int_en, 0), axis=-1))
        shock_dictionary['local']['er'].append(np.nancumsum(0.5 * dmass * 
                                                            vr ** 2,
                                                         axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['er_unb'].append(np.nansum(np.where(unbound_mask,
                                                                        0.5 * dmass[unbound_mask] * 
                                                                        vr[unbound_mask] ** 2,
                                                                        0), axis=-1))
        shock_dictionary['local']['eg'].append(np.nancumsum(grav_en,
                                                         axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['eg_unb'].append(np.nansum(np.where(unbound_mask,
                                                                         grav_en,
                                                                        0), axis=-1))
        shock_dictionary['local']['nu_heat'].append(np.nancumsum(nu_heat * dV,
                                                         axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['nu_heat_unb'].append(np.nansum(np.where(unbound_mask,
                                                                         nu_heat * dV,
                                                                        0), axis=-1))
    elif simulation.dim == 3:
        itheta = np.arange(vr.shape[1])
        iphi = np.arange(vr.shape[0])
        shock_dictionary['local']['m'].append(np.nancumsum(dmass,
                                                            axis=-1)[iphi,
                                                                    itheta,
                                                                    indx_shock])
        shock_dictionary['local']['m_unb'].append(np.nansum(np.where(unbound_mask,
                                                                                dmass,
                                                                                0), axis=-1))
        shock_dictionary['local']['eunb'].append(np.nansum(np.where(unbound_mask,
                                                                        mhd_ene * dV,
                                                                        0), axis=-1))
        shock_dictionary['local']['ei'].append(np.nancumsum(int_en,
                                                            axis=-1)[iphi,
                                                                     itheta,
                                                                     indx_shock])
        shock_dictionary['local']['ei_unb'].append(np.nansum(np.where(unbound_mask,
                                                                        int_en, 0), axis=-1))
        shock_dictionary['local']['er'].append(np.nancumsum(0.5 * dmass * 
                                                            vr ** 2,
                                                            axis=-1)[iphi, itheta,
                                                            indx_shock])
        shock_dictionary['local']['er_unb'].append(np.nansum(np.where(unbound_mask,
                                                                        0.5 * dmass[unbound_mask] * 
                                                                        vr[unbound_mask] ** 2,
                                                                        0), axis=-1))
        shock_dictionary['local']['eg'].append(np.nancumsum(grav_en,
                                                            axis=-1)[iphi, itheta,
                                                            indx_shock])
        shock_dictionary['local']['eg_unb'].append(np.nansum(np.where(unbound_mask,
                                                                            grav_en,
                                                                        0), axis=-1))
        shock_dictionary['local']['nu_heat'].append(np.nancumsum(nu_heat * dV,
                                                            axis=-1)[iphi, itheta,
                                                            indx_shock])
        shock_dictionary['local']['nu_heat_unb'].append(np.nansum(np.where(unbound_mask,
                                                                            nu_heat * dV,
                                                                        0), axis=-1))

def __compute_2D_global_local_quantities(shock_dictionary: dict,
                                       simulation: Simulation,
                                       file_name: str,
                                       shock_mask: np.array[bool],
                                       unbound_mask: np.array[bool],
                                       indx_shock: np.array[int],
                                       dV: aerray,
                                       Lx: aerray,
                                       Ly: aerray,
                                       Lz: aerray,
                                       dmass: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    shock_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    shock_mask : np.array[bool]
        _description_
    indx_shock : np.array[int]
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
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    vtot = np.sqrt(vr ** 2 + vt ** 2 + vp ** 2)
    rr = vr / vtot
    rt = vt / vtot
    rp = vp / vtot
    alpha_ej = np.sqrt(np.nansum(dmass[unbound_mask] * rr[unbound_mask]) ** 2 + 
                       np.nansum(dmass[unbound_mask] * rt[unbound_mask]) ** 2 +
                       np.nansum(dmass[unbound_mask] * rp[unbound_mask]) ** 2) / \
                shock_dictionary['local']['m_unb'][-1]
    Pej = np.sqrt(np.nansum((dmass[unbound_mask] * vr[unbound_mask])) ** 2 +
                  np.nansum((dmass[unbound_mask] * vt[unbound_mask])) ** 2 +
                  np.nansum((dmass[unbound_mask] * vp[unbound_mask])) ** 2) / \
                  np.nansum(dmass[unbound_mask] * vtot[unbound_mask])
    
    eth = (0.5 * vt ** 2 * dV).to(u.erg)
    eph = (0.5 * vp ** 2 * dV).to(u.erg)
    jx = np.nansum(Lx[shock_mask])
    jy = np.nansum(Ly[shock_mask])
    jz = np.nansum(Lz[shock_mask])
    jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
    jx_unb = np.nansum(Lx[unbound_mask])
    jy_unb = np.nansum(Ly[unbound_mask])
    jz_unb = np.nansum(Lz[unbound_mask])
    jtot_unb = np.sqrt(jx_unb ** 2 + jy_unb ** 2 + jz_unb ** 2)
    shock_dictionary['global']['et'].append(np.nansum(eth[shock_mask]))
    shock_dictionary['global']['ep'].append(np.nansum(eph[shock_mask]))
    shock_dictionary['global']['et_unb'].append(np.nansum(eth[unbound_mask]))
    shock_dictionary['global']['ep_unb'].append(np.nansum(eph[unbound_mask]))
    shock_dictionary['global']['jx'].append(jx)
    shock_dictionary['global']['jy'].append(jy)
    shock_dictionary['global']['jz'].append(jz)
    shock_dictionary['global']['jtot'].append(jtot)
    shock_dictionary['global']['jx_unb'].append(jx_unb)
    shock_dictionary['global']['jy_unb'].append(jy_unb)
    shock_dictionary['global']['jz_unb'].append(jz_unb)
    shock_dictionary['global']['jtot_unb'].append(jtot_unb)
    shock_dictionary['global']['alpha_ej'].append(alpha_ej)
    shock_dictionary['global']['alpha_mom'].append(Pej)
    shock_dictionary['local']['et_unb'].append(np.nansum(np.where(unbound_mask,
                                                                  eth, 0), axis=-1))
    shock_dictionary['local']['ep_unb'].append(np.nansum(np.where(unbound_mask,
                                                                  eph, 0), axis=-1))
    if simulation.dim == 2:
        itheta = np.arange(eth.shape[0])
        rperp2 = (simulation.cell.radius(simulation.ghost)[None, :] *
                          np.sin(simulation.cell.theta(simulation.ghost))[:, None]) ** 2
        shock_dictionary['global']['I'].append(np.nansum((dmass.to(u.g) * 
                                                        rperp2)[shock_mask]))
        shock_dictionary['global']['I_unb'].append(np.nansum((dmass.to(u.g) * 
                                                                rperp2)[unbound_mask]))
        shock_dictionary['local']['et'].append(np.nancumsum(eth,axis=-1)[itheta,
                                                            indx_shock])
        shock_dictionary['local']['ep'].append(np.nancumsum(eph,axis=-1)[itheta,
                                                            indx_shock])

    elif simulation.dim == 3:
        iphi = np.arange(eth.shape[0])[:, None]
        itheta = np.arange(eth.shape[1])[None, :]
        shock_dictionary['local']['et'].append(np.nancumsum(eth,axis=-1)[iphi,
                                                                         itheta,
                                                                    indx_shock])
        shock_dictionary['local']['ep'].append(np.nancumsum(eph,axis=-1)[iphi,
                                                                         itheta,
                                                                    indx_shock])
    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        ebr = br ** 2 * 0.5 / c.mu0 * dV
        ebt = bt ** 2 * 0.5 / c.mu0 * dV
        ebp = bp ** 2 * 0.5 / c.mu0 * dV
        shock_dictionary['global']['ebr'].append(np.nansum(ebr[shock_mask]))
        shock_dictionary['global']['ebt'].append(np.nansum(ebt[shock_mask]))
        shock_dictionary['global']['ebp'].append(np.nansum(ebp[shock_mask]))
        shock_dictionary['global']['ebr_unb'].append(np.nansum(ebr[unbound_mask]))
        shock_dictionary['global']['ebt_unb'].append(np.nansum(ebt[unbound_mask]))
        shock_dictionary['global']['ebp_unb'].append(np.nansum(ebp[unbound_mask]))
        shock_dictionary['local']['ebr_unb'].append(np.nansum(ebr[unbound_mask], axis=-1))
        shock_dictionary['local']['ebt_unb'].append(np.nansum(ebt[unbound_mask], axis=-1))
        shock_dictionary['local']['ebp_unb'].append(np.nansum(ebp[unbound_mask], axis=-1))
        if simulation.dim == 2:
            shock_dictionary['local']['ebr'].append(np.nancumsum(ebr, axis=-1)[itheta,
                                                                    indx_shock])
            shock_dictionary['local']['ebt'].append(np.nancumsum(ebt, axis=-1)[itheta,
                                                                    indx_shock])
            shock_dictionary['local']['ebp'].append(np.nancumsum(ebp, axis=-1)[itheta,
                                                                    indx_shock])
        elif simulation.dim == 3:
            shock_dictionary['local']['ebr'].append(np.nancumsum(ebr,
                                                                axis=-1)[iphi,
                                                                    itheta,
                                                                    indx_shock])
            shock_dictionary['local']['ebt'].append(np.nancumsum(ebt,
                                                                axis=-1)[iphi,
                                                                    itheta,
                                                                    indx_shock])
            shock_dictionary['local']['ebp'].append(np.nancumsum(ebp,
                                                                axis=-1)[iphi,
                                                                    itheta,
                                                                    indx_shock])

def __compute_3D_global_local_quantities(shock_dictionary: dict,
                                         simulation: Simulation,
                                         file_name: str,
                                         shock_mask: np.array[bool],
                                         unbound_mask: np.array[bool],
                                         indx_shock: np.array[int],
                                         dmass: aerray,
                                         dV: aerray) -> None:
    """
    _summary_

    Parameters
    ----------
    shock_dictionary : dict
        _description_
    simulation : Simulation
        _description_
    file_name : str
        _description_
    shock_mask : np.array[bool]
        _description_
    indx_shock : np.array[int]
        _description_
    dmass : aerray
        _description_
    dV : aerray
        _description_
    """
    jtot = shock_dictionary['global']['jtot'][-1].value
    nx = shock_dictionary['global']['jx'][-1] / jtot
    ny = shock_dictionary['global']['jy'][-1] / jtot
    nz = shock_dictionary['global']['jz'][-1] / jtot
    n = np.array([nx.value, ny.value, nz.value])
    if np.isnan(n).any():
        n = np.array([0, 0., 1.])
    jtot_unb = shock_dictionary['global']['jtot_unb'][-1].value
    nx_unb = shock_dictionary['global']['jx_unb'][-1] / jtot_unb
    ny_unb = shock_dictionary['global']['jy_unb'][-1] / jtot_unb
    nz_unb = shock_dictionary['global']['jz_unb'][-1] / jtot_unb
    n_unb = np.array([nx_unb.value, ny_unb.value, nz_unb.value])
    if np.isnan(n_unb).any():
        n_unb = np.array([0, 0., 1.])
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    iphi = np.arange(vr.shape[0])[:, None]
    itheta = np.arange(vr.shape[1])[None, :]
    X = simulation.grid.cartesian_grid()
    x, y, z = X[0] * u.cm, X[1] * u.cm, X[2] * u.cm
    shock_dictionary['global']['cmass_x'].append(np.nansum(x[shock_mask] * dmass[shock_mask]) / 
                                               shock_dictionary['global']['mass'][-1])
    shock_dictionary['global']['cmass_y'].append(np.nansum(y[shock_mask] * dmass[shock_mask]) / 
                                               shock_dictionary['global']['mass'][-1])
    shock_dictionary['global']['cmass_z'].append(np.nansum(z[shock_mask] * dmass[shock_mask]) / 
                                               shock_dictionary['global']['mass'][-1])
    shock_dictionary['global']['cmass_x_unb'].append(np.nansum(x[unbound_mask] * dmass[unbound_mask]) / 
                                                shock_dictionary['global']['mass_unb'][-1])
    shock_dictionary['global']['cmass_y_unb'].append(np.nansum(y[unbound_mask] * dmass[unbound_mask]) / 
                                                shock_dictionary['global']['mass_unb'][-1])
    shock_dictionary['global']['cmass_z_unb'].append(np.nansum(z[unbound_mask] * dmass[unbound_mask]) / 
                                                shock_dictionary['global']['mass_unb'][-1])
    if n[2] >= 0.995:
        vrot = vp
        vturb = vt ** 2 + vr ** 2
        mask = shock_mask
        R = (simulation.cell.radius(simulation.ghost)[None, None, :] * 
                 np.sin(simulation.cell.theta(simulation.ghost))[None, :, None])
    else:
        X_perp = X - np.tensordot(n, X, axes=(0, 0))[None, ...] * \
            n[:, None, None, None]
        R = np.linalg.norm(X_perp, axis=0) + 1e-12
        mask = ((R > simulation.grid.radius[0].value) & (shock_mask))
        e_R = X_perp / (R)
        e_rot = np.cross(n[:, None, None, None], e_R, axis=0)
        vel_unit = vr.unit
        vx, vy, vz = simulation.grid.spherical_to_cartesian(
            vr, vt, vp
        )
        v = np.array([vx.value, vy.value, vz.value])
        vrot = np.nansum(v * e_rot, axis=0)
        vturb = np.nansum((v - vrot * e_rot) ** 2, axis=0) * vel_unit ** 2
        vrot = vrot * vel_unit

    if n_unb[2] >= 0.995:
        vrot_unb = vp
        vturb_unb = vt ** 2 + vr ** 2
        unb_mask = unbound_mask
        R_unb = (simulation.cell.radius(simulation.ghost)[None, None, :] * 
                    np.sin(simulation.cell.theta(simulation.ghost))[None, :, None])
    else:
        X_perp = X - np.tensordot(n, X, axes=(0, 0))[None, ...] * \
            n[:, None, None, None]
        R_unb = np.linalg.norm(X_perp, axis=0) + 1e-12
        unb_mask = ((R > simulation.grid.radius[0].value) & (unbound_mask))
        e_R = X_perp / (R_unb)
        e_rot_unb = np.cross(n[:, None, None, None], e_R, axis=0)
        vel_unit = vr.unit
        vx, vy, vz = simulation.grid.spherical_to_cartesian(
            vr, vt, vp
        )
        v = np.array([vx.value, vy.value, vz.value])
        vrot_unb = np.nansum(v * e_rot_unb, axis=0)
        vturb_unb = np.nansum((v - vrot_unb * e_rot_unb) ** 2, axis=0) * vel_unit ** 2
        vrot_unb = vrot_unb * vel_unit
    shock_dictionary['global']['I'].append(np.nansum((dmass.to(u.g) * 
                                                    R ** 2)[shock_mask]))
    shock_dictionary['global']['erot'].append(np.nansum(0.5 * (vrot ** 2 * 
                                                             dmass)[mask]).to(u.erg))
    shock_dictionary['global']['eres'].append(np.nansum(0.5 * (vturb *
                                                             dmass)[mask]).to(u.erg))
    shock_dictionary['global']['I_unb'].append(np.nansum((dmass.to(u.g) * 
                                                    R ** 2)[unbound_mask]))
    shock_dictionary['global']['erot_unb'].append(np.nansum(0.5 * (vrot_unb ** 2 * 
                                                             dmass)[unb_mask]).to(u.erg))
    shock_dictionary['global']['eres_unb'].append(np.nansum(0.5 * (vturb_unb *
                                                             dmass)[unb_mask]).to(u.erg))
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
        shock_dictionary['global']['ebpol'].append(np.nansum(0.5 * (bpol ** 2 * 
                                                                     dV)[mask]) 
                                                 / c.mu0)
        shock_dictionary['global']['ebtor'].append(np.nansum(0.5 * (btor ** 2 * 
                                                                    dV)[mask]) 
                                                         / c.mu0)
        if n_unb[2] >= 0.995:
            bpol = np.sqrt(br ** 2 + bt ** 2)
            btor = bp
        else:
            bunit = br.unit
            bx, by, bz = simulation.grid.spherical_to_cartesian(
            br, bt, bp
            )
            b = np.array([bx.value, by.value, bz.value])
            btor = np.nansum(b * e_rot_unb, axis=0)
            bpol = np.sqrt(np.nansum((b - btor * e_rot_unb) ** 2, axis=0)) * bunit
            btor = btor * bunit
        shock_dictionary['global']['ebpol_unb'].append(np.nansum(0.5 * (bpol ** 2 * 
                                                                     dV)[unb_mask]) 
                                                 / c.mu0)
        shock_dictionary['global']['ebtor_unb'].append(np.nansum(0.5 * (btor ** 2 * 
                                                                    dV)[unb_mask]) 
                                                         / c.mu0)