from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.units.constants import constants as c
from AeViz.utils.math_utils import function_average, IDL_derivative
from AeViz.utils.files.string_utils import merge_strings

def declare_profiles_dictionary(dim: int,
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
    profiles_dictionary = {
        'inertia_profile': [],
        'radial_velocity_profile': [],
        'rho_profile': [],
        'temperature_profile': [],
        'gas_pressure_profile': [],
        'entropy_profile': [],
        'Ye_profile': [],
        'soundspeed_profile': [],
        'internal_energy_profile': [],
        'gravitational_potential_profile': [],
        'convective_flux_profile': [],
    }
    if dim > 1:
        profiles_dictionary['theta_velocity_profile'] = []
        profiles_dictionary['phi_velocity_profile'] = []
        profiles_dictionary['radial_velocity_rms_profile'] = []
        profiles_dictionary['theta_velocity_rms_profile'] = []
        profiles_dictionary['phi_velocity_rms_profile'] = []
        profiles_dictionary['angular_momentum_x_profile'] = []
        profiles_dictionary['angular_momentum_y_profile'] = []
        profiles_dictionary['angular_momentum_z_profile'] = []
        profiles_dictionary['angular_momentum_tot_profile'] = []
        profiles_dictionary['Solberg_Høiland_profile'] = []
        profiles_dictionary['omega_profile'] = []
        profiles_dictionary['Rossby_number'] = []
        profiles_dictionary['BV_frequency_profile'] = []
        profiles_dictionary['epicyclic_frequency_profile'] = []
        if magdim > 1:
            profiles_dictionary['magnetic_field_radial_profile'] = []
            profiles_dictionary['magnetic_field_theta_profile'] = []
            profiles_dictionary['magnetic_field_phi_profile'] = []
            profiles_dictionary['magnetic_field_total_profile'] = []
            
    return profiles_dictionary

def compute_profiles_postprocessing(profiles_dictionary: dict,
                                    simulation: Simulation,
                                    file_name: str,
                                    dV: aerray,
                                    dOmega: aerray,
                                    jx: aerray | None,
                                    jy: aerray | None,
                                    jz: aerray | None,
                                    jtot: aerray | None,
                                    I: aerray) -> None:
    __compute_1D_global_local_quantities(profiles_dictionary,
                                         simulation,
                                         file_name,
                                         dOmega,
                                         I)
    if simulation.dim > 1:
        __compute_2D_global_local_quantities(profiles_dictionary,
                                             simulation,
                                             file_name,
                                             dOmega,
                                             jx,
                                             jy,
                                             jz,
                                             jtot)
    if simulation.dim > 2:
        __compute_3D_global_local_quantities()


def __compute_1D_global_local_quantities(profiles_dictionary: dict,
                                        simulation: Simulation,
                                        file_name: str,
                                        dOmega: aerray,
                                        I: aerray) -> None:
    
    I_prof = np.nansum(I, axis=range(simulation.dim - 1)).to(u.g / u.cm ** 2)
    I_prof.set(name='inertia_moment_profile', label=r'$I(r)$',
               cmap='cividis', limits=[1e13, 1e30], log=True)
    vr =   function_average(simulation.radial_velocity(file_name), simulation.dim,
                            'Omega', dOmega)
    rho =  function_average(simulation.rho(file_name), simulation.dim,
                            'Omega', dOmega)
    t =    function_average(simulation.temperature(file_name), simulation.dim,
                            'Omega', dOmega)
    p =    function_average(simulation.gas_pressure(file_name), simulation.dim,
                            'Omega', dOmega)
    s =    function_average(simulation.entropy(file_name), simulation.dim,
                            'Omega', dOmega)
    ye =   function_average(simulation.Ye(file_name), simulation.dim,
                            'Omega', dOmega)
    cs =   function_average(simulation.soundspeed(file_name), simulation.dim,
                            'Omega', dOmega)
    eint = function_average(simulation.internal_energy(file_name), simulation.dim,
                            'Omega', dOmega)
    gpot = function_average(simulation.gravitational_potential(file_name), simulation.dim,
                            'Omega', dOmega)
    fconv = simulation.convective_flux(file_name)

    profiles_dictionary['inertia_profile'].append(I_prof)
    profiles_dictionary['radial_velocity_profile'].append(vr)
    profiles_dictionary['rho_profile'].append(rho)
    profiles_dictionary['temperature_profile'].append(t)
    profiles_dictionary['gas_pressure_profile'].append(p)
    profiles_dictionary['entropy_profile'].append(s)
    profiles_dictionary['Ye_profile'].append(ye)
    profiles_dictionary['soundspeed_profile'].append(cs)
    profiles_dictionary['internal_energy_profile'].append(eint)
    profiles_dictionary['gravitational_potential_profile'].append(gpot)    
    profiles_dictionary['convective_flux_profile'].append(fconv)

def __compute_2D_global_local_quantities(profiles_dictionary: dict,
                                        simulation: Simulation,
                                        file_name: str,
                                        dOmega: aerray,
                                        jx: aerray,
                                        jy:aerray,
                                        jz: aerray,
                                        jtot: aerray) -> None:
    
    def rms_quantity(qt, simulation, file_name, dOmega):
        quant = getattr(simulation, qt)(file_name, diff=True, norm=False)
        nm, lb, cb, lm, lg = quant.name, quant.label, quant.cbar, quant.limits, quant.log
        nm = 'profile_' + nm
        lb = merge_strings(r'$\sqrt{$', lb, r'$^2}$')
        qt_rms = np.sqrt(function_average(quant ** 2, simulation.dim, 'Omega',
                                          dOmega))
        qt_rms.set(name=nm, label=lb, cbar=cb, limits=lm, log=lg)
        return qt_rms

    vt = function_average(simulation.theta_velocity(file_name), simulation.dim,
                            'Omega', dOmega)
    vp = function_average(simulation.phi_velocity(file_name), simulation.dim,
                            'Omega', dOmega)
    vr_rms = rms_quantity('radial_velocity', simulation, file_name, dOmega)
    vt_rms = rms_quantity('theta_velocity', simulation, file_name, dOmega)
    vp_rms = rms_quantity('phi_velocity', simulation, file_name, dOmega)
    jx.set(name='angular_momentum_x_profile', label=r'$J_x(r)$',
           cmap='spring', limits=[-1e45, 1e45], log=True)
    jy.set(name='angular_momentum_y_profile', label=r'$J_y(r)$',
            cmap='summer', limits=[-1e45, 1e45], log=True)
    jz.set(name='angular_momentum_z_profile', label=r'$J_z(r)$',
            cmap='autumn', limits=[1e45, 1e49], log=True)
    jtot.set(name='angular_momentum_tot_profile', label=r'$J_{\rm tot}(r)$',
             cmap='spring', limits=[-1e45, 1e49], log=True)
    profiles_dictionary['theta_velocity_profile'].append(vt)
    profiles_dictionary['phi_velocity_profile'].append(vp)
    profiles_dictionary['radial_velocity_rms_profile'].append(vr_rms)
    profiles_dictionary['theta_velocity_rms_profile'].append(vt_rms)
    profiles_dictionary['phi_velocity_rms_profile'].append(vp_rms)
    profiles_dictionary['angular_momentum_x_profile'].append(jx)
    profiles_dictionary['angular_momentum_y_profile'].append(jy)
    profiles_dictionary['angular_momentum_z_profile'].append(jz)
    profiles_dictionary['angular_momentum_tot_profile'].append(jtot)
    kappa = simulation.epicyclic_frequency(file_name)
    bv_freq = simulation.BV_frequency(file_name)
    so_crit = function_average(kappa + bv_freq, simulation.dim,
                                'Omega', dOmega)
    kappa   = function_average(kappa, simulation.dim,
                                'Omega', dOmega)
    bv_freq = function_average(bv_freq, simulation.dim,
                                'Omega', dOmega)
    omega   = function_average(simulation.omega(file_name), simulation.dim,
                                'Omega', dOmega)
    rossby  = function_average(simulation.Rossby_number(file_name), simulation.dim,
                                'Omega', dOmega)
    profiles_dictionary['Solberg_Høiland_profile'].append(so_crit)
    profiles_dictionary['omega_profile'].append(omega)
    profiles_dictionary['Rossby_number'].append(rossby)
    profiles_dictionary['BV_frequency_profile'].append(bv_freq)
    profiles_dictionary['epicyclic_frequency_profile'].append(kappa)

    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        btot = np.sqrt(br ** 2 + bt ** 2 + bp ** 2)
        btot.set(name='total_magnetic_field_profile',
                 label=r'$|\mathbf{B}|$',
                 limits=[1e8, 1e15],
                 log=True,
                 cmap='inferno')
        btot = function_average(btot, simulation.dim, 'Omega', dOmega)
        br = function_average(br, simulation.dim, 'Omega', dOmega)
        bt = function_average(bt, simulation.dim, 'Omega', dOmega)
        bp = function_average(bp, simulation.dim, 'Omega', dOmega)
        profiles_dictionary['magnetic_field_radial_profile'].append(br)
        profiles_dictionary['magnetic_field_theta_profile'].append(bt)
        profiles_dictionary['magnetic_field_phi_profile'].append(bp)
        profiles_dictionary['magnetic_field_total_profile'].append(btot)

def __compute_3D_global_local_quantities(profiles_dictionary: dict,
                                        simulation: Simulation,
                                        file_name: str,
                                        dOmega: aerray,
                                        n: np.array,
                                        X_perp: np.array,
                                        eR: np.array,
                                        grid_unit: u
                                        ) -> None:
    ephi = np.cross(n[:, None, None, :], eR, axis=0)
    
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    omg_og = simulation.omega(file_name)

    r = simulation.cell.radius(simulation.ghost)

    mask_z = (np.ones_like(vr.value) * n[2] > 0.995)

    vx, vy, vz = simulation.grid.spherical_to_cartesian(vr, vt, vp)
    v_cart = np.array([vx.value, vy.value, vz.value])

    vrot = np.where(mask_z, vp, np.nansum(v_cart * ephi, axis=0) * vr.unit)
    vres = np.where(mask_z, np.sqrt(vr ** 2 + vt ** 2),
                    np.sqrt(np.nansum((v_cart - (vrot.value[None, ...] * ephi)) ** 2,
                              axis=0)) * vr.unit)
    omg = np.where(mask_z, omg_og, vrot / (X_perp * grid_unit))

    domega_dr = IDL_derivative(r, omg)

    kappa2 = 2.0 * omg * (2.0 * omg + r[None, None, :] * domega_dr)
    so = kappa2 + bv_freq
    


