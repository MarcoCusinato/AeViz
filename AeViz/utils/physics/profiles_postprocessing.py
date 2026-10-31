from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.units.constants import constants as c
from AeViz.utils.math_utils import IDL_derivative
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
        'gravitational_potential_profile': []
    }
    if dim > 1:
        profiles_dictionary['theta_velocity_profile'] = []
        profiles_dictionary['phi_velocity_profile'] = []
        profiles_dictionary['radial_velocity_rms_profile'] = []
        profiles_dictionary['theta_velocity_rms_profile'] = []
        profiles_dictionary['phi_velocity_rms_profile'] = []
        profiles_dictionary['convective_flux_profile'] = []
        profiles_dictionary['angular_momentum_x_profile'] = []
        profiles_dictionary['angular_momentum_y_profile'] = []
        profiles_dictionary['angular_momentum_z_profile'] = []
        profiles_dictionary['angular_momentum_tot_profile'] = []
        profiles_dictionary['Solberg_Hoiland_profile'] = []
        profiles_dictionary['omega_profile'] = []
        profiles_dictionary['Rossby_number'] = []
        profiles_dictionary['BV_frequency_profile'] = []
        profiles_dictionary['epicyclic_frequency_profile'] = []
        profiles_dictionary['vorticity_magnitude_profile'] = []
        profiles_dictionary['vorticity_x_profile'] = []
        profiles_dictionary['vorticity_y_profile'] = []
        profiles_dictionary['vorticity_z_profile'] = []
        profiles_dictionary['vortensity_magnitude_profile'] = []
        profiles_dictionary['vortensity_x_profile'] = []
        profiles_dictionary['vortensity_y_profile'] = []
        profiles_dictionary['vortensity_z_profile'] = []
        profiles_dictionary['enstrophy_profile'] = []
        if magdim > 1:
            profiles_dictionary['magnetic_field_radial_profile'] = []
            profiles_dictionary['magnetic_field_theta_profile'] = []
            profiles_dictionary['magnetic_field_phi_profile'] = []
            profiles_dictionary['magnetic_field_total_profile'] = []
    if dim > 2:
        profiles_dictionary['rotational_velocity_profile'] = []
        profiles_dictionary['residual_velocity_profile'] = []
            
    return profiles_dictionary

def create_profiles_metadata(profile_dictionary:dict) -> dict:
    meta = {
        key: dict(
            name=value[0].name,
            label=value[0].label,
            limits=value[0].limits,
            cmap=value[0].cmap,
            log=value[0].log
         ) for key, value in profile_dictionary.items()
    }
    return meta

def compute_profiles_postprocessing(profiles_dictionary: dict,
                                    simulation: Simulation,
                                    file_name: str,
                                    jx: aerray | None,
                                    jy: aerray | None,
                                    jz: aerray | None,
                                    jtot: aerray | None,
                                    I: aerray,
                                    n: np.array,
                                    X_perp: np.array,
                                    eR: np.array,
                                    grid_unit: u) -> None:
    __compute_1D_profiles(profiles_dictionary,
                          simulation,
                          file_name,
                          I)
    if simulation.dim > 1:
        __compute_2D_profiles(profiles_dictionary,
                              simulation,
                              file_name,
                              jx,
                              jy,
                              jz,
                              jtot)
    if simulation.dim > 2:
        __compute_3D_profiles(profiles_dictionary,
                              simulation,
                              file_name,
                              n,
                              X_perp,
                              eR,
                              grid_unit)

def __compute_1D_profiles(profiles_dictionary: dict,
                                        simulation: Simulation,
                                        file_name: str,
                                        I: aerray) -> None:
    
    I_prof = np.nansum(I, axis=tuple(range(
        simulation.dim - 1))).to(u.g * u.cm ** 2)
    I_prof.set(name='inertia_moment_profile', label=r'$I(r)$',
               cmap='cividis', limits=[1e13, 1e30], log=True)
    vr =   simulation.grid.omega_average(simulation.radial_velocity(file_name))
    rho =  simulation.grid.omega_average(simulation.rho(file_name))
    t =    simulation.grid.omega_average(simulation.temperature(file_name))
    p =    simulation.grid.omega_average(simulation.gas_pressure(file_name))
    s =    simulation.grid.omega_average(simulation.entropy(file_name))
    ye =   simulation.grid.omega_average(simulation.Ye(file_name))
    cs =   simulation.grid.omega_average(simulation.soundspeed(file_name))
    eint = simulation.grid.omega_average(simulation.internal_energy(file_name))
    gpot = simulation.grid.omega_average(simulation.gravitational_potential(file_name))

    profiles_dictionary['inertia_profile'].append(I_prof[..., None])
    profiles_dictionary['radial_velocity_profile'].append(vr[..., None])
    profiles_dictionary['rho_profile'].append(rho[..., None])
    profiles_dictionary['temperature_profile'].append(t[..., None])
    profiles_dictionary['gas_pressure_profile'].append(p[..., None])
    profiles_dictionary['entropy_profile'].append(s[..., None])
    profiles_dictionary['Ye_profile'].append(ye[..., None])
    profiles_dictionary['soundspeed_profile'].append(cs[..., None])
    profiles_dictionary['internal_energy_profile'].append(eint[..., None])
    profiles_dictionary['gravitational_potential_profile'].append(gpot[..., None])

def __compute_2D_profiles(profiles_dictionary: dict,
                                        simulation: Simulation,
                                        file_name: str,
                                        jx: aerray,
                                        jy:aerray,
                                        jz: aerray,
                                        jtot: aerray) -> None:
    
    def rms_quantity(qt, simulation, file_name):
        quant = getattr(simulation, qt)(file_name, diff=True, norm=False)
        nm, lb, cb, lm, lg = quant.name, quant.label, quant.cbar, quant.limits, quant.log
        nm = 'profile_' + nm
        lb = merge_strings(r'$\sqrt{$', lb, r'$^2}$')
        qt_rms = np.sqrt(simulation.grid.omega_average(quant ** 2))
        qt_rms.set(name=nm, label=lb, cbar=cb, limits=lm, log=lg)
        return qt_rms

    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    fconv = simulation.convective_flux(file_name)
    rho = simulation.rho(file_name)
    vr_rms = rms_quantity('radial_velocity', simulation, file_name)
    vt_rms = rms_quantity('theta_velocity', simulation, file_name)
    vp_rms = rms_quantity('phi_velocity', simulation, file_name)
    jx.set(name='angular_momentum_x_profile', label=r'$J_x(r)$',
           cmap='spring', limits=[-1e45, 1e45], log=True)
    jy.set(name='angular_momentum_y_profile', label=r'$J_y(r)$',
            cmap='summer', limits=[-1e45, 1e45], log=True)
    jz.set(name='angular_momentum_z_profile', label=r'$J_z(r)$',
            cmap='autumn', limits=[1e45, 1e49], log=True)
    jtot.set(name='angular_momentum_tot_profile', label=r'$J_{\rm tot}(r)$',
             cmap='spring', limits=[-1e45, 1e49], log=True)
    #compute the vorticity, entrophy and vortensity
    curl_v = simulation.grid.curl([vr, vt, vp])
    enstrophy = 0.5 * (curl_v[0] ** 2 + curl_v[1] ** 2 + curl_v[2] ** 2)
    enstrophy.set(name='enstrophy', label=r'$\varepsilon$',
                  cmap='inferno', limits=[1e0, 1e8], log=True)
    vort_mag = np.sqrt(curl_v[0] ** 2 + curl_v[1] ** 2 + curl_v[2] ** 2)
    vort_mag.set(name='vorticity_magnitude', label=r'$|\omega|$',
                  cmap='magma', limits=[1e0, 1e4], log=True)
    vort_x = simulation.grid.project(curl_v, 'x')
    vort_y = simulation.grid.project(curl_v, 'y')
    vort_z = simulation.grid.project(curl_v, 'z')
    vort_x.set(name='vorticity_x', label=r'$\omega_x$',
              cmap='bwr', limits=[-1e3, 1e3], log=True)
    vort_y.set(name='vorticity_y', label=r'$\omega_y$',
              cmap='bwr', limits=[-1e3, 1e3], log=True)
    vort_z.set(name='vorticity_z', label=r'$\omega_z$',
              cmap='bwr', limits=[-1e3, 1e3], log=True)
    vortens_mag = vort_mag / rho
    vortens_mag.set(name='vortensity_magnitude', label=r'$|\xi|$',
                    cmap='gray', limits=[1e17, 1e13], log=True)
    vortensity_x = vort_x / rho
    vortensity_y = vort_y / rho
    vortensity_z = vort_z / rho
    vortensity_x.set(name='vortensity_x', label=r'$\xi_x$',
                     cmap='seismic', limits=[-1e14, 1e14], log=True)
    vortensity_y.set(name='vortensity_y', label=r'$\xi_y$',
                     cmap='seismic', limits=[-1e14, 1e14], log=True)
    vortensity_z.set(name='vortensity_z', label=r'$\xi_z$',
                     cmap='seismic', limits=[-1e14, 1e14], log=True)

    profiles_dictionary['theta_velocity_profile'].append(simulation.grid.omega_average(vt)[..., None])
    profiles_dictionary['phi_velocity_profile'].append(simulation.grid.omega_average(vp)[..., None])
    profiles_dictionary['radial_velocity_rms_profile'].append(vr_rms[..., None])
    profiles_dictionary['theta_velocity_rms_profile'].append(vt_rms[..., None])
    profiles_dictionary['phi_velocity_rms_profile'].append(vp_rms[..., None])
    profiles_dictionary['convective_flux_profile'].append(fconv[..., None])
    profiles_dictionary['angular_momentum_x_profile'].append(jx[..., None])
    profiles_dictionary['angular_momentum_y_profile'].append(jy[..., None])
    profiles_dictionary['angular_momentum_z_profile'].append(jz[..., None])
    profiles_dictionary['angular_momentum_tot_profile'].append(jtot[..., None])
    profiles_dictionary['vorticity_magnitude_profile'].append(simulation.grid.omega_average(vort_mag)[..., None])
    profiles_dictionary['vorticity_x_profile'].append(simulation.grid.omega_average(vort_x)[..., None])
    profiles_dictionary['vorticity_y_profile'].append(simulation.grid.omega_average(vort_y)[..., None])
    profiles_dictionary['vorticity_z_profile'].append(simulation.grid.omega_average(vort_z)[..., None])
    profiles_dictionary['vortensity_magnitude_profile'].append(simulation.grid.omega_average(vortens_mag)[..., None])
    profiles_dictionary['vortensity_x_profile'].append(simulation.grid.omega_average(vortensity_x)[..., None])
    profiles_dictionary['vortensity_y_profile'].append(simulation.grid.omega_average(vortensity_y)[..., None])
    profiles_dictionary['vortensity_z_profile'].append(simulation.grid.omega_average(vortensity_z)[..., None])
    profiles_dictionary['enstrophy_profile'].append(simulation.grid.omega_average(enstrophy)[..., None])

    if simulation.dim == 2:
        kappa = simulation.epicyclic_frequency(file_name)
        bv_freq = simulation.BV_frequency(file_name)
        so_crit = kappa + bv_freq
        so_crit.set(name='soldberg hoiland criterion',
            label=r'$\omega^2_{\rm SH}$',
            limits=bv_freq.limits,
            cmap='seismic',
            log=True)
        kappa   = simulation.grid.omega_average(kappa)
        bv_freq = simulation.grid.omega_average(bv_freq)
        omega   = simulation.grid.omega_average(simulation.omega(file_name))
        rossby  = simulation.grid.omega_average(simulation.Rossby_number(file_name))
        profiles_dictionary['omega_profile'].append(omega[..., None])
        profiles_dictionary['Solberg_Hoiland_profile'].append(simulation.grid.omega_average(so_crit)[..., None])
        profiles_dictionary['Rossby_number'].append(rossby[..., None])
        profiles_dictionary['BV_frequency_profile'].append(bv_freq[..., None])
        profiles_dictionary['epicyclic_frequency_profile'].append(kappa[..., None])

    if simulation.evolved_qts['magdim'] > 0:
        br, bt, bp = simulation.magnetic_fields(file_name, comp='all')
        btot = np.sqrt(br ** 2 + bt ** 2 + bp ** 2)
        btot.set(name='total_magnetic_field_profile',
                 label=r'$|\mathbf{B}|$',
                 limits=[1e8, 1e15],
                 log=True,
                 cmap='inferno')
        btot = simulation.grid.omega_average(btot)
        br =   simulation.grid.omega_average(br)
        bt =   simulation.grid.omega_average(bt)
        bp =   simulation.grid.omega_average(bp)
        profiles_dictionary['magnetic_field_radial_profile'].append(br[..., None])
        profiles_dictionary['magnetic_field_theta_profile'].append(bt[..., None])
        profiles_dictionary['magnetic_field_phi_profile'].append(bp[..., None])
        profiles_dictionary['magnetic_field_total_profile'].append(btot[..., None])

def __compute_3D_profiles(profiles_dictionary: dict,
                          simulation: Simulation,
                          file_name: str,
                          n: np.array,
                          X_perp: np.array,
                          eR: np.array,
                          grid_unit: u
                          ) -> None:
    ephi = np.cross(n[:, None, None, :], eR, axis=0)
    
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    vconv = simulation.convective_velocity(file_name)
    rho = simulation.rho(file_name)
    kappa_og = simulation.epicyclic_frequency(file_name)
    bv_freq = simulation.BV_frequency(file_name)
    so_crit_og = kappa_og + bv_freq
    so_crit.set(name='soldberg hoiland criterion',
                label=r'$\omega^2_{\rm SH}$',
                limits=bv_freq.limits,
                cmap='seismic',
                log=True)
    kappa_og   = simulation.grid.omega_average(kappa_og)
    omega_og   = simulation.omega(file_name)
    rossby_og  = simulation.grid.omega_average(simulation.Rossby_number(file_name))
    profiles_dictionary['omega_profile'].append(simulation.grid.omega_average(omega_og)[..., None])
    profiles_dictionary['Solberg_Hoiland_profile'].append(simulation.grid.omega_average(so_crit_og)[..., None])
    profiles_dictionary['Rossby_number'].append(rossby_og[..., None])
    profiles_dictionary['BV_frequency_profile'].append(simulation.grid.omega_average(bv_freq)[..., None])
    profiles_dictionary['epicyclic_frequency_profile'].append(kappa_og[..., None])

    ## Follow the shell axis
    r = simulation.grid.radius
    H = 1 / np.abs(IDL_derivative(r, rho) / rho * r)
    mask_z = (np.ones_like(vr.value) * n[2] > 0.995)

    vx, vy, vz = simulation.grid.spherical_to_cartesian(vr, vt, vp)
    v_cart = np.array([vx.value, vy.value, vz.value])
    vrot = np.where(mask_z, vp, np.nansum(v_cart * ephi, axis=0) * vr.unit)
    vrot.set(name  =vp.name,
             label =r'$v_{\rm rot}$',
             limits=vp.limits,
             cmap  =vp.cmap,
             log   =vp.log)
    vres = np.where(mask_z, np.sqrt(vr ** 2 + vt ** 2),
                    np.sqrt(np.nansum((v_cart - (vrot.value[None, ...] * ephi)) ** 2,
                              axis=0)) * vr.unit)
    vres.set(name  =vconv.name,
             label =r'$v_{\rm rot}$',
             limits=vconv.limits,
             cmap  =vconv.cmap,
             log   =vconv.log)
    omg = np.where(mask_z, omega_og, vrot / (X_perp * grid_unit))
    omg.set(name  =omega_og.name,
            label =omega_og.label,
            limits=omega_og.limits,
            cmap  =omega_og.cmap,
            log   =omega_og.log)
    ross = vconv / (r * omg * H) * u.radian
    ross.set(name='Rossby', label=r'$Ro$',
            cmap='RdYlBu_r', log=False, limits=[0, 1])
    
    domega_dr = IDL_derivative(r, omg)
    kappa2 = 2.0 * omg * (2.0 * omg + r[None, None, :] * domega_dr)
    kappa2.set(name  =kappa_og.name,
               label =kappa_og.label,
               limits=kappa_og.limits,
               cmap  =kappa_og.cmap,
               log   =kappa_og.log)
    so_crit = kappa2 + bv_freq
    so_crit.set(name='soldberg hoiland criterion',
                label=r'$\omega^2_{\rm SH}$',
                limits=bv_freq.limits,
                cmap='seismic',
                log=True)
    profiles_dictionary['omega_axis_profile'].append(simulation.grid.omega_average(omg)[..., None])
    profiles_dictionary['Solberg_Hoiland_axis_profile'].append(simulation.grid.omega_average(so_crit)[..., None])
    profiles_dictionary['Rossby_axis_profile'].append(simulation.grid.omega_average(ross)[..., None])
    profiles_dictionary['epicyclic_frequency_axis_profile'].append(simulation.grid.omega_average(kappa2)[..., None])
    profiles_dictionary['rotational_velocity_profile'].append(simulation.grid.omega_average(vrot)[..., None])
    profiles_dictionary['residual_velocity_profile'].append(simulation.grid.omega_average(vres)[..., None])


