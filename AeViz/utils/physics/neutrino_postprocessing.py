from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import neutrino_sphere_radii
from AeViz.units.constants import constants as c
from AeViz.cell.cell_methods.spherical_methods import dOmega




def compute_neutrino_postprocessing(neutrino_dictionary: dict,
                                    simulation: Simulation,
                                    file_name: str,
                                    dOmega: aerray,
                                    Wh: tuple[tuple[aerray]] | None,
                                    index500: int,
                                    index_theta: int | None,
                                    indices_theta_avg: tuple[int] | None,
                                    index_phi: int | None,
                                    PNS_mass: aerray | None) -> None:
    __reduced_neutrino_output(neutrino_dictionary,
                              simulation,
                              file_name,
                              index500,
                              index_theta,
                              index_phi)
    __full_neutrino_output(neutrino_dictionary,
                           simulation,
                           file_name,
                           index500,
                           index_theta,
                           indices_theta_avg,
                           dOmega,
                           index_phi,
                           Wh,
                           PNS_mass
                           )

def __full_neutrino_output(neutrino_dictionary: dict,
                           simulation: Simulation,
                           file_name: str,
                           index500: int,
                           indexTheta: int | None,
                           indicesThetaAvg: tuple[int] | None,
                           dOmega: aerray,
                           indexPhi: int | None,
                           Wh: tuple[tuple[aerray]] | None,
                           PNS_mass: aerray | None) -> None:
    ene = simulation.cell.E_nu()
    while ene.ndim < simulation.dim:
        ene = ene[None, ...]
    nuer, nuar, nuxr = neutrino_sphere_radii(simulation,
                                             file_name)
    neutrino_dictionary['nue']['r'].append(nuer)
    neutrino_dictionary['nue']['rmin'].append(np.nanmin(nuer))
    neutrino_dictionary['nue']['rmax'].append(np.nanmax(nuer))
    neutrino_dictionary['nue']['ravg'].append(
        simulation.grid.radii_angular_average(nuer))
    neutrino_dictionary['nua']['r'].append(nuar)
    neutrino_dictionary['nua']['rmin'].append(np.nanmin(nuar))
    neutrino_dictionary['nua']['rmax'].append(np.nanmax(nuar))
    neutrino_dictionary['nua']['ravg'].append(
        simulation.grid.radii_angular_average(nuar))
    neutrino_dictionary['nux']['r'].append(nuxr)
    neutrino_dictionary['nux']['rmin'].append(np.nanmin(nuxr))
    neutrino_dictionary['nux']['rmax'].append(np.nanmax(nuxr))
    neutrino_dictionary['nux']['ravg'].append(
        simulation.grid.radii_angular_average(nuxr))
    Lnu_lab_frame = (simulation.labframe_neutrino_fluxes(file_name) *
                     simulation.grid.radius[index500] ** 2)
    for inu, nu in enumerate(['nue', 'nua', 'nux']):
        lumlab = Lnu_lab_frame[inu][..., index500, :] if simulation.dim == 1 \
            else Lnu_lab_frame[inu][..., index500, :, 0]
        num_nu = np.nansum(lumlab / ene, axis=-1)
        lnu = np.nansum(lumlab, axis=-1)
        l2nu = np.nansum(lumlab * ene, axis=-1)
        if simulation.dim == 1:
            neutrino_dictionary[nu]['Llab_avg'].append((4 * np.pi * lnu).to(u.erg /
                                                                        u.s))
            neutrino_dictionary[nu]['Elab_avg'].append((lnu / num_nu).to(u.MeV))
            neutrino_dictionary[nu]['E2lab_avg'].append((l2nu / num_nu).to(
                u.MeV ** 2))
        else:
            angle = 0
            for it in indexTheta:
                if simulation.dim == 3:
                    neutrino_dictionary[nu][f'Llab_{angle}'].append(
                        (4 * np.pi * lnu[indexPhi, it]).to(u.erg / u.s))
                    neutrino_dictionary[nu][f'Elab_{angle}'].append(
                        (lnu[indexPhi, it] / num_nu[indexPhi, it]).to(u.MeV))
                    neutrino_dictionary[nu][f'E2lab_{angle}'].append(
                        (l2nu[indexPhi, it] / num_nu[indexPhi, it]).to(
                            u.MeV ** 2))
                else:
                    neutrino_dictionary[nu][f'Llab_{angle}'].append(
                        (4 * np.pi * lnu[it]).to(u.erg / u.s))
                    neutrino_dictionary[nu][f'Elab_{angle}'].append(
                        (lnu[it] / num_nu[it]).to(u.MeV))
                    neutrino_dictionary[nu][f'E2lab_{angle}'].append(
                        (l2nu[it] / num_nu[it]).to(
                            u.MeV ** 2))
                angle += 10
            num_nu *= dOmega
            lnu *= dOmega
            l2nu *= dOmega
            j1, j2, j3 = indicesThetaAvg
            neutrino_dictionary[nu][f'Llab_avg'].append(
                np.nansum(lnu).to(u.erg / u.s))
            neutrino_dictionary[nu][f'Elab_avg'].append(
                (np.nansum(lnu) / np.nansum(num_nu)).to(u.MeV))
            neutrino_dictionary[nu][f'E2lab_avg'].append(
                (np.nansum(l2nu) / np.nansum(num_nu)).to(
                    u.MeV ** 2))
            neutrino_dictionary[nu][f'Llab_pol'].append(
                (4 * np.pi * np.nansum(lnu[..., :j1]) / 
                    np.nansum(dOmega[..., :j1])).to(u.erg / u.s))
            neutrino_dictionary[nu][f'Elab_pol'].append(
                (np.nansum(lnu[..., :j1]) / np.nansum(num_nu[..., :j1])).to(u.MeV))
            neutrino_dictionary[nu][f'E2lab_pol'].append(
                (np.nansum(l2nu[..., :j1]) / np.nansum(num_nu[..., :j1])).to(
                    u.MeV ** 2))
            neutrino_dictionary[nu][f'Llab_eq'].append(
                (4 * np.pi * np.nansum(lnu[..., j2:j3]) / 
                    np.nansum(dOmega[..., j2:j3])).to(u.erg / u.s))
            neutrino_dictionary[nu][f'Elab_eq'].append(
                (np.nansum(lnu[..., j2:j3]) / 
                    np.nansum(num_nu[..., j2:j3])).to(u.MeV))
            neutrino_dictionary[nu][f'E2lab_eq'].append(
                (np.nansum(l2nu[..., j2:j3]) / 
                    np.nansum(num_nu[..., j2:j3])).to(
                    u.MeV ** 2))
            __compute_kick_ang_mom(neutrino_dictionary[nu],
                                   simulation,
                                   np.nansum(Lnu_lab_frame[inu], axis=-2),
                                   dOmega,
                                   PNS_mass,
                                   simulation.grid.radius[index500])
            if simulation.dim == 2:
                ## lnu already multiplied by dOmega
                gws = __compute_neutrino_contribution_gws_2D(lnu, Wh)
                neutrino_dictionary[nu]['hdot_eqpl'].append(gws)
            elif simulation.dim == 3:
                gweq, gwpol, gwene = __compute_neutrino_contribution_gws_3D(lnu,
                                                                            Wh,
                                                                            dOmega)
                neutrino_dictionary[nu]['hdot_eqpl'].append(gweq[0])
                neutrino_dictionary[nu]['hdot_eqcr'].append(gweq[1])
                neutrino_dictionary[nu]['hdot_polpl'].append(gwpol[0])
                neutrino_dictionary[nu]['hdot_polcr'].append(gwpol[1])
                neutrino_dictionary[nu]['l2pl'].append(gwene[0])
                neutrino_dictionary[nu]['l2cr'].append(gwene[1])

def __reduced_neutrino_output(neutrino_dictionary: dict,
                              simulation: Simulation,
                              file_name: str,
                              index500: int,
                              indexTheta: int | None,
                              indexPhi: int | None) -> None:
    equivalent_factor = np.pi * 4 * simulation.grid.radius[index500] ** 2
    nulums = simulation.neutrino_luminosity(file_name)
    for fl, lum in zip(['nue', 'nua', 'nux'], nulums):
        if indexPhi is not None:
            lum = lum[indexPhi, :, :]
        if indexTheta is None:
            neutrino_dictionary[fl]['Leq'].append(equivalent_factor * lum[index500])
        else:
            angle = 0
            for it in indexTheta:
                neutrino_dictionary[fl][f'Leq_{angle}'].append(
                    equivalent_factor * lum[it, index500])
                angle += 10

def __compute_neutrino_contribution_gws_2D(neuflux: aerray,
                                           wh: aerray):
    wh = wh[0][0]
    gw = (2 * c.G / (c.c ** 4) * np.nansum(wh * neuflux)).to(u.cm / u.s)
    return gw

def __compute_neutrino_contribution_gws_3D(neuflux: aerray,
                                           wh: tuple[tuple[aerray]],
                                           dOmega_obs: aerray):
    weq, wpol, wtd = wh

    gwpleq = (2 * c.G / (c.c ** 4) * np.nansum(weq[0] * neuflux)).to(u.cm / u.s)
    gwcreq = (2 * c.G / (c.c ** 4) * np.nansum(weq[1] * neuflux)).to(u.cm / u.s)

    gwplpol = (2 * c.G / (c.c ** 4) * np.nansum(wpol[0] * neuflux)).to(u.cm / u.s)
    gwcrpol = (2 * c.G / (c.c ** 4) * np.nansum(wpol[1] * neuflux)).to(u.cm / u.s)

    l2pl = np.nansum(dOmega_obs * (wtd[0] * neuflux[None, None, ...]) ** 2)
    l2cr = np.nansum(dOmega_obs * (wtd[1] * neuflux[None, None, ...]) ** 2)
    return ((gwpleq, gwcreq), (gwplpol, gwcrpol), (l2pl, l2cr)) 

def __compute_kick_ang_mom(neutrino_dictionary: dict,
                           simulation: Simulation,
                           nufluxes: aerray,
                           dOmega:aerray,
                           PNS_mass: aerray,
                           r500: aerray):
    nufluxes = nufluxes * dOmega[..., None]
    ## taking only the radial flux since we are far enogh and the neutrinos are free-streaming
    fx, fy, fz = simulation.grid.spherical_to_cartesian(
        nufluxes[..., 0][..., None],
        np.zeros_like(nufluxes[..., 0][..., None]),
        np.zeros_like(nufluxes[..., 0][..., None])
    )
    accz =  (- np.nansum(fz) / (PNS_mass * c.c)).to(u.cm / u.s ** 2)
    accx =  (- np.nansum(fx) / (PNS_mass * c.c)).to(u.cm / u.s ** 2) \
        if simulation.dim == 3 else 0 * accz.unit
    accy =  (- np.nansum(fy) / (PNS_mass * c.c)).to(u.cm / u.s ** 2) \
        if simulation.dim == 3 else 0 * accz.unit
    
    neutrino_dictionary['accx'].append(accx.to(u.cm / u.s ** 2))
    neutrino_dictionary['accy'].append(accy.to(u.cm / u.s ** 2))
    neutrino_dictionary['accz'].append(accz.to(u.cm / u.s ** 2))

    ## compute the neutrino contribution to the ang mom
    ## as 0, +fphi, -ftheta
    jx, jy, jz = simulation.grid.spherical_to_cartesian(
        np.zeros_like(nufluxes[..., 0][..., None]),
        nufluxes[..., 2][..., None],
        -nufluxes[..., 1][..., None]
    )
    jz = r500 / c.c * np.nansum(jz)
    jx = r500 / c.c * np.nansum(jx) if simulation.dim == 3 else 0 * jz.unit
    jy = r500 / c.c * np.nansum(jy) if simulation.dim == 3 else 0 * jz.unit
    neutrino_dictionary['jx'].append(jx)
    neutrino_dictionary['jy'].append(jy)
    neutrino_dictionary['jz'].append(jz)

def compute_weights(theta: aerray,
                    phi: aerray,
                    alpha: aerray,
                    beta: aerray,
                    dim: int) -> tuple[aerray]:
    if dim == 1:
        whpl, whcr = None, None
    elif dim == 2:
        whpl = 0.5 * (2 * np.abs(np.cos(theta)) - 1)
        whcr = 0
    else:
        theta = theta[None, :]
        phi = phi[:, None]
        if alpha.size > 1 and beta.size > 1:
            alpha = alpha[:, None, None, None]
            beta = beta[None, :, None, None]
            theta = theta[None, None, ...]
            phi = phi[None, None, ...]
        N = (((np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * 
              np.sin(theta) * np.cos(beta) -
             np.cos(theta) * np.sin(beta)) ** 2 + np.sin(theta) ** 2 *
             (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)) ** 2)
        Dpl = ((1 + (np.cos(phi) * np.cos(alpha) + np.sin(phi) *
                     np.sin(alpha)) * np.sin(theta) * np.sin(beta) +
                np.cos(theta) * np.cos(beta)) * 
                (((np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * 
                  np.sin(theta) * np.cos(beta) - np.cos(theta) * np.sin(beta)) ** 2 - 
                np.sin(theta) ** 2 * (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)) ** 2))
        Dcr = ((1 + (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) *
                np.sin(theta) * np.sin(beta) + np.cos(theta) * np.cos(beta)) * 
               2 * ((np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * 
                    np.sin(theta) * np.cos(beta)- np.cos(theta) * np.sin(beta)) * 
                np.sin(theta) * (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)))
        whpl = Dpl / N
        whcr = Dcr / N
    return (whpl, whcr)

def build_obs_grid(nalpha: int):
    assert (nalpha % 2 == 0), "Number of alphas must be even."
    dangle = 2 * np.pi / nalpha * u.radian
    beta = (np.arange(int(nalpha / 2)) + 0.5) * dangle
    beta_sinb = np.sin(beta)[:, None, None, None] * dangle
    dOmega_obs = (np.ones(int(nalpha / 2))[None, :, None, None] * dangle) * beta_sinb
    return dOmega_obs
## POINT B2