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
                                    index_phi: int | None) -> None:
    __reduced_neutrino_output(neutrino_dictionary,
                              simulation,
                              file_name,
                              index500,
                              index_theta,
                              index_phi)
    

def __full_neutrino_output(neutrino_dictionary: dict,
                           simulation: Simulation,
                           file_name: str,
                           index500: int,
                           indexTheta: int | None,
                           indices_theta_avg: tuple[int] | None,
                           dOmega: aerray,
                           indexPhi: int | None,
                           Wh: tuple[tuple[aerray]] | None) -> None:
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
                     simulation.radius[index500] ** 2)
    for inu, nu in enumerate(['nue', 'nua', 'nux']):
        num_nu = np.nansum(Lnu_lab_frame[inu][..., index500, :] / ene, axis=-1)
        lnu = np.nansum(Lnu_lab_frame[inu][..., index500, :], axis=-1)
        l2nu = np.nansum(Lnu_lab_frame[inu][..., index500, :] * ene, axis=-1)
        if indexTheta is None:
            neutrino_dictionary[nu]['Llab_avg'].append((4 * np.pi * lnu).to(u.erg /
                                                                        u.s))
            neutrino_dictionary[nu]['Elab_avg'].append((lnu / num_nu).to(u.MeV))
            neutrino_dictionary[nu]['E2lab_avg'].append((l2nu / num_nu).to(
                u.MeV ** 2))
        else:
            angle = 0
            for it in indexTheta:
                if indexPhi:
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
                num_nu[inu] *= dOmega
                lnu[inu] *= dOmega
                l2nu[inu] *= dOmega
                j1, j2, j3 = indices_theta_avg
                neutrino_dictionary[nu][f'Llab_avg'].append(
                    np.nansum(lnu).to(u.erg / u.s))
                neutrino_dictionary[nu][f'Elab_avg'].append(
                    np.nansum(lnu / num_nu).to(u.MeV))
                neutrino_dictionary[nu][f'E2lab_avg'].append(
                    np.nansum(l2nu / num_nu).to(
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
        if simulation.dim == 2:
            gws = __compute_neutrino_contribution_gws_2D(lnu, Wh)
            neutrino_dictionary[nu]['hdot'].append(gws)
        elif simulation.dim == 3
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
        if indexPhi:
            lum = lum[indexPhi, :, :]
        if indexTheta is None:
            neutrino_dictionary[fl]['Leq'] = equivalent_factor * lum[index500]
        else:
            angle = 0
            for it in indexTheta:
                neutrino_dictionary[fl][f'Leq_{angle}'] = \
                    equivalent_factor * lum[it, index500]
                angle += 10

def __compute_neutrino_contribution_gws_2D(neuflux: aerray,
                                           wh: aerray):
    wh = wh[0][0]
    gw = (2 * c.G / (c.c ** 4) * np.nansum(wh * neuflux)).to(u.cm / u.s)
    return gw

def __compute_neutrino_contribution_gws_3D(neuflux: aerray,
                                           wh: tuple[tuple[aerray]],
                                           dOmega: aerray):
    weq, wpol, wtd = wh

    gwpleq = (2 * c.G / (c.c ** 4) * np.nansum(weq[0] * neuflux)).to(u.cm / u.s)
    gwcreq = (2 * c.G / (c.c ** 4) * np.nansum(weq[1] * neuflux)).to(u.cm / u.s)

    gwplpol = (2 * c.G / (c.c ** 4) * np.nansum(wpol[0] * neuflux)).to(u.cm / u.s)
    gwcrpol = (2 * c.G / (c.c ** 4) * np.nansum(wpol[1] * neuflux)).to(u.cm / u.s)

    l2pl = np.nansum(dOmega * (wtd[0] * neuflux) ** 2)
    l2cr = np.nansum(dOmega * (wtd[1] * neuflux) ** 2)
    return ((gwpleq, gwcreq), (gwplpol, gwcrpol), (l2pl, l2cr)) 

def compute_weights(theta: aerray, phi: aerray, alpha: aerray, beta: aerray, dim: int) -> tuple[aerray]:
    if dim == 1:
        whpl, whcr = None, None
    elif dim == 2:
        whpl = np.sin(theta) * np. pi * (2 * np.abs(np.cos(theta)) - 1)
        whcr = 0
    else:
        theta = theta[None, :]
        phi = phi[:, None]
        if alpha.size > 1 and beta.size > 1:
            alpha = alpha[None, :]
            beta = beta[:, None]
        N = (((np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * np.sin(theta) * np.cos(beta) -
             np.cos(theta) * np.sin(beta)) ** 2 + np.sin(theta) ** 2 *
             (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)) ** 2)
        Dpl = ((1 + (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * np.sin(theta) * np.sin(beta) +
                np.cos(theta) * np.cos(beta)) * 
                (((np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * 
                  np.sin(alpha) * np.cos(beta) - np.cos(theta) * np.sin(beta)) ** 2 - 
                np.sin(theta) ** 2 * (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)) ** 2))
        Dcr = ((1 + (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha)) * np.sin(theta) * np.sin(beta) +
                np.cos(theta) * np.cos(beta)) * 2 * ((np.cos(phi) * np.cos(alpha) + np.sin(theta) * np.sin(alpha)) *
                                                     np.sin(theta) * np.cos(beta)) * np.sin(theta) *
                                                     (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)))
        whpl = Dpl / N
        whcr = Dcr / N
    return (whpl, whcr)