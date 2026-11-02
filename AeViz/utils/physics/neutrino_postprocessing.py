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
                                    dV: aerray,) -> None:
    index500 = np.argmax(simulation.grid.radius >= (500 * u.km))
    index_theta = None
    index_phi = None
    if simulation.dim > 1:
        index_theta = [np.argmax(simulation.grid.theta >= (dd * u.deg))
                       for dd in (np.arange(19) * 10)]
        index_theta[-1] = -1 if index_theta[-1] == 0 else index_theta[-1]
    if simulation.dim > 2:
        index_phi = 0

def __full_neutrino_output(neutrino_dictionary: dict,
                           simulation: Simulation,
                           file_name: str,
                           index500: int,
                           indexTheta: int | None,
                           indices_theta_avg: tuple[int] | None,
                           dOmega: aerray,
                           indexPhi: int | None) -> None:
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
    Lnu_lab_frame = (simulation.labframe_neutrino_luminosity(file_name) *
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