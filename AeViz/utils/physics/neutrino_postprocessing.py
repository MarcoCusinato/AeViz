from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import neutrino_sphere_radii
from AeViz.units.constants import constants as c


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
                           file_name: str) -> None:
    ene = simulation.cell.E_nu()
    while ene.ndim < simulation.dim:
        ene = ene[None, ...]
    nuer, nuar, nuxr = neutrino_sphere_radii(simulation,
                                             file_name)
    Lnu_lab_frame = simulation.labframe_neutrino_luminosity(file_name)

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