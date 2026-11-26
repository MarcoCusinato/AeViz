from __future__ import annotations

import numpy as np
import scipy.special as sp

from AeViz.simulation import Simulation
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
from AeViz.units import aerray
from AeViz.utils.files.string_utils import apply_symbol, merge_strings


def declare_spherical_harmonics_dicionary(dim):
    sph_harmonics_dictionary = {"decomposition": []}
    if dim > 2:
        sph_harmonics_dictionary["power_harm"] = []
    return sph_harmonics_dictionary


def create_metadata_spherical_harmonics(quantity):
    name = quantity.name + "spherical_harmonic_decomposition"
    namesq = quantity.name + "spherical_harmonic_decomposition_summed"
    label = quantity.label
    cmap = "YlGn"
    if "_" in label:
        labelsq = merge_strings(apply_symbol(label), r"$^{\ell}$")
        label = merge_strings(apply_symbol(label), r"$^{\ell m}$")
    else:
        labelsq = merge_strings(apply_symbol(label), r"$_{\ell m}$")
        label = merge_strings(apply_symbol(label), r"$_{\ell m}$")
    meta = {
        "decomposition": {
            "name": name,
            "label": label,
            "cmap": cmap,
            "limits": [-1, 1],
            "log": False,
        },
        "power_harm": {
            "name": namesq,
            "label": labelsq,
            "cmap": cmap,
            "limits": [-1, 1],
            "log": False,
        },
    }
    return meta


def compute_spherical_harmonics_decomposition(
    sph_harmonics_dictionary: dict,
    SpH: SphericalHarmonics,
    quantity: aerray,
    simulation: Simulation,
    dOmega: aerray,
    lmax: int,
) -> None:
    if simulation.dim == 2:
        __compute_2D_spherical_harmonics(
            sph_harmonics_dictionary,
            simulation,
            quantity,
            dOmega,
            SpH,
            lmax,
        )
    else:
        __compute_3D_spherical_harmonics(
            sph_harmonics_dictionary,
            simulation,
            quantity,
            dOmega,
            SpH,
            lmax,
        )


def __compute_2D_spherical_harmonics(
    sph_harmonics_dictionary: dict,
    simulation: Simulation,
    quantity: aerray,
    dOmega: aerray,
    SpH: SphericalHarmonics,
    lmax: int,
):
    """
    Computes the spherical harmonics decomposition of the density up to
    a value of l equal to lmax summed for all ms.

    Parameters
    ----------
    simulation : Simulation
        smiulation of which to compute the spherical harmonics
        decomposition
    file_name : str
        name of the timestep containing the local output.
    quantity: aerray
        quantity to decompose.
    dOmega : aerray
        solid angle element
    SpH : SphericalHarmonics
        class containing the methods to compute the spherical harmonics
    lmax : int
        maximum l to consider

    Returns
    -------
    aerray
        contains the decomposition in spherical harmonics of the
        timestep
    """
    dec_arr = np.zeros((lmax + 1, quantity.shape[-1]))
    for l in range(lmax + 1):
        m = 0
        Ylm = SpH.Ylm_norm(
            m, l, simulation.grid.theta, simulation.cell.phi(simulation.ghost)
        )
        if Ylm.ndim < quantity.ndim:
            Ylm = Ylm[..., None]
        harm = np.sum(quantity * Ylm * dOmega, axis=(0))
        dec_arr[l, :] = harm
    sph_harmonics_dictionary["decomposition"].append(dec_arr[..., None] * quantity.unit)


def __compute_3D_spherical_harmonics(
    sph_harmonics_dictionary: dict,
    simulation: Simulation,
    quantity: aerray,
    dOmega: aerray,
    SpH: SphericalHarmonics,
    lmax: int,
):
    """
    Computes the spherical harmonics decomposition of the density up to
    a value of l equal to lmax summed for all ms.

    Parameters
    ----------
    simulation : Simulation
        smiulation of which to compute the spherical harmonics
        decomposition
    file_name : str
        name of the timestep containing the local output.
    quantity: aerray
        quantity to decompose.
    dOmega : aerray
        solid angle element
    SpH : SphericalHarmonics
        class containing the methods to compute the spherical harmonics
    lmax : int
        maximum l to consider

    Returns
    -------
    aerray
        contains the decomposition in spherical harmonics of the
        timestep
    """
    dec_arr = np.zeros((int(sp.factorial(5)) + 1, quantity.shape[-1]))
    power_arr = np.zeros((lmax + 1, quantity.shape[-1]))
    harm_index = 0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            Ylm = SpH.Ylm_norm(m, l, simulation.grid.theta, simulation.grid.phi)
            if Ylm.ndim < quantity.ndim:
                Ylm = Ylm[..., None]
            harm = np.sum(quantity * Ylm * dOmega, axis=(0, 1))
            if harm_index < dec_arr.shape[0]:
                dec_arr[harm_index, :] = harm
            power_arr[l, :] += harm**2
            harm_index += 1
    sph_harmonics_dictionary["power_harm"].append(
        np.sqrt(power_arr)[..., None] * quantity.unit
    )
    sph_harmonics_dictionary["decomposition"].append(dec_arr[..., None] * quantity.unit)
