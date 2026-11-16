from __future__ import annotations  # noqa: N999

import numpy as np

from AeViz.simulation.simulation import Simulation
from AeViz.units import aerray, u


def declare_GW_dictionary(dim: int) -> dict:
    """
    Declare the neutrino post-processing dictionary.

    Parameters
    ----------
    dim : int
        Simulation dimension (1, 2 or 3).

    Returns
    -------
    dict
        contains the integrated strain in each region and the radial profile
    """
    gws_dictionary = {}
    if dim == 2:
        gws_dictionary["ne220"] = []
        gws_dictionary["full_hpleq"] = []
        gws_dictionary["ne220_corr"] = []
        gws_dictionary["core_hpleq"] = []
        gws_dictionary["pns_hpleq"] = []
        gws_dictionary["innercore_hpleq"] = []
        gws_dictionary["shock_hpleq"] = []
        gws_dictionary["pns_env_hpleq"] = []
        gws_dictionary["sonic_env_hpleq"] = []
        gws_dictionary["out_pns_hpleq"] = []
        gws_dictionary["shocked_hpleq"] = []

    return gws_dictionary


def GWs_metadata() -> dict:
    """
    Metadata for the gravitational waves postprocessing quantities.

    Returns
    -------
    dict
        Metadata dictionary.
    """
    meta = {}
    meta["ne220"] = {
        "name": "AE220",
        "label": r"$h_{\rm +, eq}(r, t)$",
        "limits": [-5, 5],
        "cmap": "seismic",
        "log": False,
    }
    meta["full_hpleq"] = {
        "name": "full_strain",
        "label": r"$h_{\rm +, eq}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["ne220_corr"] = {
        "name": "correction",
        "label": r"$h_{\rm corr}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["core_hpleq"] = {
        "name": "core_strain",
        "label": r"$h_{\rm +, eq, core}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["pns_hpleq"] = {
        "name": "pns_strain",
        "label": r"$h_{\rm +, eq, PNS}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["innercore_hpleq"] = {
        "name": "innercore_strain",
        "label": r"$h_{\rm +, eq, inn}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["shock_hpleq"] = {
        "name": "shock_strain",
        "label": r"$h_{\rm +, eq, shock}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["pns_env_hpleq"] = {
        "name": "pns_envelope_strain",
        "label": r"$h_{\rm +, eq, PNS env}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["sonic_env_hpleq"] = {
        "name": "sonic_envelope_strain",
        "label": r"$h_{\rm +, eq, se}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["out_pns_hpleq"] = {
        "name": "inn_pns_strain",
        "label": r"$h_{\rm +, eq, PNS-inn}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }
    meta["shocked_hpleq"] = {
        "name": "shocked_strain",
        "label": r"$h_{\rm +, eq, shock env}$",
        "limits": [-100, 100],
        "cmap": "seismic",
        "log": False,
    }

    return meta


def compute_gws_postprocessing(
    gws_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    dOmega: aerray,
    dV: aerray,
    cor_rad: aerray,
    pns_rad: aerray,
    inn_rad: aerray,
    sho_rad: aerray,
):
    if simulation.dim == 2:
        __compute_2D_gws_postprocessing(
            gws_dictionary,
            simulation,
            file_name,
            dV,
            dOmega,
            cor_rad,
            pns_rad,
            inn_rad,
            sho_rad,
        )


def __compute_2D_gws_postprocessing(
    gws_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    dV: aerray,
    dOmega: aerray,
    cor_rad: aerray,
    pns_rad: aerray,
    inn_rad: aerray,
    sho_rad: aerray,
):
    radius = simulation.grid.radius
    ctheta = np.cos(simulation.grid.theta)[:, None]
    rho = simulation.rho(file_name)
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    ne220 = (
        dV
        * radius
        * rho
        * (vr * (3 * ctheta**2 - 1) - 3 * vt * ctheta * np.sqrt(1 - ctheta**2))
    )
    correction = __compute_2D_Zha_correction(dOmega, ctheta, radius, rho, vr)
    gws_dictionary["ne220"].append(np.sum(ne220, axis=0)[..., None])
    gws_dictionary["full_hpleq"].append(np.sum(ne220))
    gws_dictionary["ne220_corr"].append(np.sum(correction, axis=0)[..., None])
    gws_dictionary["core_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, None, cor_rad, 0 * u.cm, 0 * u.cm
        )
    )
    gws_dictionary["pns_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, None, pns_rad, 0 * u.cm, 0 * u.cm
        )
    )
    gws_dictionary["innercore_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, None, inn_rad, 0 * u.cm, 20 * u.km
        )
    )
    gws_dictionary["shock_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, None, sho_rad, 0 * u.cm, 0 * u.cm
        )
    )
    gws_dictionary["pns_env_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, cor_rad, pns_rad, 0 * u.cm, 0 * u.cm
        )
    )
    gws_dictionary["sonic_env_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, cor_rad, inn_rad, 0 * u.cm, 20 * u.km
        )
    )
    gws_dictionary["out_pns_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, pns_rad, inn_rad, 0 * u.cm, 20 * u.km
        )
    )
    gws_dictionary["shocked_hpleq"].append(
        __compute_2D_partial_strain(
            radius, ne220, correction, pns_rad, sho_rad, 0 * u.cm, 0 * u.cm
        )
    )


def __compute_2D_Zha_correction(
    dOmega: aerray,
    ctheta: aerray,
    radius: aerray,
    rho: aerray,
    vr: aerray,
) -> aerray:
    """
    Computes the Zha correction into the strains computed on spherica shells
    """
    P2 = 0.5 * (3 * ctheta**2 - 1) * dOmega[:, None]
    prod = radius[None, :] ** 4 * rho * vr * P2
    return prod


def __compute_2D_Zha_surface_correction(
    r1_ind: int,
    r2_ind: int,
    Zha_corr: aerray,
) -> aerray:
    if r1_ind is None:
        r1 = 0 * Zha_corr.unit
    else:
        r1 = Zha_corr[np.arange(Zha_corr.shape[0]), r1_ind]
        r1 = r1.sum()
    if r2_ind is None:
        r2 = 0 * Zha_corr.unit
    else:
        r2 = Zha_corr[np.arange(Zha_corr.shape[0]), r2_ind]
        r2 = r2.sum()
    return r2 - r1


def __compute_2D_partial_strain(
    radius: aerray,
    ne220: aerray,
    correction: aerray,
    rad1: aerray,
    rad2: aerray,
    r1_add: aerray,
    r2_add: aerray,
):
    if rad1 is None:
        mask1 = None
        rad1_ind = None
    else:
        rad1 += r1_add
        mask1 = radius <= rad1[..., None]
        rad1_ind = np.argmax(radius[None, :] >= rad1[..., None], axis=-1)

    rad2 += r2_add
    mask2 = radius <= rad2[..., None]
    rad2_ind = np.argmax(radius[None, :] >= rad2[..., None], axis=-1)

    if mask1 is not None:
        mask = mask2 & np.logical_not(mask1)
    else:
        mask = mask2

    corr = __compute_2D_Zha_surface_correction(rad1_ind, rad2_ind, correction)

    return np.sum(ne220[mask]) + corr
