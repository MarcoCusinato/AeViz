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
    gws_dictionary["matter_strain"] = []
    gws_dictionary["full_hpleq"] = []
    gws_dictionary["corr"] = []
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
    meta["matter_strain"] = {
        "name": "matter_strain",
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
    meta["corr"] = {
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
    dOmega: aerray | None,
    dV: aerray,
    cor_rad: aerray,
    pns_rad: aerray,
    inn_rad: aerray,
    sho_rad: aerray,
    gradient: aerray | None,
    Ylm: aerray | None,
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
    elif simulation.dim == 3:
        __compute_3D_gws_postprocessing(
            gws_dictionary,
            simulation,
            file_name,
            dV,
            gradient,
            Ylm,
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
    gws_dictionary["matter_strain"].append(np.sum(ne220, axis=0)[..., None])
    gws_dictionary["full_hpleq"].append(np.sum(ne220))
    gws_dictionary["corr"].append(np.sum(correction, axis=0)[..., None])
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


def __compute_3D_gws_postprocessing(
    gws_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    dV: aerray,
    gradient: aerray,
    Ylm: aerray,
    cor_rad: aerray,
    pns_rad: aerray,
    inn_rad: aerray,
    sho_rad: aerray,
):
    radius = simulation.grid.radius
    inn_rad += 20 * u.km
    mask_co = radius <= cor_rad[..., None]
    ind_co = np.argmax(radius[None, None, :] >= cor_rad[..., None], axis=-1)
    mask_pn = radius <= pns_rad[..., None]
    ind_pn = np.argmax(radius[None, None, :] >= pns_rad[..., None], axis=-1)
    mask_in = radius <= inn_rad[..., None]
    ind_in = np.argmax(radius[None, None, :] >= inn_rad[..., None], axis=-1)
    mask_sho = radius <= sho_rad[..., None]
    ind_sho = np.argmax(radius[None, None, :] >= sho_rad[..., None], axis=-1)
    rho = simulation.rho(file_name)
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    vp = simulation.phi_velocity(file_name)
    rho_vr = rho * vr
    dmass = rho * dV
    corr = []
    Qradial = []
    Qfull = []
    Qcor = []
    Qpns = []
    Qinn = []
    Qsho = []
    Qpne = []
    Qsoe = []
    Qopn = []
    Qshk = []
    for i in range(5):
        Qdot = dmass * (
            vr * gradient[i][..., 0]
            + vt * gradient[i][..., 1]
            + vp * gradient[i][..., 2]
        )
        Qcorr = rho_vr * Ylm[i]
        Qradial.append(np.sum(Qdot, axis=(0, 1))[..., None])
        corr.append(np.sum(Qcorr, axis=(0, 1))[..., None])
        Qfull.append(np.sum(Qdot)[..., None])
        Qcor.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                None,
                mask_co,
                None,
                ind_co,
            )
        )
        Qpns.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                None,
                mask_pn,
                None,
                ind_pn,
            )
        )
        Qinn.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                None,
                mask_in,
                None,
                ind_in,
            )
        )
        Qsho.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                None,
                mask_sho,
                None,
                ind_sho,
            )
        )
        Qpne.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                mask_co,
                mask_pn,
                ind_co,
                ind_pn,
            )
        )
        Qsoe.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                mask_co,
                mask_in,
                ind_co,
                ind_in,
            )
        )
        Qopn.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                mask_pn,
                mask_in,
                ind_pn,
                ind_in,
            )
        )
        Qshk.append(
            __compute_3D_partial_Q(
                Qdot,
                Qcorr,
                mask_pn,
                mask_sho,
                ind_pn,
                ind_sho,
            )
        )
    gws_dictionary["matter_strain"].append(
        np.concatenate(
            Qradial,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["full_hpleq"].append(
        np.concatenate(
            Qfull,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["corr"].append(
        np.concatenate(
            corr,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["core_hpleq"].append(
        np.concatenate(
            Qcor,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["pns_hpleq"].append(
        np.concatenate(
            Qpns,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["innercore_hpleq"].append(
        np.concatenate(
            Qinn,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["shock_hpleq"].append(
        np.concatenate(
            Qsho,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["pns_env_hpleq"].append(
        np.concatenate(
            Qpne,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["sonic_env_hpleq"].append(
        np.concatenate(
            Qsoe,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["out_pns_hpleq"].append(
        np.concatenate(
            Qopn,
            axis=-1,
        )[..., None]
    )
    gws_dictionary["shocked_hpleq"].append(
        np.concatenate(
            Qshk,
            axis=-1,
        )[..., None]
    )


def __compute_3D_surface_Q_correction(
    correction: aerray,
    r1_ind: int,
    r2_ind: int,
) -> aerray:
    th_ind = np.arange(correction.shape[0])[:, None]
    ph_ind = np.arange(correction.shape[1])[None, :]
    if r1_ind is None:
        r1 = 0 * correction.unit
    else:
        r1 = correction[
            th_ind,
            ph_ind,
            r1_ind,
        ]
        r1 = r1.sum()
    if r2_ind is None:
        r2 = 0 * correction.unit
    else:
        r2 = correction[
            th_ind,
            ph_ind,
            r2_ind,
        ]
        r2 = r2.sum()
    return r2 - r1


def __compute_3D_partial_Q(
    Qdot: aerray,
    correction: aerray,
    mask1: np.array[np.bool],
    mask2: np.array[np.bool],
    r1_ind: np.array[np.int32],
    r2_ind: np.array[np.int32],
) -> aerray:
    if mask1 is None:
        mask = mask2
    else:
        mask = mask2 & np.logical_not(mask1)
    corr = __compute_3D_surface_Q_correction(
        correction,
        r1_ind,
        r2_ind,
    )
    corrected_Q = np.sum(Qdot[mask]) - corr
    return corrected_Q[..., None]
