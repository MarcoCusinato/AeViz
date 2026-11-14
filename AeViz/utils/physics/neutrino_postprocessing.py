from __future__ import annotations

from copy import deepcopy

import numpy as np

from AeViz.simulation import Simulation
from AeViz.units import aerray, u
from AeViz.units.constants import constants as c
from AeViz.utils.physics.radii_utils import neutrino_sphere_radii


def declare_neutrino_dictionary(dim: int, angles: list[int] | None = None):
    """
    Declare the neutrino post-processing dictionary.

    Parameters
    ----------
    dim : int
        Simulation dimension (1, 2 or 3).
    angles : list[int] | None
        Polar angles (degrees) of Leq_{a}, Llab_{a}, Elab_{a}, E2lab_{a}.
        Must be None in 1D and given in 2D/3D.

    Returns
    -------
    dict
        {'nue': {...}, 'nua': {...}, 'nux': {...}}; every leaf is a list.
        'nux' is the sum of the four heavy-lepton species.
        acc* and j* are rates acting on the PNS (cm/s^2, g cm^2/s^2).
        l2pl, l2cr are in cm^2/s^2 sr (Eq. above, B2).
    """
    nu_spec = {
        "r": [],
        "rmin": [],
        "rmax": [],
        "ravg": [],
        "Llab_avg": [],
        "Elab_avg": [],
        "E2lab_avg": [],
    }

    if dim == 1:
        nu_spec["Leq"] = []
    if dim > 1:
        angl = 0
        for i in range(len(angles)):
            nu_spec[f"Leq_{angl}"] = []
            nu_spec[f"Llab_{angl}"] = []
            nu_spec[f"Elab_{angl}"] = []
            nu_spec[f"E2lab_{angl}"] = []
            angl += 10
        nu_spec["Llab_pol"] = []
        nu_spec["Elab_pol"] = []
        nu_spec["E2lab_pol"] = []
        nu_spec["Llab_eq"] = []
        nu_spec["Elab_eq"] = []
        nu_spec["E2lab_eq"] = []
        nu_spec["accx"] = []
        nu_spec["accy"] = []
        nu_spec["accz"] = []
        nu_spec["jx"] = []
        nu_spec["jy"] = []
        nu_spec["jz"] = []
        nu_spec["hdot_eqpl"] = []
    if dim > 2:
        nu_spec["hdot_eqcr"] = []
        nu_spec["hdot_polpl"] = []
        nu_spec["hdot_polcr"] = []
        nu_spec["l2pl"] = []
        nu_spec["l2cr"] = []
    neutrino_dictionary = {nu: deepcopy(nu_spec) for nu in ("nue", "nua", "nux")}
    return neutrino_dictionary


def neutrino_quantity_metadata(angles: list[int] | None = None) -> dict:
    """
    Metadata for neutrino post-processing quantities.

    The same metadata applies to every species. The token 'NU' in a label
    must be replaced by the species symbol in meta['species'][nu]['tex'].
    Angle-resolved keys (Llab_{a}, ...) are generated from the templates.
    """

    nu_labels = {"nue": r"\nu_e", "nua": r"\overline{\nu}_e", "nux": r"\nu_x"}
    meta = {}
    for nu in ("nue", "nua", "nux"):
        lb = nu_labels[nu]
        meta[nu] = {
            "rmin": {
                "name": f"neutrinosphere_radius_min_{nu}",
                "label": r"$R_{" + lb + r"}^{\rm min}$",
                "limits": [0, 1e2],
                "cmap": "plasma",
                "log": False,
            },
            "rmax": {
                "name": f"neutrinosphere_radius_max_{nu}",
                "label": r"$R_{" + lb + r"}^{\rm max}$",
                "limits": [0, 1e2],
                "cmap": "plasma",
                "log": False,
            },
            "ravg": {
                "name": f"neutrinosphere_radius_avg_{nu}",
                "label": r"$\langle R_{" + lb + r"}\rangle$",
                "limits": [0, 1e2],
                "cmap": "plasma",
                "log": False,
            },
            "r": {
                "name": f"neutrinosphere_radius_surface_{nu}",
                "label": r"$R^{" + lb + r"}(\theta,\phi)$",
                "limits": [5, 1e2],
                "cmap": "plasma",
                "log": False,
            },
            "Leq": {
                "name": f"equivalent_luminosity_{nu}",
                "label": r"$L_{" + lb + r"}^{\rm eq}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            },
            "Llab_avg": {
                "name": f"lab_luminosity_average_{nu}",
                "label": r"$L_{" + lb + r"}^{\rm lab}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            },
            "Elab_avg": {
                "name": f"lab_mean_energy_average_{nu}",
                "label": r"$\langle E^{" + lb + r"}\rangle_{\rm lab}$",
                "limits": [0, 30],
                "cmap": "viridis",
                "log": False,
            },
            "E2lab_avg": {
                "name": f"lab_mean_square_energy_average_{nu}",
                "label": r"$\langle E^{2," + lb + r"}\rangle_{\rm lab}$",
                "limits": [0, 200],
                "cmap": "viridis",
                "log": False,
            },
            "Llab_pol": {
                "name": f"lab_luminosity_polar_{nu}",
                "label": r"$L_{" + lb + r"}^{\rm lab,pol}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            },
            "Elab_pol": {
                "name": f"lab_mean_energy_polar_{nu}",
                "label": r"$\langle E^{" + lb + r"}\rangle_{\rm lab,pol}$",
                "limits": [0, 30],
                "cmap": "viridis",
                "log": False,
            },
            "E2lab_pol": {
                "name": f"lab_mean_square_energy_polar_{nu}",
                "label": r"$\langle E^{2," + lb + r"}\rangle_{\rm lab,pol}$",
                "limits": [0, 200],
                "cmap": "viridis",
                "log": False,
            },
            "Llab_eq": {
                "name": f"lab_luminosity_equatorial_{nu}",
                "label": r"$L_{" + lb + r"}^{\rm lab,eq}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            },
            "Elab_eq": {
                "name": f"lab_mean_energy_equatorial_{nu}",
                "label": r"$\langle E^{" + lb + r"}\rangle_{\rm lab,eq}$",
                "limits": [0, 30],
                "cmap": "viridis",
                "log": False,
            },
            "E2lab_eq": {
                "name": f"lab_mean_square_energy_equatorial_{nu}",
                "label": r"$\langle E^{2," + lb + r"}\rangle_{\rm lab,eq}$",
                "limits": [0, 200],
                "cmap": "viridis",
                "log": False,
            },
            "accx": {
                "name": f"neutrino_kick_acceleration_x_{nu}",
                "label": r"$a_{" + lb + r"}^{x}$",
                "limits": [-1e8, 1e8],
                "cmap": "coolwarm",
                "log": False,
            },
            "accy": {
                "name": f"neutrino_kick_acceleration_y_{nu}",
                "label": r"$a_{" + lb + r"}^{y}$",
                "limits": [-1e8, 1e8],
                "cmap": "coolwarm",
                "log": False,
            },
            "accz": {
                "name": f"neutrino_kick_acceleration_z_{nu}",
                "label": r"$a_{" + lb + r"}^{z}$",
                "limits": [-1e8, 1e8],
                "cmap": "coolwarm",
                "log": False,
            },
            "jx": {
                "name": f"neutrino_torque_x_{nu}",
                "label": r"$\dot{J}_{" + lb + r"}^{x}$",
                "limits": [-1e49, 1e49],
                "cmap": "coolwarm",
                "log": False,
            },
            "jy": {
                "name": f"neutrino_torque_y_{nu}",
                "label": r"$\dot{J}_{" + lb + r"}^{y}$",
                "limits": [-1e49, 1e49],
                "cmap": "coolwarm",
                "log": False,
            },
            "jz": {
                "name": f"neutrino_torque_z_{nu}",
                "label": r"$\dot{J}_{" + lb + r"}^{z}$",
                "limits": [-1e49, 1e49],
                "cmap": "coolwarm",
                "log": False,
            },
            "hdot_eqpl": {
                "name": f"neutrino_gw_rate_equatorial_plus_{nu}",
                "label": r"$\dot{(rh)}_{" + lb + r"}^{+,\rm eq}$",
                "limits": [-2e3, 2e3],
                "cmap": "seismic",
                "log": False,
            },
            "hdot_eqcr": {
                "name": f"neutrino_gw_rate_equatorial_cross_{nu}",
                "label": r"$\dot{(rh)}_{" + lb + r"}^{\times,\rm eq}$",
                "limits": [-2e3, 2e3],
                "cmap": "seismic",
                "log": False,
            },
            "hdot_polpl": {
                "name": f"neutrino_gw_rate_polar_plus_{nu}",
                "label": r"$\dot{(rh)}_{" + lb + r"}^{+,\rm pol}$",
                "limits": [-2e3, 2e3],
                "cmap": "seismic",
                "log": False,
            },
            "hdot_polcr": {
                "name": f"neutrino_gw_rate_polar_cross_{nu}",
                "label": r"$\dot{(rh)}_{" + lb + r"}^{\times,\rm pol}$",
                "limits": [-2e3, 2e3],
                "cmap": "seismic",
                "log": False,
            },
            "l2pl": {
                "name": f"neutrino_gw_l2_plus_{nu}",
                "label": r"$\mathcal{L}^{2," + lb + r"}_{+}$",
                "limits": [1e100, 1e106],
                "cmap": "Oranges",
                "log": True,
            },
            "l2cr": {
                "name": f"neutrino_gw_l2_cross_{nu}",
                "label": r"$\mathcal{L}^{2," + lb + r"}_{\times}$",
                "limits": [1e100, 1e106],
                "cmap": "Oranges",
                "log": True,
            },
        }
        angl = 0
        for i in range(len(angles)):
            meta[nu][f"Leq_{angl}"] = {
                "name": f"equivalent_luminosity_{nu}_{angl}",
                "label": r"$L_{" + lb + r"}^{\rm eq," + str(angl) + r"^\circ}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            }
            meta[nu][f"Llab_{angl}"] = {
                "name": f"lab_luminosity_polar_{nu}_{angl}",
                "label": r"$L_{" + lb + r"}^{\rm lab," + str(angl) + r"^\circ}$",
                "limits": [0, 1e53],
                "cmap": "inferno",
                "log": False,
            }
            meta[nu][f"Elab_{angl}"] = {
                "name": f"lab_mean_energy_equatorial_{nu}_{angl}",
                "label": r"$\langle E^{"
                + lb
                + r"}\rangle_{\rm lab,"
                + str(angl)
                + r"^\circ}$",
                "limits": [0, 30],
                "cmap": "viridis",
                "log": False,
            }
            meta[nu][f"E2lab_{angl}"] = {
                "name": f"lab_mean_square_energy_equatorial_{nu}_{angl}",
                "label": r"$\langle E^{2,"
                + lb
                + r"}\rangle_{\rm lab,"
                + str(angl)
                + r"^\circ}$",
                "limits": [0, 200],
                "cmap": "viridis",
                "log": False,
            }
            angl += 10
    return meta


def compute_neutrino_postprocessing(
    neutrino_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    dOmega: aerray,
    Wh: tuple[tuple[aerray]] | None,
    dOmega_obs: aerray | None,
    index500: int,
    index_theta: int | None,
    indices_theta_avg: tuple[int] | None,
    index_phi: int | None,
    PNS_mass: aerray | None,
) -> None:
    __reduced_neutrino_output(
        neutrino_dictionary, simulation, file_name, index500, index_theta, index_phi
    )
    __full_neutrino_output(
        neutrino_dictionary,
        simulation,
        file_name,
        index500,
        index_theta,
        indices_theta_avg,
        dOmega,
        index_phi,
        Wh,
        dOmega_obs,
        PNS_mass,
    )
    return True


def __full_neutrino_output(
    neutrino_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    index500: int,
    indexTheta: int | None,
    indicesThetaAvg: tuple[int] | None,
    dOmega: aerray,
    indexPhi: int | None,
    Wh: tuple[tuple[aerray]] | None,
    dOmega_obs: aerray | None,
    PNS_mass: aerray | None,
) -> None:
    ene = simulation.cell.E_nu()
    while ene.ndim < simulation.dim:
        ene = ene[None, ...]
    nuer, nuar, nuxr = neutrino_sphere_radii(simulation, file_name)
    neutrino_dictionary["nue"]["r"].append(nuer[..., None])
    neutrino_dictionary["nue"]["rmin"].append(np.nanmin(nuer))
    neutrino_dictionary["nue"]["rmax"].append(np.nanmax(nuer))
    neutrino_dictionary["nue"]["ravg"].append(
        simulation.grid.radii_angular_average(nuer)
    )
    neutrino_dictionary["nua"]["r"].append(nuar[..., None])
    neutrino_dictionary["nua"]["rmin"].append(np.nanmin(nuar))
    neutrino_dictionary["nua"]["rmax"].append(np.nanmax(nuar))
    neutrino_dictionary["nua"]["ravg"].append(
        simulation.grid.radii_angular_average(nuar)
    )
    neutrino_dictionary["nux"]["r"].append(nuxr[..., None])
    neutrino_dictionary["nux"]["rmin"].append(np.nanmin(nuxr))
    neutrino_dictionary["nux"]["rmax"].append(np.nanmax(nuxr))
    neutrino_dictionary["nux"]["ravg"].append(
        simulation.grid.radii_angular_average(nuxr)
    )
    Lnu_lab_frame = simulation.labframe_neutrino_fluxes(file_name)
    Lnu_lab_frame = [
        Lnu_lab_frame[i] * simulation.grid.radius[index500] ** 2 for i in range(3)
    ]
    for inu, nu in enumerate(["nue", "nua", "nux"]):
        lumlab = (
            Lnu_lab_frame[inu][..., index500, :]
            if simulation.dim == 1
            else Lnu_lab_frame[inu][..., index500, :, 0]
        )
        num_nu = np.nansum(lumlab / ene, axis=-1)
        lnu = np.nansum(lumlab, axis=-1)
        l2nu = np.nansum(lumlab * ene, axis=-1)
        if simulation.dim == 1:
            neutrino_dictionary[nu]["Llab_avg"].append(
                (4 * np.pi * lnu).to(u.erg / u.s)
            )
            neutrino_dictionary[nu]["Elab_avg"].append((lnu / num_nu).to(u.MeV))
            neutrino_dictionary[nu]["E2lab_avg"].append((l2nu / num_nu).to(u.MeV**2))
        else:
            angle = 0
            for it in indexTheta:
                if simulation.dim == 3:
                    neutrino_dictionary[nu][f"Llab_{angle}"].append(
                        (4 * np.pi * lnu[indexPhi, it]).to(u.erg / u.s)
                    )
                    neutrino_dictionary[nu][f"Elab_{angle}"].append(
                        (lnu[indexPhi, it] / num_nu[indexPhi, it]).to(u.MeV)
                    )
                    neutrino_dictionary[nu][f"E2lab_{angle}"].append(
                        (l2nu[indexPhi, it] / num_nu[indexPhi, it]).to(u.MeV**2)
                    )
                else:
                    neutrino_dictionary[nu][f"Llab_{angle}"].append(
                        (4 * np.pi * lnu[it]).to(u.erg / u.s)
                    )
                    neutrino_dictionary[nu][f"Elab_{angle}"].append(
                        (lnu[it] / num_nu[it]).to(u.MeV)
                    )
                    neutrino_dictionary[nu][f"E2lab_{angle}"].append(
                        (l2nu[it] / num_nu[it]).to(u.MeV**2)
                    )
                angle += 10
            num_nu *= dOmega
            lnu *= dOmega
            l2nu *= dOmega
            j1, j2, j3 = indicesThetaAvg
            neutrino_dictionary[nu]["Llab_avg"].append(np.nansum(lnu).to(u.erg / u.s))
            neutrino_dictionary[nu]["Elab_avg"].append(
                (np.nansum(lnu) / np.nansum(num_nu)).to(u.MeV)
            )
            neutrino_dictionary[nu]["E2lab_avg"].append(
                (np.nansum(l2nu) / np.nansum(num_nu)).to(u.MeV**2)
            )
            neutrino_dictionary[nu]["Llab_pol"].append(
                (4 * np.pi * np.nansum(lnu[..., :j1]) / np.nansum(dOmega[..., :j1])).to(
                    u.erg / u.s
                )
            )
            neutrino_dictionary[nu]["Elab_pol"].append(
                (np.nansum(lnu[..., :j1]) / np.nansum(num_nu[..., :j1])).to(u.MeV)
            )
            neutrino_dictionary[nu]["E2lab_pol"].append(
                (np.nansum(l2nu[..., :j1]) / np.nansum(num_nu[..., :j1])).to(u.MeV**2)
            )
            neutrino_dictionary[nu]["Llab_eq"].append(
                (
                    4
                    * np.pi
                    * np.nansum(lnu[..., j2:j3])
                    / np.nansum(dOmega[..., j2:j3])
                ).to(u.erg / u.s)
            )
            neutrino_dictionary[nu]["Elab_eq"].append(
                (np.nansum(lnu[..., j2:j3]) / np.nansum(num_nu[..., j2:j3])).to(u.MeV)
            )
            neutrino_dictionary[nu]["E2lab_eq"].append(
                (np.nansum(l2nu[..., j2:j3]) / np.nansum(num_nu[..., j2:j3])).to(
                    u.MeV**2
                )
            )
            __compute_kick_ang_mom(
                neutrino_dictionary[nu],
                simulation,
                np.nansum(Lnu_lab_frame[inu], axis=-2)[..., index500, :],
                dOmega,
                PNS_mass,
                simulation.grid.radius[index500],
            )
            if simulation.dim == 2:
                ## lnu already multiplied by dOmega
                gws = __compute_neutrino_contribution_gws_2D(lnu, Wh)
                neutrino_dictionary[nu]["hdot_eqpl"].append(gws)
            elif simulation.dim == 3:
                gweq, gwpol, gwene = __compute_neutrino_contribution_gws_3D(
                    lnu, Wh, dOmega_obs
                )
                neutrino_dictionary[nu]["hdot_eqpl"].append(gweq[0])
                neutrino_dictionary[nu]["hdot_eqcr"].append(gweq[1])
                neutrino_dictionary[nu]["hdot_polpl"].append(gwpol[0])
                neutrino_dictionary[nu]["hdot_polcr"].append(gwpol[1])
                neutrino_dictionary[nu]["l2pl"].append(gwene[0])
                neutrino_dictionary[nu]["l2cr"].append(gwene[1])


def __reduced_neutrino_output(
    neutrino_dictionary: dict,
    simulation: Simulation,
    file_name: str,
    index500: int,
    indexTheta: int | None,
    indexPhi: int | None,
) -> None:
    equivalent_factor = np.pi * 4 * simulation.grid.radius[index500] ** 2
    nulums = simulation.neutrino_luminosity(file_name)
    for fl, lum in zip(["nue", "nua", "nux"], nulums):
        if indexPhi is not None:
            lum = lum[indexPhi, :, :]
        if indexTheta is None:
            neutrino_dictionary[fl]["Leq"].append(equivalent_factor * lum[index500])
        else:
            angle = 0
            for it in indexTheta:
                neutrino_dictionary[fl][f"Leq_{angle}"].append(
                    equivalent_factor * lum[it, index500]
                )
                angle += 10


def __compute_neutrino_contribution_gws_2D(neuflux: aerray, wh: aerray):
    wh = wh[0][0]
    gw = (2 * c.G / (c.c**4) * np.nansum(wh * neuflux)).to(u.cm / u.s)
    return gw


def __compute_neutrino_contribution_gws_3D(
    neuflux: aerray, wh: tuple[tuple[aerray]], dOmega_obs: aerray
):
    weq, wpol, wtd = wh

    gwpleq = (2 * c.G / (c.c**4) * np.nansum(weq[0] * neuflux)).to(u.cm / u.s)
    gwcreq = (2 * c.G / (c.c**4) * np.nansum(weq[1] * neuflux)).to(u.cm / u.s)

    gwplpol = (2 * c.G / (c.c**4) * np.nansum(wpol[0] * neuflux)).to(u.cm / u.s)
    gwcrpol = (2 * c.G / (c.c**4) * np.nansum(wpol[1] * neuflux)).to(u.cm / u.s)

    l2pl = np.nansum(dOmega_obs * (wtd[0] * neuflux[None, None, ...]) ** 2)
    l2cr = np.nansum(dOmega_obs * (wtd[1] * neuflux[None, None, ...]) ** 2)
    return ((gwpleq, gwcreq), (gwplpol, gwcrpol), (l2pl, l2cr))


def __compute_kick_ang_mom(
    neutrino_dictionary: dict,
    simulation: Simulation,
    nufluxes: aerray,
    dOmega: aerray,
    PNS_mass: aerray,
    r500: aerray,
):
    nufluxes = nufluxes * dOmega[..., None]
    ## taking only the radial flux since we are far enogh and the neutrinos are free-streaming
    fx, fy, fz = simulation.grid.spherical_to_cartesian(
        nufluxes[..., 0][..., None],
        np.zeros_like(nufluxes[..., 0][..., None]),
        np.zeros_like(nufluxes[..., 0][..., None]),
    )
    accz = (-np.nansum(fz) / (PNS_mass * c.c)).to(u.cm / u.s**2)
    accx = (
        (-np.nansum(fx) / (PNS_mass * c.c)).to(u.cm / u.s**2)
        if simulation.dim == 3
        else 0 * accz.unit
    )
    accy = (
        (-np.nansum(fy) / (PNS_mass * c.c)).to(u.cm / u.s**2)
        if simulation.dim == 3
        else 0 * accz.unit
    )

    neutrino_dictionary["accx"].append(accx.to(u.cm / u.s**2))
    neutrino_dictionary["accy"].append(accy.to(u.cm / u.s**2))
    neutrino_dictionary["accz"].append(accz.to(u.cm / u.s**2))

    ## compute the neutrino contribution to the ang mom
    ## as 0, +fphi, -ftheta
    jx, jy, jz = simulation.grid.spherical_to_cartesian(
        np.zeros_like(nufluxes[..., 0][..., None]),
        nufluxes[..., 2][..., None],
        -nufluxes[..., 1][..., None],
    )
    jz = r500 / c.c * np.nansum(jz)
    jx = r500 / c.c * np.nansum(jx) if simulation.dim == 3 else 0 * jz.unit
    jy = r500 / c.c * np.nansum(jy) if simulation.dim == 3 else 0 * jz.unit
    neutrino_dictionary["jx"].append(jx)
    neutrino_dictionary["jy"].append(jy)
    neutrino_dictionary["jz"].append(jz)


def compute_weights(
    theta: aerray, phi: aerray, alpha: aerray, beta: aerray, dim: int
) -> tuple[aerray]:
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
        N = (
            (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha))
            * np.sin(theta)
            * np.cos(beta)
            - np.cos(theta) * np.sin(beta)
        ) ** 2 + np.sin(theta) ** 2 * (
            np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)
        ) ** 2
        Dpl = (
            1
            + (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha))
            * np.sin(theta)
            * np.sin(beta)
            + np.cos(theta) * np.cos(beta)
        ) * (
            (
                (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha))
                * np.sin(theta)
                * np.cos(beta)
                - np.cos(theta) * np.sin(beta)
            )
            ** 2
            - np.sin(theta) ** 2
            * (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha)) ** 2
        )
        Dcr = (
            (
                1
                + (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha))
                * np.sin(theta)
                * np.sin(beta)
                + np.cos(theta) * np.cos(beta)
            )
            * 2
            * (
                (np.cos(phi) * np.cos(alpha) + np.sin(phi) * np.sin(alpha))
                * np.sin(theta)
                * np.cos(beta)
                - np.cos(theta) * np.sin(beta)
            )
            * np.sin(theta)
            * (np.sin(phi) * np.cos(alpha) - np.cos(phi) * np.sin(alpha))
        )
        whpl = Dpl / N
        whcr = Dcr / N
    return (whpl, whcr)


def build_obs_grid(nalpha: int):
    assert nalpha % 2 == 0, "Number of alphas must be even."
    dangle = 2 * np.pi / nalpha * u.radian
    beta = (np.arange(int(nalpha / 2)) + 0.5) * dangle
    alpha = (np.arange(int(nalpha)) + 0.5) * dangle
    beta_sinb = np.sin(beta)[:, None, None, None] * dangle
    dOmega_obs = (np.ones(int(nalpha / 2))[None, :, None, None] * dangle) * beta_sinb
    return dOmega_obs, alpha, beta
