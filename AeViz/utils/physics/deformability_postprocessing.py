from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import Akima1DInterpolator

from AeViz.simulation import Simulation
from AeViz.units import aerray, u
from AeViz.units.constants import constants as c

"""
The procedure follows what is described in
https://doi.org/10.1103/PhysRevD.81.123016
The linearised metric follows equation 12
The Love number follows equation 14
"""


def declare_tidal_dictionary() -> dict:
    tidal_dictionary = {
        "core": {"tidal_deformability": [], "kappa2": []},
        "PNS": {"tidal_deformability": [], "kappa2": []},
    }
    return tidal_dictionary


def tidal_quantity_metadata() -> dict:
    meta = {
        "core": {
            "tidal_deformability": {
                "name": "tidal_deformability_core",
                "label": r"$\Lambda^\mathrm{core}$",
                "limits": [1, 1e4],
                "log": True,
                "cmap": "cividis",
            },
            "kappa2": {
                "name": "kappa2_core",
                "label": r"$k_2^\mathrm{core}$",
                "limits": [0, 0.2],
                "log": False,
                "cmap": "cividis",
            },
        },
        "PNS": {
            "tidal_deformability": {
                "name": "tidal_deformability_PNS",
                "label": r"$\Lambda^\mathrm{PNS}$",
                "limits": [1, 1e4],
                "log": True,
                "cmap": "cividis",
            },
            "kappa2": {
                "name": "kappa2_PNS",
                "label": r"$k_2^\mathrm{PNS}$",
                "limits": [0, 0.2],
                "log": False,
                "cmap": "cividis",
            },
        },
    }
    return meta


def compute_tidal_love(
    tidal_dictionary: dict,
    simulation: Simulation,
    density: aerray,
    pressure: aerray,
    soundspeed: aerray,
    entropy: aerray,
):
    ind_cor = np.argmax(entropy >= (4 * u.kBol / u.bry))
    ind_pns = np.argmax(density >= (1e11 * u.g / u.cm**3))
    mass = (
        np.cumsum(
            4 * np.pi * simulation.cell.dr_integration(simulation.ghost) * density
        ).to(u.Msun)
        / u.Msun
    )
    geo_dens = __to_geometrised_dens(density)
    geo_press = __to_geometrised_pres(pressure)
    geo_cs = __to_geometrised_vel(soundspeed)
    radius = __to_geometrised_len(simulation.radius)
    l_cr, k2_cr = __solve_tidal_love(
        mass[ind_cor] / radius[ind_cor],
        geo_press[: ind_cor + 1],
        mass[: ind_cor + 1],
        geo_cs[: ind_cor + 1],
        geo_dens[: ind_cor + 1],
        radius[: ind_cor + 1],
    )
    l_pn, k2_pn = __solve_tidal_love(
        mass[ind_pns] / radius[ind_pns],
        geo_press[: ind_pns + 1],
        mass[: ind_pns + 1],
        geo_cs[: ind_pns + 1],
        geo_dens[: ind_pns + 1],
        radius[: ind_pns + 1],
    )
    tidal_dictionary["core/kappa2"].append(k2_cr)
    tidal_dictionary["core/tidal_deformability"].append(l_cr)
    tidal_dictionary["PNS/kappa2"].append(k2_pn)
    tidal_dictionary["PNS/tidal_deformability"].append(l_pn)


def __to_geometrised_len(x: aerray) -> aerray:
    """
    Convert the radius to geometrised units

    Parameters
    ----------
    x : aerray
        radius

    Returns
    -------
    aerray
    """
    return x / c.G / c.Msol * c.c**2


def __to_geometrised_dens(x: aerray) -> aerray:
    """
    Convert the density to geometrised units

    Parameters
    ----------
    x : aerray
        density

    Returns
    -------
    aerray
    """
    return x / c.Msol * c.G**3 * c.Msol**3 / c.c**6


def __to_geometrised_vel(x: aerray) -> aerray:
    """
    Convert the velocity to geometrised units

    Parameters
    ----------
    x : aerray
        velocity

    Returns
    -------
    aerray
    """
    return x / c.c


def __to_geometrised_pres(x: aerray) -> aerray:
    """
    Convert the pressure to geometrised units

    Parameters
    ----------
    x : aerray
        pressure

    Returns
    -------
    aerray
    """
    return __to_geometrised_dens(x) / c.c**2


def __linearised_metric(
    r: aerray, y: list[aerray], p: aerray, m: aerray, cs: aerray, rho: aerray
) -> list[aerray]:
    """
    Computes the linearised metric as in equation 12 of Hinderer et al,
    2010, https://doi.org/10.1103/PhysRevD.81.123016.

    .. math::
        \frac{d\beta}{dr}
            =
            2\\left(1-\frac{2m_r}{r}\right)^{-1} H
            \\left\\{
            -2\\pi \\left[5\\epsilon + 9p + f(\\epsilon + p)\right]
            + \frac{3}{r^2}
            + 2\\left(1-\frac{2m_r}{r}\right)^{-1}
            \\left(\frac{m_r}{r^2}+4\\pi r p\right)^2
            \right\\}
            +\frac{2\beta}{r}
            \\left(1-\frac{2m_r}{r}\right)^{-1}
            \\left\\{
            -1+\frac{m_r}{r}+2\\pi r^2(\\epsilon-p)
            \right\\}.

    Parameters
    ----------
    r : aerray
        radius
    y : list[aerray]
        .. math::
            \frac{dH}{dr} \\qquad H
    p : aerray
        pressure
    m : aerray
        mass
    cs : aerray
        soundspeed
    rho : aerray
        density

    Returns
    -------
    list[aerray]
        Returns the first and second dericative of H
    """
    dH_dr, H = y
    if r == 0:
        return [y[0], 0]

    coeff = 1 / (1 - 2 * m(r) / r)

    d2H_dr2 = 2 * coeff * H * (
        -2 * np.pi * (5 * rho(r) + 9 * p(r) + (rho(r) + p(r)) / cs(r) ** 2)
        + 3 / r**2
        + 2 * coeff * (m(r) / r**2 + 4 * np.pi * r * p(r)) ** 4
    ) + 2 * dH_dr / r * coeff * (-1 + m(r) / r + 2 * np.pi * r**2 * (rho(r) - p(r)))
    return [d2H_dr2, dH_dr]


def __compute_Love_number(
    xi: aerray,
    y: aerray,
) -> aerray:
    """
    Compute the love number for a given compacness and
    R H'(R)/H(R) = y, as of equation 14 of Hinderer et al 2010.
    .. math::
        k_2 = \frac{8C^5}{5(1-2C)^2}
            \\left[
            2 + 2C(y-1) - y
            \\left\\{
            2C\\left[6 - 3y + 3C(5y-8)\right]
            + 4C^3\\left[13 - 11y + C(3y-2) + 2C^2(1+y)\right]
            + 3(1-2C)^2\\left[2 - y + 2C(y-1)\right]\\ln(1-2C)
            \right\\}
            \right]^{-1}

    Parameters
    ----------
    xi : aerray
        compactness
    y : aerray
        .. math::
            R \frac{dH}{dr}\frac{1}{H}

    Returns
    -------
    aerray
        returns the tidal love number
    """
    XI = 1 - 2 * xi
    kappa2 = (
        8
        / 5
        * xi**5
        * XI**2
        * (2 + 2 * xi * (y - 1) - y)
        / (
            2 * xi * (6 - 3 * y + 3 * xi * (5 * y - 8))
            + 4 * xi**3 * (13 - 11 * y + xi * (3 * y - 2) + 2 * xi**2 * (1 + y))
            + 3 * XI**2 * (2 - y + 2 * xi * (y - 1)) * (np.log(XI))
        )
    )
    return kappa2


def __compute_tidal_deformability(
    kappa2: aerray,
    xi: aerray,
) -> aerray:
    """
    Compute the tidal deformability from the love number

    Parameters
    ----------
    kappa2 : aerray
        tidal love number
    xi : aerray
        compactness

    Returns
    -------
    aerray
        tidal deformability of the star
    """
    return 2 / 3 * kappa2 / xi**5


def __solve_tidal_love(
    xi: aerray,
    pres: aerray,
    mass: aerray,
    soundspeed: aerray,
    dens: aerray,
    radius: aerray,
) -> list[aerray]:
    """
    Computes the tidal Love number for a given compactness xi

    Parameters
    ----------
    xi : aerray
        compactness
    pres : aerray
        pressure
    mass : aerray
        enclosed mass
    soundspeed : aerray
        speed of sound of the matter
    dens : aerray
        density
    radius : aerray
        radius

    Returns
    -------
    list[aerray]
        list contining the tidal love number and tidal deformability
    """
    pres = Akima1DInterpolator(radius, pres, method="akima", extrapolate=True)
    dens = Akima1DInterpolator(radius, dens, method="akima", extrapolate=True)
    mass = Akima1DInterpolator(radius, mass, method="akima", extrapolate=True)
    soundspeed = Akima1DInterpolator(
        radius, soundspeed, method="akima", extrapolate=True
    )
    solution = solve_ivp(
        __linearised_metric,
        [radius[0], radius[-1]],
        [2 * radius[0], radius[0] ** 2],
        method="LSODA",
        t_eval=radius,
        args=[pres, mass, soundspeed, dens],
        max_step=10000,
        atol=1e-9,
        rtol=1e-5,
    )
    y = solution.t[-1] * solution.y[0][-1] / solution.y[1][-1]
    kappa2 = __compute_Love_number(xi, y)
    tidal_d = __compute_tidal_deformability(kappa2, xi)
    return kappa2, tidal_d
