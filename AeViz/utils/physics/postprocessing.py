from __future__ import annotations

import warnings

import numpy as np
from tqdm import tqdm

from AeViz.simulation.simulation import Simulation
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
from AeViz.units import u
from AeViz.utils.files.file_utils import load_dataset, save_merge_dictionary_hdf
from AeViz.utils.physics.deformability_postprocessing import (
    compute_tidal_love,
    declare_tidal_dictionary,
    tidal_quantity_metadata,
)
from AeViz.utils.physics.gain_postprocessing import (
    compute_gain_postprocessing,
    declare_gain_dictionary,
    gain_quantity_metadata,
)
from AeViz.utils.physics.GW_postprocessing import (
    GWs_metadata,
    compute_gws_postprocessing,
    declare_GW_dictionary,
)
from AeViz.utils.physics.innercore_postprocessing import (
    compute_innercore_postprocessing,
    declare_innercore_dictionary,
    innercore_quantity_metadata,
)
from AeViz.utils.physics.neutrino_postprocessing import (
    build_obs_grid,
    compute_neutrino_postprocessing,
    compute_weights,
    declare_neutrino_dictionary,
    neutrino_quantity_metadata,
)
from AeViz.utils.physics.PNS_core_postprocessing import (
    PNScore_quantity_metadata,
    compute_PNScore_postprocessing,
    declare_PNScore_dictionary,
)
from AeViz.utils.physics.PNS_postprocessing import (
    PNS_quantity_metadata,
    compute_PNS_postprocessing,
    declare_PNS_dictionary,
)
from AeViz.utils.physics.profiles_postprocessing import (
    compute_profiles_postprocessing,
    create_profiles_metadata,
    declare_profiles_dictionary,
)
from AeViz.utils.physics.shock_postprocessing import (
    compute_shock_postprocessing,
    declare_shock_dictionary,
    shock_quantity_metadata,
)
from AeViz.utils.physics.spherical_harmonics_postprocessing import (
    compute_spherical_harmonics_decomposition,
    create_metadata_spherical_harmonics,
    declare_spherical_harmonics_dicionary,
)
from AeViz.utils.utils import check_existence, checkpoints


def supernova_postprocessing(
    simulation: Simulation, save_checkpoints: bool = True
) -> None:
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    file_list = simulation.hdf_file_list
    pns_processed = []
    if check_existence(simulation, "PNS_postprocessing.h5"):
        pns_processed = load_dataset(
            simulation.storage_path, "PNS_postprocessing.h5", "processed_hdf", False
        )
        if len(pns_processed) == len(file_list) and pns_processed[-1] == file_list[-1]:
            return
        tqdm.write(
            f"Checkpoint found for the Supernova Postprocessing, starting from {pns_processed[-1]}."
        )
    else:
        tqdm.write("Starting Supernova Postprocessing.")

    dV = simulation.cell.dVolume_integration(simulation.ghost)
    radius = simulation.cell.radius(simulation.ghost)
    grid = simulation.grid.cartesian_grid()
    dOmega = simulation.cell.dOmega(simulation.ghost)

    ## unpack the grid
    if simulation.dim == 1:
        X = grid * radius.unit
        rperp2 = radius**2
    elif simulation.dim == 2:
        rperp2 = (
            radius[None, :] * np.sin(simulation.cell.theta(simulation.ghost))[:, None]
        ) ** 2
        X, Y, Z = grid[0] * radius.unit, grid[1] * radius.unit, 0 * radius.unit
    elif simulation.dim == 3:
        rperp2 = (
            radius[None, None, :]
            * np.sin(simulation.cell.theta(simulation.ghost))[None, :, None]
        ) ** 2
        X, Y, Z = grid[0] * radius.unit, X[1] * radius.unit, grid[3] * radius.unit

    start_index = len(pns_processed)
    file_list = file_list[start_index:]
    ## Profiles declaration
    profiles_dict = declare_profiles_dictionary(
        simulation.dim, simulation.evolved_qts["magdim"]
    )
    profiles_metadata = None
    ## PNS core declaration part
    PNS_core_dict = declare_PNScore_dictionary(
        simulation.dim, simulation.evolved_qts["magdim"]
    )
    pnscore_metadata = PNScore_quantity_metadata()
    ## PNS declaration part
    PNS_dict = declare_PNS_dictionary(simulation.dim, simulation.evolved_qts["magdim"])
    pns_metadata = PNS_quantity_metadata()
    ## innercore declaration part
    innercore_dict = declare_innercore_dictionary(
        simulation.dim, simulation.evolved_qts["magdim"]
    )
    innercore_metadata = innercore_quantity_metadata()
    ## shock declaration part
    shock_dict = declare_shock_dictionary(
        simulation.dim, simulation.evolved_qts["magdim"]
    )
    shock_metadata = shock_quantity_metadata()
    ## gain declaration part
    gain_dict = declare_gain_dictionary(
        simulation.dim, simulation.evolved_qts["magdim"]
    )
    gain_metadata = gain_quantity_metadata()
    global_time = []
    check_index = 0
    pbar = tqdm(file_list, "Supernova postprocessing", colour="#008cff")
    for file_name in pbar:
        pbar.set_postfix(file=file_name)
        ## Compute the angular momentum
        dmass = simulation.rho(file_name) * dV
        time = simulation.time(file_name)
        ## Compute angular momentum
        if simulation.dim == 1:
            jx, jy, jz = None, None, None
            jtot = None
            shell_jx, shell_jy, shell_jz = None, None, None
            shell_jtot = None
        else:
            vr = simulation.radial_velocity(file_name)
            vt = simulation.theta_velocity(file_name)
            vp = simulation.phi_velocity(file_name)
            vx, vy, vz = simulation.grid.spherical_to_cartesian(
                dmass * vr, dmass * vt, dmass * vp
            )
            if simulation.dim == 2:
                jx, jy, jz = (
                    np.zeros(vz.shape) * vz.unit * X.unit,
                    np.zeros(vz.shape) * vz.unit * X.unit,
                    (vy * X - vx * Y),
                )
            else:
                jx, jy, jz = (vz * X - vy * Z), (vx * Z - vz * X), (vy * X - vx * Y)
            jtot = np.sqrt(jx**2 + jy**2 + jz**2)
            shell_jx = np.nansum(jx, axis=tuple(range(simulation.dim - 1)))
            shell_jy = np.nansum(jy, axis=tuple(range(simulation.dim - 1)))
            shell_jz = np.nansum(jz, axis=tuple(range(simulation.dim - 1)))
            shell_jtot = np.nansum(jtot, axis=tuple(range(simulation.dim - 1)))
        ## Compute inertia moment
        if simulation.dim < 3:
            Inertia = dmass * rperp2
            n, X_perp, eR = None, None, None
        else:
            nx = shell_jx / shell_jtot
            ny = shell_jy / shell_jtot
            nz = shell_jz / shell_jtot
            mm = np.isnan(nx) | np.isnan(ny) | np.isnan(nz)
            nx[mm] = 0.0
            ny[mm] = 0.0
            nz[mm] = 1.0
            n = np.array([nx, ny, nz])
            X_perp_vec = (
                grid
                - np.einsum("ir,iptr->ptr", n, grid)[None, ...] * n[:, None, None, :]
            )
            X_perp = np.linalg.norm(X_perp_vec, axis=0) + 1e-12
            eR = X_perp_vec / X_perp[None, ...]
            X_perp2 = np.where((n[2] >= 0.995)[None, None, :], rperp2, X_perp**2)
            Inertia = dmass * X_perp2
        compute_profiles_postprocessing(
            profiles_dict,
            simulation,
            file_name,
            shell_jx,
            shell_jy,
            shell_jz,
            shell_jtot,
            Inertia,
            n,
            X_perp,
            eR,
            radius.unit,
        )
        compute_PNS_postprocessing(
            PNS_dict, simulation, file_name, dV, dOmega, jx, jy, jz
        )
        compute_PNScore_postprocessing(
            PNS_core_dict, simulation, file_name, dV, dOmega, jx, jy, jz
        )
        compute_innercore_postprocessing(
            innercore_dict, simulation, file_name, dV, dOmega, jx, jy, jz
        )
        compute_shock_postprocessing(
            shock_dict,
            simulation,
            file_name,
            dV,
            jx,
            jy,
            jz,
            PNS_dict["local"]["r"][-1],
            innercore_dict["local"]["r"][-1],
        )
        compute_gain_postprocessing(
            gain_dict,
            simulation,
            file_name,
            dV,
            PNS_dict["local"]["r"][-1],
            shock_dict["local"]["r"][-1],
            jx,
            jy,
            jz,
        )
        global_time.append(time)
        pns_processed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            pbar.set_description("Saving.")
            if not profiles_metadata:
                profiles_metadata = create_profiles_metadata(profiles_dict)
            save_merge_dictionary_hdf(
                simulation,
                (global_time, profiles_dict, pns_processed),
                simulation.storage_path,
                "profiles_postprocessing.h5",
                None,
                **profiles_metadata,
            )
            save_merge_dictionary_hdf(
                simulation,
                (global_time, PNS_core_dict, pns_processed),
                simulation.storage_path,
                "PNScore_postprocessing.h5",
                ["global", "local"],
                **pnscore_metadata,
            )
            save_merge_dictionary_hdf(
                simulation,
                (global_time, PNS_dict, pns_processed),
                simulation.storage_path,
                "PNS_postprocessing.h5",
                ["global", "local"],
                **pns_metadata,
            )
            save_merge_dictionary_hdf(
                simulation,
                (global_time, innercore_dict, pns_processed),
                simulation.storage_path,
                "innercore_postprocessing.h5",
                ["global", "local"],
                **innercore_metadata,
            )
            save_merge_dictionary_hdf(
                simulation,
                (global_time, shock_dict, pns_processed),
                simulation.storage_path,
                "shock_postprocessing.h5",
                ["global", "local"],
                **shock_metadata,
            )
            save_merge_dictionary_hdf(
                simulation,
                (global_time, gain_dict, pns_processed),
                simulation.storage_path,
                "gain_postprocessing.h5",
                ["global", "local"],
                **gain_metadata,
            )
            check_index = 0
            profiles_dict = declare_profiles_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            PNS_core_dict = declare_PNScore_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            PNS_dict = declare_PNS_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            innercore_dict = declare_innercore_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            shock_dict = declare_shock_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            gain_dict = declare_gain_dictionary(
                simulation.dim, simulation.evolved_qts["magdim"]
            )
            global_time = []
            pbar.set_description("Supernova postprocessing")
        else:
            check_index += 1

    if not profiles_metadata:
        profiles_metadata = create_profiles_metadata(profiles_dict)
    save_merge_dictionary_hdf(
        simulation,
        (global_time, profiles_dict, pns_processed),
        simulation.storage_path,
        "profiles_postprocessing.h5",
        None,
        **profiles_metadata,
    )
    save_merge_dictionary_hdf(
        simulation,
        (global_time, PNS_core_dict, pns_processed),
        simulation.storage_path,
        "PNScore_postprocessing.h5",
        ["global", "local"],
        **pnscore_metadata,
    )
    save_merge_dictionary_hdf(
        simulation,
        (global_time, PNS_dict, pns_processed),
        simulation.storage_path,
        "PNS_postprocessing.h5",
        ["global", "local"],
        **pns_metadata,
    )
    save_merge_dictionary_hdf(
        simulation,
        (global_time, innercore_dict, pns_processed),
        simulation.storage_path,
        "innercore_postprocessing.h5",
        ["global", "local"],
        **innercore_metadata,
    )
    save_merge_dictionary_hdf(
        simulation,
        (global_time, shock_dict, pns_processed),
        simulation.storage_path,
        "shock_postprocessing.h5",
        ["global", "local"],
        **shock_metadata,
    )
    save_merge_dictionary_hdf(
        simulation,
        (global_time, gain_dict, pns_processed),
        simulation.storage_path,
        "gain_postprocessing.h5",
        ["global", "local"],
        **gain_metadata,
    )
    tqdm.write("Supernova postprocessing completed!")


def supernova_neutrino_postprocessing(
    simulation: Simulation, save_checkpoints: bool = True
) -> None:
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    file_list = simulation.hdf_file_list
    postprocessed = []
    if check_existence(simulation, "neutrino_postprocessing.h5"):
        postprocessed = load_dataset(
            simulation.storage_path,
            "neutrino_postprocessing.h5",
            "processed_hdf",
            False,
        )
        if len(postprocessed) == len(file_list) and postprocessed[-1] == file_list[-1]:
            return
        tqdm.write(
            f"Checkpoint found for the Neutrino Postprocessing, starting from {postprocessed[-1]}."
        )
    else:
        tqdm.write("Starting Neutrino Postprocessing.")
    if not check_existence(simulation, "PNS_postprocessing.h5"):
        warnings.warn("No PNS mass saved, run the supernova postprocessing first.")
        return
    pns_mass = load_dataset(
        simulation.storage_path,
        "PNS_postprocessing.h5",
        "global/mass",
        False,
    )

    start_index = len(postprocessed)
    end_index = min(len(file_list), len(pns_mass))
    file_list = file_list[start_index:end_index]
    index500 = np.argmax(simulation.grid.radius >= (500 * u.km))
    dOmega = simulation.cell.dOmega(simulation.ghost)
    theta = simulation.cell.theta(simulation.ghost)
    phi = simulation.cell.phi(simulation.ghost)
    index_theta = None
    index_phi = None
    weights = [[None, None], [None, None], [None, None]]
    indices_theta_avg = 0
    domega_obs = None
    if simulation.dim > 1:
        index_theta = [
            np.argmax(simulation.grid.theta >= (dd * u.deg))
            for dd in (np.arange(19) * 10)
        ]
        index_theta[-1] = -1 if index_theta[-1] == 0 else index_theta[-1]
        weights[0] = compute_weights(
            theta, phi, 0 * u.radian, np.pi / 2 * u.radian, simulation.dim
        )
        indices_theta_avg = [
            np.argmax(theta >= (6 * u.deg)),
            np.argmax(theta >= (86 * u.deg)),
            np.argmax(theta >= (94 * u.deg)),
        ]
    if simulation.dim > 2:
        index_phi = 0
        weights[1] = compute_weights(
            theta, phi, 0 * u.radian, 0 * u.radian, simulation.dim
        )
        domega_obs, alpha, beta = build_obs_grid(256)
        weights[2] = compute_weights(theta, phi, alpha, beta, simulation.dim)
    neutrino_dict = declare_neutrino_dictionary(simulation.dim, index_theta)
    neutrino_metadata = neutrino_quantity_metadata(index_theta)
    global_time = []
    reduced_time = []
    check_index = 0
    pbar = tqdm(file_list, "Neutrino postprocessing", colour="#07ff9c")
    for i, file_name in enumerate(pbar):
        pbar.set_postfix(file=file_name)
        time = simulation.time(file_name)
        global_time.append(time)
        save_loc = compute_neutrino_postprocessing(
            neutrino_dict,
            simulation,
            file_name,
            dOmega,
            weights,
            domega_obs,
            index500,
            index_theta,
            indices_theta_avg,
            index_phi,
            pns_mass[i + start_index],
        )
        if save_loc:
            reduced_time.append(time)
        postprocessed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            pbar.set_description("Saving.")
            save_merge_dictionary_hdf(
                simulation,
                (global_time, neutrino_dict, postprocessed),
                simulation.storage_path,
                "neutrino_postprocessing.h5",
                ["nue", "nua", "nux"],
                **neutrino_metadata,
            )
            check_index = 0
            neutrino_dict = declare_neutrino_dictionary(
                simulation.dim,
                index_theta,
            )
            global_time = []
            reduced_time = []
            pbar.set_description("Neutrino postprocessing")
        else:
            check_index += 1
    save_merge_dictionary_hdf(
        simulation,
        (global_time, neutrino_dict, postprocessed),
        simulation.storage_path,
        "neutrino_postprocessing.h5",
        ["nue", "nua", "nux"],
        **neutrino_metadata,
    )
    tqdm.write("Neutrino postprocessing completed!")


def supernova_gws_postprocessing(
    simulation: Simulation,
    save_checkpoints: bool = True,
) -> None:
    if simulation.dim == 1:
        return
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    file_list = simulation.hdf_file_list
    postprocessed = []
    if check_existence(simulation, "GWs_postprocessing.h5"):
        postprocessed = load_dataset(
            simulation.storage_path,
            "GWs_postprocessing.h5",
            "processed_hdf",
            False,
        )
        if len(postprocessed) == len(file_list) and postprocessed[-1] == file_list[-1]:
            return
        tqdm.write(
            "Checkpoint found for the Gravitational Waves Postprocessing, "
            f" starting from {postprocessed[-1]}."
        )
    else:
        tqdm.write("Starting Gravitational Waves Postprocessing.")
    gws_dict = declare_GW_dictionary(simulation.dim)
    gws_metadata = GWs_metadata()
    global_time = []
    check_index = 0
    dOmega = simulation.cell.dOmega(simulation.ghost)
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    gradient = None
    Ylm = None
    if simulation.dim == 3:
        theta = simulation.grid.theta
        phi = simulation.grid.phi
        radius = simulation.grid.radius
        gradient = []
        Ylm = []
        harm = SphericalHarmonics()
        for m in range(-2, 3):
            Y2m_r = (
                harm.Ylm_conj(m, 2, theta, phi)[..., None] * radius[None, None, :] ** 2
            )
            gradient.append(simulation.grid.gradient(Y2m_r, True))
            Ylm.append(Y2m_r * radius[None, None, ...] ** 2 * dOmega[..., None])

    core_radius = load_dataset(
        simulation.storage_path,
        "PNScore_postprocessing.h5",
        "local/r",
        False,
    )
    pns_radius = load_dataset(
        simulation.storage_path,
        "PNS_postprocessing.h5",
        "local/r",
        False,
    )
    innercore_radius = load_dataset(
        simulation.storage_path,
        "innercore_postprocessing.h5",
        "local/r",
        False,
    )
    shock_radius = load_dataset(
        simulation.storage_path,
        "shock_postprocessing.h5",
        "local/r",
        False,
    )
    start_index = len(postprocessed)
    end_index = min(
        len(file_list),
        core_radius.shape[-1],
        pns_radius.shape[-1],
        innercore_radius.shape[-1],
        shock_radius.shape[-1],
    )
    file_list = file_list[start_index:end_index]
    pbar = tqdm(file_list, "GWs postprocessing", colour="#ffa407")
    for i, file_name in enumerate(pbar):
        pbar.set_postfix(file=file_name)
        time = simulation.time(file_name)
        global_time.append(time)
        compute_gws_postprocessing(
            gws_dict,
            simulation,
            file_name,
            dOmega,
            dV,
            core_radius[..., i + start_index],
            pns_radius[..., i + start_index],
            innercore_radius[..., i + start_index],
            shock_radius[..., i + start_index],
            gradient,
            Ylm,
        )
        postprocessed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            pbar.set_description("Saving.")
            save_merge_dictionary_hdf(
                simulation,
                (global_time, gws_dict, postprocessed),
                simulation.storage_path,
                "GWs_postprocessing.h5",
                None,
                **gws_metadata,
            )
            check_index = 0
            gws_dict = declare_GW_dictionary(simulation.dim)
            global_time = []
            pbar.set_description("GWs postprocessing")
        else:
            check_index += 1

    save_merge_dictionary_hdf(
        simulation,
        (global_time, gws_dict, postprocessed),
        simulation.storage_path,
        "GWs_postprocessing.h5",
        None,
        **gws_metadata,
    )


def tidal_deformability_postprocesing(
    simulation: Simulation,
    save_checkpoints: bool = True,
) -> None:
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    file_list = simulation.hdf_file_list
    postprocessed = []
    if check_existence(simulation, "tidal_postprocessing.h5"):
        postprocessed = load_dataset(
            simulation.storage_path,
            "tidal_postprocessing.h5",
            "processed_hdf",
            False,
        )
        if (len(postprocessed) == len(file_list)) and (
            postprocessed[-1] == file_list[-1]
        ):
            return
        tqdm.write(
            "Checkpoint found for the Tidal Deformabilty Postprocessing, "
            f" starting from {postprocessed[-1]}."
        )
    else:
        tqdm.write("Starting Tidal Deformabilty Postprocessing.")

    if not check_existence(simulation, "profiles_postprocessing.h5"):
        warnings.warn("No profiles saved, run the supernova postprocessing first.")
        return
    rho = load_dataset(
        simulation.storage_path,
        "profiles_postprocessing.h5",
        "rho_profile",
        True,
    )
    tm = rho.time
    rho = rho.data
    cs = load_dataset(
        simulation.storage_path,
        "profiles_postprocessing.h5",
        "soundspeed_profile",
        False,
    )
    s = load_dataset(
        simulation.storage_path,
        "profiles_postprocessing.h5",
        "entropy_profile",
        False,
    )
    p = load_dataset(
        simulation.storage_path,
        "profiles_postprocessing.h5",
        "gas_pressure_profile",
        False,
    )
    start_index = len(postprocessed)
    end_index = min(len(file_list), len(tm))
    file_list = file_list[start_index:end_index]
    pbar = tqdm(
        file_list,
        "Tidal deformability postprocessing",
        colour="#c507ff",
    )
    global_time = []
    tidal_dict = declare_tidal_dictionary()
    tidal_meta = tidal_quantity_metadata()
    check_index = 0
    for i, file_name in enumerate(pbar):
        pbar.set_postfix(file=file_name)
        time = tm[i + start_index]
        global_time.append(time)
        compute_tidal_love(
            tidal_dict,
            simulation,
            rho[:, i + start_index],
            p[:, i + start_index],
            cs[:, i + start_index],
            s[:, i + start_index],
        )
        postprocessed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            pbar.set_description("Saving.")
            save_merge_dictionary_hdf(
                simulation,
                (global_time, tidal_dict, postprocessed),
                simulation.storage_path,
                "tidal_postprocessing.h5",
                ["PNS", "core"],
                **tidal_meta,
            )
            check_index = 0
            tidal_dict = declare_tidal_dictionary()
            global_time = []
            pbar.set_description("Tidal deformability postprocessing")
        else:
            check_index += 1

    save_merge_dictionary_hdf(
        simulation,
        (global_time, tidal_dict, postprocessed),
        simulation.storage_path,
        "tidal_postprocessing.h5",
        ["PNS", "core"],
        **tidal_meta,
    )


def spherical_harmonics_postprocessing(
    simulation: Simulation,
    quantity_name: str,
    save_checkpoints: bool = True,
) -> None:
    comp = None
    rad = None
    isradius = False
    if simulation.dim == 1:
        return
    out_name = f"{quantity_name}.h5"
    if "-" in quantity_name:
        isradius = True
        splits = quantity_name.split("-")
        if len(splits) > 3:
            raise ValueError(f"Too many values, expected at most 3, got {len(splits)}.")
        quantity_name = splits[0]
        splits = splits[1:]
        if len(splits) == 1:
            splits = splits[0]
            if splits in ["avg", "min", "max", "full"]:
                rad = "full"
            else:
                comp = splits
        elif len(splits) == 2:
            if splits[0] in ["avg", "min", "max", "full"]:
                rad = "full"
                comp = splits[1]
            elif splits[1] in ["avg", "min", "max", "full"]:
                rad = "full"
                comp = splits[0]
            else:
                raise TypeError("Quantity not recognised.")
        if comp is not None:
            out_name = f"{quantity_name}_{comp}.h5"
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    file_list = simulation.hdf_file_list
    postprocessed = []
    if check_existence(simulation, out_name):
        postprocessed = load_dataset(
            simulation.storage_path,
            out_name,
            "processed_hdf",
            False,
        )
        if (len(postprocessed) == len(file_list)) and (
            postprocessed[-1] == file_list[-1]
        ):
            return
        tqdm.write(
            f"Checkpoint found for the {quantity_name} spherical harmonis decomposition, "
            f" starting from {postprocessed[-1]}."
        )
    else:
        tqdm.write(f"Starting spherical harmonics decomposition of {quantity_name}.")
    assert hasattr(simulation, quantity_name), (
        f"The selected quantity ({quantity_name}) is not present in the simulation."
    )
    deco_quantity = None
    spher_meta = None
    if isradius:
        try:
            if comp is not None and rad is not None:
                deco_quantity = getattr(simulation, quantity_name)(comp=comp, rad=rad)
            elif rad is not None and comp is None:
                deco_quantity = getattr(simulation, quantity_name)(rad=rad)
            else:
                deco_quantity = getattr(simulation, quantity_name)(
                    comp=comp, rad="full"
                )
            tm_deco = deco_quantity.time
            deco_quantity = deco_quantity.data
            spher_meta = create_metadata_spherical_harmonics(deco_quantity)
        except IndexError:
            isradius = False
            deco_quantity = None

    start_index = len(postprocessed)
    end_index = (
        min(len(file_list), deco_quantity.shape[-1])
        if deco_quantity is not None
        else len(file_list)
    )
    file_list = file_list[start_index:end_index]
    pbar = tqdm(
        file_list,
        f"Spherical harmonics decomposition of {quantity_name}",
        colour="#ff0739",
    )
    global_time = []
    sph_harm_dict = declare_spherical_harmonics_dicionary(simulation.dim)
    check_index = 0
    dOmega = simulation.cell.dOmega(simulation.ghost)
    spher_harm = SphericalHarmonics()
    lmax = int((simulation.theta / 2 - 1) * 1.11)

    for i, file_name in enumerate(pbar):
        pbar.set_postfix(file=file_name)
        if isradius:
            time = tm_deco[i + start_index]
            dd_qt = deco_quantity[i + start_index]
        else:
            time = simulation.time(file_name)
            if comp is not None:
                dd_qt = getattr(simulation, quantity_name)(file_name, comp=comp)
            else:
                dd_qt = getattr(simulation, quantity_name)(file_name)
        global_time.append(time)
        compute_spherical_harmonics_decomposition(
            sph_harm_dict,
            spher_harm,
            dd_qt,
            simulation,
            dOmega,
            lmax,
        )
        postprocessed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            pbar.set_description("Saving.")
            if spher_meta is None:
                spher_meta = create_metadata_spherical_harmonics(dd_qt)
            save_merge_dictionary_hdf(
                simulation,
                (global_time, spher_harm, postprocessed),
                simulation.storage_path,
                out_name,
                None,
                **spher_meta,
            )
            check_index = 0
            spher_harm = declare_spherical_harmonics_dicionary()
            global_time = []
            pbar.set_description(
                f"Spherical harmonics decomposition of {quantity_name}",
            )
        else:
            check_index += 1
    if spher_meta is None:
        spher_meta = create_metadata_spherical_harmonics(dd_qt)
    save_merge_dictionary_hdf(
        simulation,
        (global_time, spher_harm, postprocessed),
        simulation.storage_path,
        out_name,
        None,
        **spher_meta,
    )
