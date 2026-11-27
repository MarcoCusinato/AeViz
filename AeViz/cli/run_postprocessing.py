import argparse
import os
import warnings

import numpy as np
from tqdm import tqdm

from AeViz.simulation.simulation import Simulation
from AeViz.utils.physics.postprocessing import (
    spherical_harmonics_postprocessing,
    supernova_gws_postprocessing,
    supernova_neutrino_postprocessing,
    supernova_postprocessing,
    tidal_deformability_postprocesing,
)
from AeViz.utils.utils import restart_from


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sim-name",
        type=str,
        required=True,
        help="Name of the simulation of which postprocessing ' \
                            'will be done.",
    )
    parser.add_argument(
        "--sim-path",
        type=str,
        default=None,
        required=False,
        help="Path of the simulation.",
    )
    parser.add_argument(
        "--tob",
        default=None,
        type=float,
        required=False,
        help="Sets the time of bounce manually and stores it into \
                            a file",
    )
    parser.add_argument(
        "--spherical-harmonics",
        action="store_true",
        default=False,
        help="Enables the computation of the spherical harmonics",
    )
    parser.add_argument(
        "--show-warnings",
        action="store_true",
        help="Do not suppress numerical warnings",
        default=False,
    )
    parser.add_argument(
        "--restart-from",
        default=None,
        type=str,
        help="The name of the timestep from which to restart the postprocessing.",
    )
    args = parser.parse_args()
    sim = Simulation(args.sim_name, args.sim_path)
    sim.neglect_new()  ## Probably not needed, here just for legacy reasons.
    if args.tob:
        sim.tob = args.tob
        np.savetxt(
            os.path.join(sim.storage_path, "tob.dat"),
            np.array([args.tob]),
        )
    if args.restart_from is not None:
        restart_from(
            sim,
            args.restart_from,
            [
                "PNS_postprocessing.h5",
                "gain_postprocessing.h5",
                "innercore_postprocessing.h5",
                "profiles_postprocessing.h5",
                "PNScore_postprocessing.h5",
                "shock_postprocessing.h5",
            ],
        )
    postprocessing = [
        supernova_postprocessing,
        supernova_neutrino_postprocessing,
        supernova_gws_postprocessing,
        tidal_deformability_postprocesing,
    ]

    pbar = tqdm(
        postprocessing,
        f"Postprocessing of {sim.simulation_name}",
        colour="#aaff00",
    )
    sph_harms = ["rho"]
    if args.show_warnings:
        for pp in pbar:
            pp(sim, True)
        pbar_harmonics = tqdm(
            sph_harms,
            f"Spherical harmonics for {sim.simulation_name}",
            colour="#e8f405",
        )
        for sp in pbar_harmonics:
            spherical_harmonics_postprocessing(sim, sp, True)
    else:
        with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            warnings.filterwarnings("ignore", message="Mean of empty slice")
            for pp in pbar:
                pp(sim, True)
            pbar_harmonics = tqdm(
                sph_harms,
                f"Spherical harmonics for {sim.simulation_name}",
                colour="#e8f405",
            )
            for sp in pbar_harmonics:
                spherical_harmonics_postprocessing(sim, sp, True)


if __name__ == "__main__":
    main()
