from __future__ import annotations


def declare_GW_dictionary(dim: int) -> dict:
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
