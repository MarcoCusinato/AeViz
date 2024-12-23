from AeViz.simulation.methods import *
from AeViz.utils.spherical_harmonics_radial import (calculate_rho_decomposition,
                                                    get_sph_profile,
                                                    get_sph_profiles_r,
                                                    get_data_for_barcode)
## -----------------------------------------------------------------
## SPHERICAL HARMONICS
## -----------------------------------------------------------------

@smooth
def rho_spherical_harmonics(self, l=0, m=None, zero_norm=True,
                            rhomin=None, rhomax=None, r=None,
                            save_checkpoints=True, **kwargs):
    if m is None:
        calculate_rho_decomposition(self, save_checkpoints, msum=True)
    else:
        calculate_rho_decomposition(self, save_checkpoints)
    if [rhomin, rhomax, r] == [None, None, None]:
        radius = self.cell.radius(self.ghost)
        time, rlm = get_sph_profile(self, l, m)
        if zero_norm:
            if m is None:
                _, r00 = get_sph_profile(self, 0)
            else:
                _, r00 = get_sph_profile(self, 0, 0)
            rlm /= r00
        return time, radius, rlm
    else:
        return get_sph_profiles_r(self, l=l, m=m, rhomin=rhomin,
                                    rhomax=rhomax, r=r, zero_norm=zero_norm)

def data_for_barcode(self, lmin=None, lmax=None, msum=False,
                        r=None, rhomin=None, rhomax=None, save_checkpoints=True,
                        zero_norm=True, **kwargs):
    """
    Returns the data for the barcode plot. The data is the spherical
    harmonics decomposition of the density.
    """
    calculate_rho_decomposition(self, save_checkpoints, msum=True)
    return get_data_for_barcode(self, lmax=lmax, lmin=lmin, msum=msum,
                                r=r, rhomin=rhomin, rhomax=rhomax,
                                zero_norm=zero_norm)