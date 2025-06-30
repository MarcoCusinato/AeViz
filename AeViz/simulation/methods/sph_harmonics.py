from AeViz.simulation.methods import *
from AeViz.utils.physics.spherical_harmonics_radial import (calculate_rho_decomposition,
                                                    get_sph_profile,
                                                    get_sph_profiles_r,
                                                    get_data_for_barcode)
"""
Functions to handle spherical harmonics decomposition of the density
data from a simulation in spherical coordinates.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## SPHERICAL HARMONICS
## -----------------------------------------------------------------

@smooth
def rho_spherical_harmonics(self, l=0, m=None, zero_norm=True,
                            rhomin=None, rhomax=None, r=None,
                            save_checkpoints=True, mode='radius', **kwargs):
    if m is None:
        calculate_rho_decomposition(self, save_checkpoints, msum=True)
        if zero_norm:
            lb = r'$\tilde{\rho}_{%d}/\tilde{\rho}_{00}$' % l
        else:
            lb = r'$\tilde{\rho}_{%d}$' % l
    else:
        calculate_rho_decomposition(self, save_checkpoints)
        if zero_norm:
            lb = r'$\tilde{\rho}_{%d,%d}/\tilde{\rho}_{00}$' % (l, m)
        else:
            lb = r'$\tilde{\rho}_{%d,%d}$' % (l, m)
    if [rhomin, rhomax, r] == [None, None, None]:
        radius = self.cell.radius(self.ghost)
        time, rlm = get_sph_profile(self, l, m)
        if zero_norm:
            if m is None:
                _, r00 = get_sph_profile(self, 0)
            else:
                _, r00 = get_sph_profile(self, 0, 0)
            rlm /= r00
        time = aerray(time, u.s, name='time', label=r'$t-t_\mathrm{b}$',
                      cmap=None, log=False,
                 limits=[-0.005, time.max()])
        rlm = aerray(rlm, u.dimensionless_unscaled, 'rho_sph_harmonics',
                     lb, 'YlGn',
                     [rlm.min() * 1.1, rlm.max() * 0.9], False)
        return aeseries(rlm, time=time, radius=radius)
    else:
        time, rlm = get_sph_profiles_r(self, l=l, m=m, rhomin=rhomin,
                                    rhomax=rhomax, r=r, zero_norm=zero_norm,
                                    mode=mode)
        time = aerray(time, u.s, name='time', label=r'$t-t_\mathrm{b}$',
                      cmap=None, log=False,
                 limits=[-0.005, time.max()])
        rlm = aerray(rlm, u.dimensionless_unscaled, 'rho_sph_harmonics',
                     lb, 'YlGn',
                     [rlm.min() * 1.1, rlm.max() * 0.9], False)
        return aeseries(rlm, time=time)

def data_for_barcode(self, lmin=None, lmax=None, msum=False,
                        r=None, rhomin=None, rhomax=None, save_checkpoints=True,
                        zero_norm=True, mode:Literal['radius', 'mass']='radius',
                        **kwargs):
    """
    Returns the data for the barcode plot. The data is the spherical
    harmonics decomposition of the density.
    """
    calculate_rho_decomposition(self, save_checkpoints, msum=True)
    time, Y, data = get_data_for_barcode(self, lmax=lmax, lmin=lmin, msum=msum,
                                r=r, rhomin=rhomin, rhomax=rhomax,
                                zero_norm=zero_norm, mode=mode)
    time = aerray(time, u.s, name='time', label=r'$t-t_\mathrm{b}$',
                  cmap=None, log=False,
             limits=[-0.005, time.max()])
    if msum:
        Y = aerray(Y, u.dimensionless_unscaled, name='l', label=r'$l$',
                   cmap=None, log=False, limits=[Y.min(), Y.max()])
    else:
        Y = aerray(Y, u.dimensionless_unscaled, name='lm', label=r'$2l+1$',
                   cmap=None, log=False, limits=[Y.min(), Y.max()])
    if zero_norm:
        data = aerray(data, u.dimensionless_unscaled, name='rho_sph_harmonics',
                  label=r'$\tilde{\rho}/\tilde{\rho}_{00}$', cmap='seismic',
                  log=False, limits=[data.min() * 1.1, data.max() * 0.9])
    else:
        data = aerray(data, u.dimensionless_unscaled, name='rho_sph_harmonics',
                  label=r'$\tilde{\rho}$', cmap='seismic',
                  log=False, limits=[data.min() * 1.1, data.max() * 0.9])
    return aeseries(data, time=time, Y=Y)