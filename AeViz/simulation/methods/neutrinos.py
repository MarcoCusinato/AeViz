from AeViz.simulation.methods import *

"""
Functions to handle neutrino data from a simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## NEUTRINO DATA
## -----------------------------------------------------------------

## ENERGY DEPENDENT

@smooth
@hdf_isopen
def neutrino_energy_density(self, file_name, **kwargs):
    nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['neutrino/e'][..., 0]), self.dim)
    nu_ene[..., 2] /= 4
    return (aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens',
                  r'$E_{\nu_\mathrm{e}}$', 'viridis', [1e28, 1e32], True),
            aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens',
                  r'$E_{\overline{\nu}_\mathrm{e}}$', 'magma',
                  [1e28, 1e32], True),
            aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens',
                  r'$E_{\nu_\mathrm{x}}$', 'plasma', [1e28, 1e32], True)
            )

@smooth
@hdf_isopen
def neutrino_momenta(self, file_name, **kwargs):
    """
    In the comoving rest frame of the fluid are equal to the
    neutrino energy fluxes
    """
    nu_flux = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['neutrino/e'][..., 1:]), self.dim)
    ## Insert a new axis to be consistent with the other dimensions
    if self.dim == 1:
        nu_flux[..., 2] /= 4
        return (aerray(nu_flux[..., 0], u.erg / u.s / u.cm ** 2, 'nue_fdens',
                  r'$F_{\nu_\mathrm{e}}$', 'PiYG_r', [-1e40, 1e40], True),
                aerray(nu_flux[..., 1], u.erg / u.s / u.cm ** 2, 'nua_fdens',
                  r'$F_{\overline{\nu}_\mathrm{e}}$', 'PuOr_r',
                  [-1e40, 1e40], True),
                aerray(nu_flux[..., 2], u.erg / u.s / u.cm ** 2, 'nux_fdens',
                  r'$F_{\nu_\mathrm{x}}$', 'RdYlBu_r', [-1e40, 1e40], True)
                )
    nu_flux[..., 2, :] /= 4
    return (aerray(nu_flux[..., 0, :], u.erg / u.s / u.cm ** 2, 'nue_fdens',
                  [r'$F_{\nu_\mathrm{e}, r}$', r'$F_{\nu_\mathrm{e}, \theta}$',
                   r'$F_{\nu_\mathrm{e}, \phi}$'],
                  ['PiYG_r', 'PRGn_r', 'BrBG_r'],
                  [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_flux[..., 1, :], u.erg / u.s / u.cm ** 2, 'nua_fdens',
                [r'$F_{\overline{\nu}_\mathrm{e}, r}$', r'$F_{\overline{\nu}_\mathrm{e}, \theta}$',
                 r'$F_{\overline{\nu}_\mathrm{e}, \phi}$'],
                ['PuOr_r', 'RdGy_r', 'RdBu_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_flux[..., 2, :], u.erg / u.s / u.cm ** 2, 'nux_fdens',
                [r'$F_{\nu_\mathrm{x}, r}$', r'$F_{\nu_\mathrm{x}, \theta}$',
                 r'$F_{\nu_\mathrm{x}, \phi}$'],
                ['RdYlBu_r', 'RdYlGn_r', 'Spectral_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True)
            )

@smooth
@hdf_isopen
def neutrino_momenta_opacities(self, file_name, **kwargs):
    nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['neutrino/oe'][..., 1:]), self.dim)
    if self.dim == 1:
        nu_opac = nu_opac[..., None]
    
    if self.dim == 1:
        return (aerray(nu_opac[..., 0], u.erg / u.s / u.cm ** 3, 'nue_kappa',
                  r'$\kappa_{\nu_\mathrm{e}}$', 'PiYG_r', [-1e40, 1e40], True),
                aerray(nu_opac[..., 1], u.erg / u.s / u.cm ** 3, 'nua_fdens',
                  r'$\kappa_{\overline{\nu}_\mathrm{e}}$', 'PuOr_r',
                  [-1e40, 1e40], True),
                aerray(nu_opac[..., 2], u.erg / u.s / u.cm ** 3, 'nux_kappa',
                  r'$\kappa_{\nu_\mathrm{x}}$', 'RdYlBu_r', [-1e40, 1e40], True)
                )
    return (aerray(nu_opac[..., 0, :], u.erg / u.s / u.cm ** 3, 'nue_kappa',
                  [r'$\kappa_{\nu_\mathrm{e}, r}$', r'$\kappa_{\nu_\mathrm{e}, \theta}$',
                   r'$\kappa_{\nu_\mathrm{e}, \phi}$'],
                  ['PiYG_r', 'PRGn_r', 'BrBG_r'],
                  [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_opac[..., 1, :], u.erg / u.s / u.cm ** 3, 'nua_kappa',
                [r'$\kappa_{\overline{\nu}_\mathrm{e}, r}$',
                 r'$\kappa_{\overline{\nu}_\mathrm{e}, \theta}$',
                 r'$\kappa_{\overline{\nu}_\mathrm{e}, \phi}$'],
                ['PuOr_r', 'RdGy_r', 'RdBu_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_opac[..., 2, :], u.erg / u.s / u.cm ** 3, 'nux_kappa',
                [r'$\kappa_{\nu_\mathrm{x}, r}$', r'$\kappa_{\nu_\mathrm{x}, \theta}$',
                 r'$\kappa_{\nu_\mathrm{x}, \phi}$'],
                ['RdYlBu_r', 'RdYlGn_r', 'Spectral_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True)
            )

@smooth
def neutrino_number_density(self, file_name, **kwargs):
    edens = list(self.neutrino_energy_density(file_name, **kwargs))
    de = self.cell.E_nu().to(u.erg)
    edens = [ed / de for ed in edens]
    [ed.set(label=s1, name=s2, limits=s3, cmap=s4, log=True) for (s1, s2, s3, s4) 
             in zip([r'$N_{\nu_\mathrm{e}}$',
                     r'$N_{\overline{\nu}_\mathrm{e}}$',
                     r'$N_{\nu_\mathrm{x}}$'],
                    ['Nnue', 'Nnua', 'Nnux'],
                    [[1e33, 1e36], [1e31, 1e34], [1e32, 1e35]],
                    ['gnuplot', 'gnuplot_2', 'CMRmap']) for ed in edens]
    return tuple(edens)

@smooth
def neutrino_mean_energy(self, file_name, **kwargs):
    """
    Average neutrino energy per cell so 
    <e> = sum_w E_nu(w) / sum_w N_nu(w),
        with w the center of the neutrino bin
    """
    edens = list(self.neutrino_energy_density(file_name, **kwargs))
    num_den = list(self.neutrino_number_density(file_name, **kwargs))
    mean_ene = [(ed.sum(axis=-1)/nd.sum(axis=-1)).to(u.MeV) for 
                (ed, nd) in zip(edens, num_den)]
    [me.set(label=s1, name=s2, limits=s3, cmap=s4) for (s1, s2, s3, s4)
     in zip([r'$\langle E_{\nu_\mathrm{e}}\rangle$',
                r'$\langle E_{\overline{\nu}_\mathrm{e}}\rangle$',
                r'$\langle E_{\nu_\mathrm{x}}\rangle$'],
            ['Enue', 'Enua', 'Enux'],
            [[0, 120], [0, 90], [0, 100]],
            ['ocean', 'gist_earth', 'terrain']) for me in mean_ene]
    return tuple(mean_ene)

## GREY

@smooth
@hdf_isopen
def neutrino_energy_density_grey(self, file_name, **kwargs):
    nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['/neutrinogrey/egrey'][..., 0]), self.dim)
    nu_ene[..., 2] /= 4
    return (aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens',
                  r'$E_{\nu_\mathrm{e}}$', 'viridis', [1e28, 1e32], True),
            aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens',
                  r'$E_{\overline{\nu}_\mathrm{e}}$', 'magma',
                  [1e28, 1e32], True),
            aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens',
                  r'$E_{\nu_\mathrm{x}}$', 'plasma', [1e28, 1e32], True)
            )

@smooth
@hdf_isopen
def neutrino_momenta_grey(self, file_name, **kwargs):
    """
    In the comoving rest frame of the fluid are equal to the 
    neutrino energy fluxes
    """
    nu_flux =  self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['/neutrinogrey/egrey'][..., 1:]), self.dim)
    if self.dim == 1:
        nu_flux[..., 2] /= 4
        return (aerray(nu_flux[..., 0], u.erg / u.s / u.cm ** 2, 'nue_fdens',
                  r'$F_{\nu_\mathrm{e}}$', 'PiYG_r', [-1e40, 1e40], True),
                aerray(nu_flux[..., 1], u.erg / u.s / u.cm ** 2, 'nua_fdens',
                  r'$F_{\overline{\nu}_\mathrm{e}}$', 'PuOr_r',
                  [-1e40, 1e40], True),
                aerray(nu_flux[..., 2], u.erg / u.s / u.cm ** 2, 'nux_fdens',
                  r'$F_{\nu_\mathrm{x}}$', 'RdYlBu_r', [-1e40, 1e40], True)
                )
    nu_flux[..., 2, :] /= 4
    return (aerray(nu_flux[..., 0, :], u.erg / u.s / u.cm ** 2, 'nue_fdens',
                  [r'$F_{\nu_\mathrm{e}, r}$', r'$F_{\nu_\mathrm{e}, \theta}$',
                   r'$F_{\nu_\mathrm{e}, \phi}$'],
                  ['PiYG_r', 'PRGn_r', 'BrBG_r'],
                  [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_flux[..., 1, :], u.erg / u.s / u.cm ** 2, 'nua_fdens',
                [r'$F_{\overline{\nu}_\mathrm{e}, r}$', r'$F_{\overline{\nu}_\mathrm{e}, \theta}$',
                 r'$F_{\overline{\nu}_\mathrm{e}, \phi}$'],
                ['PuOr_r', 'RdGy_r', 'RdBu_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_flux[..., 2, :], u.erg / u.s / u.cm ** 2, 'nux_fdens',
                [r'$F_{\nu_\mathrm{x}, r}$', r'$F_{\nu_\mathrm{x}, \theta}$',
                 r'$F_{\nu_\mathrm{x}, \phi}$'],
                ['RdYlBu_r', 'RdYlGn_r', 'Spectral_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True)
            )