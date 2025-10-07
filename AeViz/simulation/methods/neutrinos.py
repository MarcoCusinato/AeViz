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

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_energy_density(self, file_name, **kwargs):
    """
    Neutrino energy density.
    Now with NOTRINO case!
    """
    if not kwargs['notrino']:
        nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['neutrino/e'][..., 0]), self.dim)
    else:
        nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['notrino/notrino_e'][...]), self.dim)
        # Adding a fictitious to account for missing energy bin and thus
        # (hopefully) not having to change a lot in the library
        nu_ene = np.expand_dims(nu_ene, axis = -2)
    nu_ene[..., 2] /= 4
    return (aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens',
                  r'$E_{\nu_\mathrm{e}}$', 'viridis', [1e28, 1e32], True),
            aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens',
                  r'$E_{\overline{\nu}_\mathrm{e}}$', 'magma',
                  [1e28, 1e32], True),
            aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens',
                  r'$E_{\nu_\mathrm{x}}$', 'plasma', [1e28, 1e32], True)
            )

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_momenta(self, file_name, **kwargs):
    """
    In the comoving rest frame of the fluid are equal to the
    neutrino energy fluxes.
    Now with NOTRINO case!
    """
    if not kwargs['notrino']:
        nu_flux = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['neutrino/e'][..., 1:]), self.dim)
    else:
        nu_flux = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['notrino/notrino_f'][...]), self.dim)
        # Adding a fictitious to account for missing energy bin and thus
        # (hopefully) not having to change a lot in the library
        nu_flux = np.expand_dims(nu_flux, axis = -2)
    ## Insert a new axis to be consistent with the other dimensions
    if self.dim == 1 or kwargs['notrino']:
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

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_energy_opacity(self, file_name, **kwargs):
    """
    Neutrino opacity ("temporal" component of the four-vector).
    This is the absorption/emission opacity.
    Now with NOTRINO case!
    """
    if not kwargs['notrino']:
        nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['neutrino/oe'][..., 0]), self.dim)
    else:
        nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['notrino/notrino_kae'][...]), self.dim)
        # Adding a fictitious to account for missing energy bin and thus
        # (hopefully) not having to change a lot in the library
        nu_opac = np.expand_dims(nu_opac, axis = -2)
    return (aerray(nu_opac[..., 0], u.erg / u.s / u.cm ** 3, 'nue_kappa_ae',
                r'$\kappa_{\mathrm{ae}, \nu_\mathrm{e}}$', 'PiYG_r', [-1e40, 1e40], True),
            aerray(nu_opac[..., 1], u.erg / u.s / u.cm ** 3, 'nua_kappa_ae',
                r'$\kappa_{\mathrm{ae}, \overline{\nu}_\mathrm{e}}$', 'PuOr_r',
                [-1e40, 1e40], True),
            aerray(nu_opac[..., 2], u.erg / u.s / u.cm ** 3, 'nux_kappa_ae',
                r'$\kappa_{\mathrm{ae}, \nu_\mathrm{x}}$', 'RdYlBu_r', [-1e40, 1e40], True)
            )
    

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_momenta_opacities(self, file_name, **kwargs):
    """
    Neutrino opacities ("spatial" components of the four-vector).
    This is the transport opacity.
    Now with NOTRINO case!
    """
    if not kwargs['notrino']:
        nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['neutrino/oe'][..., 1:]), self.dim)
    else:
        notrino = True
        nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['notrino/notrino_ktr'][...]), self.dim)
        # Adding a fictitious to account for missing energy bin and thus
        # (hopefully) not having to change a lot in the library
        nu_opac = np.expand_dims(nu_opac, axis = -2)
    #if self.dim == 1:
    #    nu_opac = nu_opac[..., None]
    
    if self.dim == 1 or kwargs['notrino']:
        return (aerray(nu_opac[..., 0], u.erg / u.s / u.cm ** 3, 'nue_kappa_tr',
                  r'$\kappa_{\mathrm{tr}, \nu_\mathrm{e}}$', 'PiYG_r', [-1e40, 1e40], True),
                aerray(nu_opac[..., 1], u.erg / u.s / u.cm ** 3, 'nua_kappa_tr',
                  r'$\kappa_{\mathrm{tr}, \overline{\nu}_\mathrm{e}}$', 'PuOr_r',
                  [-1e40, 1e40], True),
                aerray(nu_opac[..., 2], u.erg / u.s / u.cm ** 3, 'nux_kappa_tr',
                  r'$\kappa_{\mathrm{tr}, \nu_\mathrm{x}}$', 'RdYlBu_r', [-1e40, 1e40], True)
                )
    return (aerray(nu_opac[..., 0, :], u.erg / u.s / u.cm ** 3, 'nue_kappa_tr',
                  [r'$\kappa^r_{\mathrm{tr}, \nu_\mathrm{e}}$', \
                   r'$\kappa^{\theta}_{\mathrm{tr}, \nu_\mathrm{e}}$',
                   r'$\kappa^{\phi}_{\mathrm{tr}, \nu_\mathrm{e}}$'],
                  ['PiYG_r', 'PRGn_r', 'BrBG_r'],
                  [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_opac[..., 1, :], u.erg / u.s / u.cm ** 3, 'nua_kappa_tr',
                [r'$\kappa^r_{\mathrm{tr}, \overline{\nu}_\mathrm{e}}$',
                 r'$\kappa^{\theta}_{\mathrm{tr}, \overline{\nu}_\mathrm{e}}$',
                 r'$\kappa^{\phi}_{\mathrm{tr}, \overline{\nu}_\mathrm{e}}$'],
                ['PuOr_r', 'RdGy_r', 'RdBu_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True),
            aerray(nu_opac[..., 2, :], u.erg / u.s / u.cm ** 3, 'nux_kappa_tr',
                [r'$\kappa^r_{\mathrm{tr}, \nu_\mathrm{x}}$', \
                 r'$\kappa^{\theta}_{\mathrm{tr}, \nu_\mathrm{x}}$',
                 r'$\kappa^{\phi}_{\mathrm{tr}, \nu_\mathrm{x}}$'],
                ['RdYlBu_r', 'RdYlGn_r', 'Spectral_r'],
                [[-1e40, 1e40], [-1e39, 1e39], [-1e39, 1e39]], True)
            )

@get_grid
@smooth
@hdf_isopen
def neutrino_absorption_opacity(self, file_name, **kwargs):
    """
    Neutrino absorption/emission opacity.
    """
    return neutrino_energy_opacity(self, file_name, **kwargs)

@get_grid
@smooth
@hdf_isopen
def neutrino_transport_opacities(self, file_name, **kwargs):
    """
    Neutrino transport opacities.
    """
    return neutrino_momenta_opacities(self, file_name, **kwargs)


@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_scattering_opacities(self, file_name, **kwargs):
    """
    Neutrino scattering opacities. 
    Remember the relation
        ktr = ksc + kae.
    """

    ktr = neutrino_transport_opacities(self, file_name, **kwargs)
    kae = neutrino_absorption_opacity(self, file_name, **kwargs)
    
    if self.dim == 1 or kwargs['notrino']:
        ksc_nue = ktr[0] - kae[0]
        ksc_nue.name = 'nue_kappa_sc'
        ksc_nue.label = r'$\kappa_{\mathrm{sc}, \nu_\mathrm{e}}$'
        ksc_nua = ktr[1] - kae[1]
        ksc_nua.name = 'nua_kappa_sc'
        ksc_nua.label = r'$\kappa^r_{\mathrm{sc}, \overline{\nu}_\mathrm{e}}$'
        ksc_nux = ktr[2] - kae[2]
        ksc_nux.name = 'nux_kappa_sc'
        ksc_nux.label = r'$\kappa_{\mathrm{sc}, \nu_\mathrm{x}}$'
    else:
        ksc_nue = ktr[0] - kae[0][..., None]
        ksc_nue.name = 'nue_kappa_sc'
        ksc_nue.label = r'$\kappa_{\mathrm{sc}, \nu_\mathrm{e}}$'
        ksc_nua = ktr[1] - kae[1][..., None]
        ksc_nua.name = 'nua_kappa_sc'
        ksc_nua.label = r'$\kappa^r_{\mathrm{sc}, \overline{\nu}_\mathrm{e}}$'
        ksc_nux = ktr[2] - kae[2][..., None]
        ksc_nux.name = 'nux_kappa_sc'
        ksc_nux.label = r'$\kappa_{\mathrm{sc}, \nu_\mathrm{x}}$'
    return (ksc_nue, ksc_nua, ksc_nux)

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_number_density(self, file_name, **kwargs):
    # If notrino is used, the number density is returned in  the output,
    # and there is no need to compute it.
    if kwargs['notrino']:
        ndens = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['notrino/notrino_n'][...]), 
            self.dim)
        ndens = np.expand_dims(ndens, axis = -2)
        return (aerray(ndens[..., 0], u.cm ** (-3), 'Nnue',
                  r'$N_{\nu_\mathrm{e}}$', 'gnuplot', [1e33, 1e36], True),
                aerray(ndens[..., 1], u.cm ** (-3), 'Nnua',
                  r'$N_{\overline{\nu}_\mathrm{e}}$', 'gnuplot_2',
                  [1e31, 1e34], True),
                aerray(ndens[..., 2], u.cm ** (-3), 'Nnux',
                  r'$N_{\nu_\mathrm{x}}$', 'CMRmap', [1e32, 1e35], True)
                )
    else:
        edens = list(self.neutrino_energy_density(file_name))
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

@get_grid
@smooth
def neutrino_mean_energy(self, file_name,
                         comp:Literal['all', 'nue', 'nua', 'nux']='all',
                         **kwargs):
    """
    Average neutrino energy per cell so 
    <e> = sum_w E_nu(w) / sum_w N_nu(w),
        with w the center of the neutrino bin
    """
    edens = list(self.neutrino_energy_density(file_name))
    num_den = list(self.neutrino_number_density(file_name))
    mean_ene = [(ed.sum(axis=-1)/nd.sum(axis=-1)).to(u.MeV) for 
                (ed, nd) in zip(edens, num_den)]
    [me.set(label=s1, name=s2, limits=s3, cmap=s4) for (s1, s2, s3, s4)
     in zip([r'$\langle E_{\nu_\mathrm{e}}\rangle$',
                r'$\langle E_{\overline{\nu}_\mathrm{e}}\rangle$',
                r'$\langle E_{\nu_\mathrm{x}}\rangle$'],
            ['Enue', 'Enua', 'Enux'],
            [[0, 120], [0, 90], [0, 100]],
            ['ocean', 'gist_earth', 'terrain']) for me in mean_ene]
    if comp == 'all':
        return tuple(mean_ene)
    elif comp == 'nue':
        return mean_ene[0]
    elif comp == 'nua':
        return mean_ene[1]
    elif comp == 'nux':
        return mean_ene[2]
    
@hdf_isopen
@notrino_used
def neutrino_luminosity(self, file_name, 
                        comp:Literal['all', 'nue', 'nua', 'nux']='all',
                        **kwargs):
    # get radial flux
    if self.dim == 1 or kwargs['notrino']:
        flux = list(self.neutrino_momenta_grey(file_name, **kwargs))
    else:
        flux_nu = self.neutrino_momenta_grey(file_name, **kwargs)
        flux = [flux_nu[0][...,0], flux_nu[1][...,0], flux_nu[2][...,0]]
    #dA = self.cell.ax(self.ghost)

    lum = [f for f in flux]
    [l.set(name=s1, label=s2, cmap=s3, limits=[-1e49, 1e49]) 
     for l, s1, s2, s3 in zip(
        lum,
        ['Lum_nue', 'Lum_nua', 'Lum_nux'],
        [r'$L_{\nu_\mathrm{e}}$', r'$L_{\nu_\mathrm{a}}$', r'$L_{\nu_\mathrm{x}}$'],
        ['ocean', 'gist_earth', 'terrain']
    )]

    if comp == 'nue':
        return lum[0]
    elif comp == 'nua':
        return lum[1]
    elif comp == 'nux':
        return lum[2]
    elif comp == 'all':
        return lum
    
## GREY

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_energy_density_grey(self, file_name,
                                 comp:Literal['all', 'nue', 'nua', 'nux']='all',
                                 **kwargs):
    if not kwargs['notrino']:
        try:
            nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
                self._Simulation__data_h5['/neutrinogrey/egrey'][..., 0]), self.dim)
            nu_ene[..., 2] /= 4
        except:
            nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
                    self._Simulation__data_h5['neutrino/e'][..., 0]), self.dim)
            nu_ene[..., 2] /= 4
            nu_ene = np.sum(nu_ene, axis=self.dim)
    else:
        nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
            self._Simulation__data_h5['notrino/notrino_e'][...]), self.dim)
    if comp == 'all':
        return (aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens',
                    r'$E_{\nu_\mathrm{e}}$', 'viridis', [1e28, 1e32], True),
                aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens',
                    r'$E_{\overline{\nu}_\mathrm{e}}$', 'magma',
                    [1e28, 1e32], True),
                aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens',
                    r'$E_{\nu_\mathrm{x}}$', 'plasma', [1e28, 1e32], True)
                )
    elif comp == 'nue':
        return aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens',
                    r'$E_{\nu_\mathrm{e}}$', 'viridis', [1e28, 1e32], True)
    elif comp == 'nua':
        return aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens',
                    r'$E_{\overline{\nu}_\mathrm{e}}$', 'magma',
                    [1e28, 1e32], True)
    elif comp == 'nux':
        aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens',
                    r'$E_{\nu_\mathrm{x}}$', 'plasma', [1e28, 1e32], True)

@get_grid
@smooth
@hdf_isopen
@notrino_used
def neutrino_momenta_grey(self, file_name, **kwargs):
    """
    In the comoving rest frame of the fluid are equal to the 
    neutrino energy fluxes
    """
    if not kwargs['notrino']:
        try:
            nu_flux =  self.ghost.remove_ghost_cells(np.squeeze(
                self._Simulation__data_h5['/neutrinogrey/egrey'][..., 1:]), self.dim)
        except:
            nu_flux = self.ghost.remove_ghost_cells(np.squeeze(
                    self._Simulation__data_h5['neutrino/e'][..., 1:]), self.dim)
            nu_flux = np.sum(nu_flux, axis=self.dim)
    else:
        nu_flux = self.ghost.remove_ghost_cells(np.squeeze(
                self._Simulation__data_h5['notrino/notrino_f'][...]), self.dim)
        
    if self.dim == 1 or kwargs['notrino']:
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

@get_grid
@smooth
def neutrino_number_density_grey(self, file_name,
                                 comp:Literal['all', 'nue', 'nua', 'nux']='all',
                                 **kwargs):
    ndens = self.neutrino_number_density(file_name)
    final_nd = []
    for nd in ndens:
        lm, lb, cm, nm, lg  = nd.limits, nd.label, nd.cmap, nd.name, nd.log
        nnd = nd.sum(axis=-1)
        nnd.set(limits=lm, label=lb, cmap=cm, name=nm, log=lg)
        final_nd.append(nnd)
    if comp == 'all':
        return final_nd
    elif comp == 'nue':
        return final_nd[0]
    elif comp == 'nua':
        return final_nd[1]
    elif comp == 'nux':
        return final_nd[2]
    
# BLACK BODY

@get_grid
@smooth
@hdf_isopen
def neutrino_number_density_BB(self, file_name, **kwargs):
    ndens = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['notrino/notrino_nBB'][...]), self.dim)
    ndens = np.expand_dims(ndens, axis = -2)
    return (aerray(ndens[..., 0], u.cm ** (-3), 'Nnue_BB',
                r'$N_{\nu_\mathrm{e}}^\mathrm{BB}$', 'gnuplot', [1e33, 1e36], 
                True),
            aerray(ndens[..., 1], u.cm ** (-3), 'Nnua_BB',
                r'$N_{\overline{\nu}_\mathrm{e}}^\mathrm{BB}$', 'gnuplot_2',
                [1e31, 1e34], True),
            aerray(ndens[..., 2], u.cm ** (-3), 'Nnux_BB',
                r'$N_{\nu_\mathrm{x}}^\mathrm{BB}$', 'CMRmap', [1e32, 1e35], 
                True)
            )
    
@get_grid
@smooth
@hdf_isopen
def neutrino_energy_density_BB(self, file_name, **kwargs):
    """
    Neutrino energy density, black body equivalent.
    Now with NOTRINO case!
    """
    nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
    self._Simulation__data_h5['notrino/notrino_eBB'][...]), self.dim)
    nu_ene = np.expand_dims(nu_ene, axis = -2)
    nu_ene[..., 2] /= 4
    return (aerray(nu_ene[..., 0], u.erg / u.cm ** 3, 'nue_edens_BB',
                  r'$E_{\nu_\mathrm{e}}^\mathrm{BB}$', 'viridis', [1e28, 1e32], 
                  True),
            aerray(nu_ene[..., 1], u.erg / u.cm ** 3, 'nua_edens_BB',
                  r'$E_{\overline{\nu}_\mathrm{e}}^\mathrm{BB}$', 'magma',
                  [1e28, 1e32], True),
            aerray(nu_ene[..., 2], u.erg / u.cm ** 3, 'nux_edens_BB',
                  r'$E_{\nu_\mathrm{x}}^\mathrm{BB}$', 'plasma', [1e28, 1e32], 
                  True)
            )
