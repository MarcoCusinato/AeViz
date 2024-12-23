from AeViz.simulation.methods import *

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
    return nu_ene

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
        nu_flux = nu_flux[..., None]
    nu_flux[..., 2, :] /= 4
    return nu_flux

@smooth
@hdf_isopen
def neutrino_momenta_opacities(self, file_name, **kwargs):
    nu_opac = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['neutrino/oe'][..., 1:]), self.dim)
    if self.dim == 1:
        nu_opac = nu_opac[..., None]
    return nu_opac

@smooth
def neutrino_number_density(self, file_name, **kwargs):
    return self.neutrino_energy_density(file_name, **kwargs) / \
        u.convert_to_erg(self.cell.E_nu()[:, None])

@smooth
def neutrino_mean_energy(self, file_name, **kwargs):
    """
    Average neutrino energy per cell so 
    <e> = sum_w E_nu(w) / sum_w N_nu(w),
        with w the center of the neutrino bin
    """
    return u.convert_to_MeV(
        self.neutrino_energy_density(file_name, **kwargs).sum(axis=-2) 
                    / self.neutrino_number_density(file_name,
                                                    **kwargs).sum(axis=-2))
## GREY

@smooth
@hdf_isopen
def neutrino_energy_density_grey(self, file_name, **kwargs):
    nu_ene = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['/neutrinogrey/egrey'][..., 0]), self.dim)
    nu_ene[..., 2] /= 4
    return nu_ene

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
        nu_flux = nu_flux[..., None]
    nu_flux[..., 2, :] /= 4
    return nu_flux