from astropy import units as u

new_units = dict()

kBol = u.def_unit('kBol', 1.3807e-16 * u.cm ** 2 / u.g * u.s ** (-2) / u.K ,
                  'Boltzmann constant', {'latex': r'$k_\mathrm{B}$'}, namespace=new_units)
bry = u.def_unit('bry', doc='Baryon', format={'latex': 'bry'},
                 namespace=new_units)
u.add_enabled_units(new_units)
setattr(u, 'kBol', kBol)
setattr(u, 'bry', bry)
u.set_enabled_equivalencies(u.dimensionless_angles())


__all__ = ['u']