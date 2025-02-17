from astropy import units as u

kBol = u.def_unit('kBol', 1.3807e-16 * u.cm ** 2 / u.g * u.s ** (-2) / u.K ,
                  'Boltzmann constant', {'latex': r'$k_\mathrm{B}$'})
bry = u.def_unit('bry', u.dimensionless_unscaled, 'Baryon', {'latex': 'bry'}, False)
