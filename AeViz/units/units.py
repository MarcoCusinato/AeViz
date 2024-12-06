class units:
    """
    Class with the various conversion units and constants in cgs.
    """
    def __init__(self):
        self.solar_mass = 1.98847e33 #g
        self.speed_light = 2.99792458e10 #cm/s
        self.MeV = 1.360218e-6 #erg
        self.km = 1e5 #cm
        self.ms = 1e-3 #s
        self.G = 6.67e-8#cm^3/(g*s^2)
        self.L_cactus = self.G * self.solar_mass / self.speed_light ** 2
        self.T_cactus = self.G * self.solar_mass / self.speed_light ** 3

    def convert_to_solar_masses(self, quantity):
        return quantity / self.solar_mass
    
    def convert_to_grams(self, quantity):
        return quantity * self.solar_mass
    
    def convert_to_grams(self, quantity):
        return quantity * self.solar_mass

    def convert_in_c_units(self, quantity):
        return quantity / self.speed_light

    def convert_to_MeV(self, quantity):
        return quantity / self.MeV

    def convert_to_erg(self, quantity):
        return quantity * self.MeV

    def convert_to_km(self, quantity):
        return quantity / self.km

    def convert_to_cm(self, quantity):
        return quantity * self.km

    def convert_to_ms(self, quantity):
        return quantity / self.ms

    def convert_to_s(self, quantity):
        return quantity * self.ms
    
    def convert_to_grav_mass(self, mass, energy):
        return mass + self.convert_to_solar_masses(
            energy / self.speed_light ** 2)
    
    ## Cactus units
    ## c=G=M_sun=1
    
    def convert_to_cactus_lenght(self, quantity):
        return quantity / self.L_cactus
    
    def convert_from_cactus_lenght(self, quantity):
        return quantity * self.L_cactus
    
    def convert_to_cactus_time(self, quantity):
        return quantity / self.T_cactus
    
    def convert_from_cactus_time(self, quantity):
        return quantity * self.T_cactus
    
    def convert_to_cactus_rho(self, quantity):
        return quantity * self.L_cactus ** 3 / self.solar_mass
    
    def convert_from_cactus_rho(self, quantity):
        return quantity * self.solar_mass / self.L_cactus ** 3
    
    def convert_to_cactus_pressure(self, quantity):
        return self.convert_to_cactus_rho(quantity) / self.speed_light ** 2
    
    def convert_from_cactus_pressure(self, quantity):
        return self.convert_from_cactus_rho(quantity) * self.speed_light ** 2
    
    def convert_to_cactus_velocity(self, quantity):
        return quantity / self.speed_light
    
    def convert_from_cactus_velocity(self, quantity):
        return quantity * self.speed_light
