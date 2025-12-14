"""
Dictionary containing all the informations for next and current
generation neutrino detectors. Most of the data have been found on the 
SNEWPY package
"""
from AeViz.units import u

class NeutrinoDetector:
    def __init__(self,
                 name,
                 location,
                 type,
                 mass,
                 ref_mass,
                 doi=None):
        self.name = name
        self.location = location
        self.det_type = type
        self.mass = mass * u.ktonne
        self.ref_mass = ref_mass * u.ktonne
        self.eff_fact = (self.mass/self.ref_mass).value
        self.rection_type = self.__infer_reaction(type)
        self.doi = doi
    
    def __infer_reaction(self, type):
        if type in ['wc100kt30prct']:
            return r'H$_2$O$/\overline{\nu}_e$'
        elif type in ['icecube', 'km3net']:
            return r'H$_2$O$/\overline{\nu}_e$ (String)'
        elif type in ['scint20kt']:
            return r'C$_n$H$_{2n}/\overline{\nu}_e$'
        elif type in ['halo1', 'halo2']:
            return r'Pb$\nu_e$'
        elif type in ['ar40kt']:
            return r'Ar/$\nu_e$'
        elif type in ['ds20']:
            return r'Ar/any $\nu$'
        elif type in ['xent', 'lz', 'pandax']:
            return r'Xe/any $\nu$'

superk = NeutrinoDetector("Super-Kamiokande",
                          "Japan",
                          "wc100kt30prct",
                          32,
                          100,
                          "10.1016/S0168-9002(03)00425-X")

hyper = NeutrinoDetector("Hyper-Kamiokande",
                         "Japan",
                         "wc100kt30prct",
                         220,
                         100,
                         "10.21468/SciPostPhysProc.17.019")

icecube = NeutrinoDetector("IceCube",
                           "Antarctica",
                           "icecube",
                           51600,
                           51600,
                           "10.21234/CPKQ-K003")

km3net = NeutrinoDetector("KM3NeT",
                          "UE",
                          "km3net",
                          69366 * 3,
                          69366 * 3,
                          "10.1016/j.nima.2014.05.090"
                          )

lvd = NeutrinoDetector("LVD",
                       "Italy",
                       "scint20kt",
                       1,
                       20,
                       "10.3204/DESY-PROC-2011-03/fulgione"
                       )

kamland = NeutrinoDetector("KamLAND",
                           "Japan",
                           "scint20kt",
                           1,
                           20,
                           "10.1103/PhysRevLett.94.081801"
                           )

borexino = NeutrinoDetector("Borexino",
                            "Italy",
                            "scint20kt",
                            0.278,
                            20,
                            "10.1103/PhysRevD.108.102005")

juno = NeutrinoDetector("JUNO",
                        "China",
                        "scint20kt",
                        20,
                        20,
                        "10.1088/0954-3899/43/3/030401")

snop = NeutrinoDetector("SNO+",
                        "Canada",
                        "scint20kt",
                        0.78,
                        20,
                        "10.48550/arXiv.0810.3694")

nova = NeutrinoDetector(r"NO$\nu$A",
                        "USA",
                        "novaFD",
                        14,
                        14,
                        "10.48550/arXiv.1309.7898")

halo = NeutrinoDetector("HALO",
                        "Canada",
                        "halo1",
                        0.240,
                        0.079,
                        "10.1088/1475-7516/2014/12/053")

halo1kt = NeutrinoDetector("HALO-1kt",
                           "Canada",
                           "halo2",
                           0.079,
                           1,
                           "10.1088/1475-7516/2014/12/053")

dune = NeutrinoDetector("DUNE",
                        "USA",
                        "ar40kt",
                        1,
                        1,
                        "10.48550/arXiv.2002.03005")

microboone = NeutrinoDetector("MicroBooNe",
                              "USA",
                              "ar40kt",
                              40,
                              40,
                              "10.48550/arXiv.1705.04894")

sbnd = NeutrinoDetector("SBND",
                        "USA",
                        "ar40kt",
                        0.09,
                        40,
                        "10.48550/arXiv.2203.05814")

baksan = NeutrinoDetector("Barksan",
                          "Russia",
                          "scint20kt",
                          0.12,
                          20,
                          "10.48550/arXiv.hep-ph/0409069")

darkside = NeutrinoDetector("DarkSide-20k",
                            "Italy",
                            "ds20",
                            0.0386,
                            0.0386,
                            "10.1038/s42005-024-01896-z")

xenonnt = NeutrinoDetector("XENONnT",
                           "Italy",
                           "xent",
                           0.006,
                           0.006,
                           "10.1103/PhysRevLett.133.191002")

lz = NeutrinoDetector("LZ",
                      "USA",
                      "lz",
                      0.007,
                      0.007,
                      "10.48550/arXiv.2512.08065")

pandax = NeutrinoDetector("PandaX-4T",
                          "China",
                          "pandax",
                          0.004,
                          0.004,
                          "10.48550/arXiv.2403.06220")