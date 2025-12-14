from AeViz.simulation.methods import *
from AeViz.utils.snewpy.detection_rates import (generate_event_rate_series,
                                                av_transformations,
                                                detectors_list)

"""
Methods to bridge AeViz with SNEWPY. So far only the event rate is
supported.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

@smooth
@derive
@sum_tob
def neutrino_event_rate(self, los:Literal['avg', 'eq', 'pol']='avg',
                        distance =(10*u.kpc), dt=(1*u.ms), 
                        transformation:av_transformations='NoTransformation',
                        detector:detectors_list='superk', tob_corrected=True,
                        **kwargs):
    return generate_event_rate_series(self, los, distance, dt, transformation,
                                      detector)