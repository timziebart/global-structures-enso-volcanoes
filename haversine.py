
import numpy as np



EARTH_RADIUS = 6371 # km
EARTH_DIAMETER = 2 * EARTH_RADIUS
HALF_EARTH_CIRCUMFERENCE = EARTH_RADIUS * np.pi # km
EARTH_CIRCUMFERENCE = 2 * HALF_EARTH_CIRCUMFERENCE

def distance(origin, destination, 
        radius = EARTH_RADIUS
        ):
    if not ( isinstance(origin, np.ndarray) and isinstance(destination, np.ndarray) ):
        origin = np.array(origin)
        destination = np.array(destination)

    assert origin.shape[-1] == destination.shape[-1] == 2, "last axis for have to be longitude and latitude"
    if origin.shape[:-2] and destination.shape[:-1]:
        assert origin.shape == destination.shape, "origin and destination need to have the same shape, i.e. same number of points"

    origin = np.rollaxis(origin, -1)
    destination = np.rollaxis(destination, -1)

    lat1, lon1 = np.radians(origin)
    lat2, lon2 = np.radians(destination)

    origin = np.rollaxis(origin, 0, origin.ndim)
    destination = np.rollaxis(destination, 0, destination.ndim)

    dlat = lat2-lat1
    dlon = lon2-lon1
    a = np.sin(dlat/2) **2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2) **2
    c = 2 * np.arcsin(np.sqrt(a))
    d = radius * c

    return d
