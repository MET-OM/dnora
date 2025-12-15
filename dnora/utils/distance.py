import geopy.distance
import numpy as np


def distance_2points(lat1, lon1, lat2, lon2) -> float:
    """Calculate distance between two points"""
    return geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).km


def min_distance(lon, lat, lon_vec, lat_vec) -> tuple[float, int]:
    """Calculates minimum distance [km] between a given point and a list of
    points given in spherical coordinates (lon/lat degrees).

    Also returns index of the found minimum.
    """
    dx = []
    for n, __ in enumerate(lat_vec):
        dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))

    return np.array(dx).min(), np.array(dx).argmin()


def min_cartesian_distance(x, y, x_vec, y_vec) -> tuple[float, int]:
    """ "Calculates minimum distance [km] between a given point and list of points given
    in cartesian coordinates [m].

    Also returns incex of found minimum"""
    dx = ((y - y_vec) ** 2 + (x - x_vec) ** 2) ** 0.5

    return dx.min() / 1000, dx.argmin()


def lon_in_km(lat: float) -> float:
    """Converts one longitude degree to km for a given latitude."""

    return distance_2points(lat, 0, lat, 1)


def domain_size_in_km(
    lon: tuple[float, float], lat: tuple[float, float]
) -> tuple[float, float]:
    """Calculates approximate size of grid in km."""

    km_x = distance_2points(
        (lat[0] + lat[1]) / 2, lon[0], (lat[0] + lat[1]) / 2, lon[1]
    )
    km_y = distance_2points(lat[0], lon[0], lat[1], lon[0])

    return km_x, km_y

def assert_lon_almost_equal(lon1, lon2, decimal=6):
    """
    Assert that two longitude arrays are almost equal,
    accounting for wraparound at ±180°.
    """
    lon1 = np.asarray(lon1)
    lon2 = np.asarray(lon2)

    diff = (lon1 - lon2 + 180) % 360 - 180
    np.testing.assert_array_almost_equal(diff, 0.0, decimal=decimal)
