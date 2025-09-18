import numpy as np


def data_covers_grid(skeleton, grid):
    """Checks if a given skeleton covers a given grid"""
    for coord in ["lon", "lat"]:
        if skeleton.edges(coord)[0] > grid.edges(coord)[0]:
            return False
        if skeleton.edges(coord)[-1] < grid.edges(coord)[-1]:
            return False
    return True


def identify_boundary_edges(boundary_mask: np.ndarray) -> list[str]:
    """Identifies which edges has some boundary points
    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """
    edges = []

    if boundary_mask.shape[1] > 2:
        n0 = 1
        n1 = -1
    else:
        n0 = 0
        n1 = None

    if np.any(boundary_mask[-1, n0:n1]):
        edges.append("N")

    if np.any(boundary_mask[0, n0:n1]):
        edges.append("S")

    if boundary_mask.shape[0] > 2:
        n0 = 1
        n1 = -1
    else:
        n0 = 0
        n1 = None

    if np.any(boundary_mask[n0:n1, -1]):
        edges.append("E")

    if np.any(boundary_mask[n0:n1, 0]):
        edges.append("W")

    return edges


def create_ordered_boundary_list(edge_list):
    """Gets all edges, but given in a continuous clockwise direction.
    If this is not possible (e.g. ['N','S']), then ampty list is returned."""
    full_list = ["N", "E", "S", "W"]
    for ind, edge in enumerate(full_list):
        if edge not in edge_list:
            full_list[ind] = ""
    full_array = np.array(full_list)

    if len(np.where(full_array == "")[0]) == 0:
        return full_array.tolist()

    ct = 0
    while (
        np.where(full_array == "")[0][-1] != len(np.where(full_array == "")[0]) - 1
    ) and ct < 5:
        full_array = np.roll(full_array, 1)
        ct += 1

    if ct > 4:
        print(
            f"No continuous boundary can be found for edges {edge_list}. Returning empy list."
        )
        return []

    full_array = full_array[full_array != ""]

    return full_array.tolist()


def get_coords_for_boundary_edges(
    edges: list, lon_edges: tuple[float, float], lat_edges: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Create coordinate vectors for clockwise running edges.
    Assumes that edges are clockwise and continuous, which is imposed by the
    function create_ordered_boundary_list.
    Empty list return empty arrays.
    """
    lon = []
    lat = []
    for edge in edges:
        if edge == "N":
            lon.append(lon_edges[0])
            lat.append(lat_edges[1])
        if edge == "S":
            lon.append(lon_edges[1])
            lat.append(lat_edges[0])
        if edge == "W":
            lon.append(lon_edges[0])
            lat.append(lat_edges[0])
        if edge == "E":
            lon.append(lon_edges[1])
            lat.append(lat_edges[1])

    if edges:
        edge = edges[-1]  # Close the loop
        if edge == "N":
            lon.append(lon_edges[1])
            lat.append(lat_edges[1])
        if edge == "S":
            lon.append(lon_edges[0])
            lat.append(lat_edges[0])
        if edge == "W":
            lon.append(lon_edges[0])
            lat.append(lat_edges[1])
        if edge == "E":
            lon.append(lon_edges[1])
            lat.append(lat_edges[0])

    return np.array(lon), np.array(lat)


def is_gridded(data: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> bool:
    if lon is None or lat is None:
        return False

    if data.shape == (len(lat), len(lon)):
        return True

    if len(data.shape) == 1 and len(lat) == data.shape[0] and len(lon) == data.shape[0]:
        return False

    raise Exception(
        f"Size of data is {data.shape} but len(lat) = {len(lat)} and len(lon) = {len(lon)}. I don't know what is going on!"
    )


def expand_area(
    lon: tuple[float, float],
    lat: tuple[float, float],
    expansion_factor: float,
    dlon: float = 0.0,
    dlat: float = 0.0,
) -> tuple[float, float, float, float]:
    """
    Expands a lon-lat bounding box with an expansion factor.
    expansion_factor = 1 does nothing, and 1.5 expands 50% both directions.

    give dlon, dlat to ensure a minimum exansion in case the expansion factor is not enough

    E.g. for lat=(60, 61)
    expansion_factor=1.2 gives (59.9, 61.1)
    expansion_factor=1.2 and dlat = 0.25 gives (59.75, 61.25)
    """

    expand_lon = (lon[1] - lon[0]) * (expansion_factor - 1) * 0.5
    expand_lat = (lat[1] - lat[0]) * (expansion_factor - 1) * 0.5

    expand_lon = np.maximum(expand_lon, dlon)
    expand_lat = np.maximum(expand_lat, dlat)

    new_lon = lon[0] - expand_lon, lon[1] + expand_lon
    new_lat = lat[0] - expand_lat, lat[1] + expand_lat

    return new_lon, new_lat


def get_coordinates_from_ds(ds, return_dict: bool = False) -> tuple:
    """Determins if an xarray dataset is cartesian (x,y) or spherical (lon,lat)
    and returns the vecotrs (None for the ones that are not defined).

    If lon, lat is defined over time, the firs instance is grabbed."""

    if "time" in ds.dims:
        ds = ds.isel(time=0)

    if hasattr(ds, "lon") and hasattr(ds, "lat"):
        lon, lat = np.squeeze(ds.lon.values), np.squeeze(ds.lat.values)
        x, y = None, None

    if hasattr(ds, "longitude") and hasattr(ds, "latitude"):
        lon, lat = np.squeeze(ds.longitude.values), np.squeeze(ds.latitude.values)
        x, y = None, None
    if hasattr(ds, "x") and hasattr(ds, "y"):
        x, y = np.squeeze(ds.x.values), np.squeeze(ds.y.values)
        lon, lat = None, None

    if all_none([x, y, lon, lat]):
        raise AttributeError(
            "Dataset doesn't have a combination of lon(gitude)/lat(itude) or x/y!"
        )

    if return_dict:
        return {"lon": lon, "lat": lat, "x": x, "y": y}
    else:
        return lon, lat, x, y


def all_none(val) -> bool:
    return not [a for a in val if a is not None]
