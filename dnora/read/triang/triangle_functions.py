import numpy as np

from pathlib import Path


def read_node(filename: str) -> tuple[np.ndarray]:
    """Reads in node and possibly bathymetty and bounday information from .node file"""
    # Read in nodes
    with open(Path(filename).with_suffix(".node")) as f:
        n_of_nodes, n_of_coords, n_of_attributes, n_of_boundaries = list(
            map(int, f.readline().split())
        )
        if n_of_coords != 2:
            raise ValueError(f"Number of coordinates needs to be 2, not {n_of_coords}!")
        if n_of_attributes not in [0, 1]:
            raise ValueError(
                f"Number of attributes needs to be 0 (no depth) or 1 (depth provided), not {n_of_attributes}!"
            )
        if n_of_boundaries not in [0, 1]:
            raise ValueError(
                f"Number of boundaires needs to be 0 (no boundary) or 1 (boundary information included), not {n_of_boundaries}!"
            )

        x = np.zeros((n_of_nodes, 1))
        y = np.zeros((n_of_nodes, 1))
        boundaries = np.zeros((n_of_nodes, 1), int)
        bathymetry = np.zeros((n_of_nodes, 1))
        ct = 0
        for line in f:
            if n_of_attributes == 0:
                if n_of_boundaries == 0:
                    __, lon, lat = map(float, line.split())
                    bnd = 0
                    depth = np.nan
                else:
                    __, lon, lat, bnd = line.split()
                    lon = float(lon)
                    lat = float(lat)
                    bnd = int(bnd)
                    depth = np.nan
            else:
                if n_of_boundaries == 0:
                    __, lon, lat, depth = map(float, line.split())
                    bnd = 0
                else:
                    __, lon, lat, depth, bnd = line.split()
                    lon = float(lon)
                    lat = float(lat)
                    depth = float(depth)
                    bnd = int(bnd)
            x[ct] = lon
            y[ct] = lat
            boundaries[ct] = bnd
            bathymetry[ct] = depth
            ct += 1

    return (x, y, bathymetry, boundaries)


def read_ele(filename: str) -> np.ndarray:
    """Reads in triangle information from a .ele file"""
    # Read in elements (trianles)
    with open(Path(filename).with_suffix(".ele")) as f:
        n_of_triangles, n_of_nodes, n_of_attributes = list(
            map(int, f.readline().split())
        )
        if n_of_attributes > 0:
            raise ValueError(
                f"Expect to not get any attributes in an .ele-file, but header specifies {n_of_attributes}! Give possible bathymetry in the .node file or in separate .bot file."
            )
        if n_of_nodes != 3:
            raise ValueError(
                f"Expect number of nodes to be 3 (triangles), not {n_of_nodes}!"
            )

        tri = np.zeros((n_of_triangles, 3), int)
        ct = 0
        for line in f:
            tri[ct, :] = list(map(int, line.split()))[1:]
            ct += 1
        if np.min(tri) == 1:
            tri = tri - 1

    return tri
