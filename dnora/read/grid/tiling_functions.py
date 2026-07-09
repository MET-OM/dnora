from pathlib import Path
import numpy as np

def read_limits(filename: str) -> None:
    tiles = []

    file1 = open(Path(__file__).with_name(filename), "r")
    lines = file1.readlines()
    start_inds = range(0, len(lines), 5)
    lons = np.zeros((len(start_inds), 2))
    lats = np.zeros((len(start_inds), 2))
    for ind, n in enumerate(start_inds):
        tiles.append(lines[n].strip())
        lons[ind, 0] = float(lines[n + 1])
        lons[ind, 1] = float(lines[n + 2])
        lats[ind, 0] = float(lines[n + 3])
        lats[ind, 1] = float(lines[n + 4])

    return tiles, lons, lats
