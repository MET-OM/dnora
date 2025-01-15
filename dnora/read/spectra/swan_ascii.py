import numpy as np
import pandas as pd


def decode_lonlat(file) -> tuple:
    print(">>> Starting decoding of LONLAT")
    line = file.readline()
    n_locations = int(line.strip().split(" ")[0])
    print(f"Found {n_locations} locations...")
    lon = np.zeros(n_locations)
    lat = np.zeros(n_locations)
    for ct in range(n_locations):
        line = file.readline()
        lon[ct] = float(line.strip().split(" ")[0])
        lat[ct] = float(line.strip().split(" ")[-1])
        print(f"LOC{ct+1}: {lon[ct]} E, {lat[ct]} W")

    print("<<< Ending decoding of LONLAT")
    return lon, lat, file


def decode_afreq(file) -> tuple:
    print(">>> Starting decoding of AFREQ")
    line = file.readline()
    n_freq = int(line.strip().split(" ")[0])
    print(f"Found {n_freq} frequencies...")
    freq = np.zeros(n_freq)
    for ct in range(n_freq):
        line = file.readline()
        freq[ct] = float(line.strip().split(" ")[0])
        if ct == 0:
            print(f"First freq: {freq[ct]}")
        elif ct == n_freq - 1:
            print(f"Last freq: {freq[ct]}")

    print("<<< Ending decoding of AFREQ")
    return freq, file


def decode_ndir(file) -> tuple:
    print(">>> Starting decoding of NDIR")
    line = file.readline()
    n_dir = int(line.strip().split(" ")[0])
    print(f"Found {n_dir} directions...")
    dirs = np.zeros(n_dir)
    for ct in range(n_dir):
        line = file.readline()
        dirs[ct] = float(line.strip().split(" ")[0])
        if ct == 0:
            print(f"First dir: {dirs[ct]}")
        elif ct == n_dir - 1:
            print(f"Last dir: {dirs[ct]}")

    print("<<< Ending decoding of NDIR")
    return dirs, file


def decode_vadens(file, n_loc, n_freq, n_dir, start_time, end_time) -> tuple:
    print(">>> Starting decoding of VaDens")
    line = file.readline()  # Unit
    line = file.readline()  # Exception value

    times = []
    specs = []
    specs = []
    line = file.readline()
    if start_time is not None:
        start_time = pd.to_datetime(start_time)
    if end_time is not None:
        end_time = pd.to_datetime(end_time)
    while line:
        time_stamp = " ".join(line.strip().split(" ")[0].split("."))
        if start_time is not None and pd.to_datetime(time_stamp) < start_time:
            line = file.readline()
            continue
        if end_time is not None and pd.to_datetime(time_stamp) > end_time:
            break

        times.append(" ".join(line.strip().split(" ")[0].split(".")))
        if pd.to_datetime(times[-1]).hour == 0:
            print(pd.to_datetime(times[-1]).strftime("%Y-%m-%d"))

        single_spec = np.zeros((n_loc, n_freq, n_dir))
        no_data_spec = []

        for k in range(n_loc):
            line = file.readline()
            if "NODATA" in line:
                no_data_spec.append(k + 1)
            if "NODATA" in line or "ZERO" in line:
                single_spec[k, :, :] = 0
                continue

            assert "FACTOR" in line, f"Expected this line to be 'FACTOR'!"

            line = file.readline()
            factor = float(line)
            for n in range(n_freq):
                line = file.readline()
                single_spec[k, n, :] = (
                    np.array([int(a) for a in line.split(" ") if a]) * factor
                )
        specs.append(single_spec)
        line = file.readline()

    for i in no_data_spec:
        print(f"Spectrum {i}/{n_loc} had NODATA. (Set to 0)")

    print("<<< Ending decoding of VaDens")
    times = pd.to_datetime(times)
    specs = aa = np.stack(specs)

    return times, specs


def read_swan_ascii_spec(filename: str, start_time: str = None, end_time: str = None):
    with open(filename, "r") as file:
        line = file.readline()
        while line:
            line = file.readline()

            if "LONLAT" in line:
                lon, lat, file = decode_lonlat(file)
            if "AFREQ" in line:
                freq, file = decode_afreq(file)
            if "NDIR" in line:
                dirs, file = decode_ndir(file)
                dirs[dirs < 0] = dirs[dirs < 0] + 360
            if "VaDens" in line:
                time, spec = decode_vadens(
                    file, len(lon), len(freq), len(dirs), start_time, end_time
                )

    inds = np.argsort(dirs)
    spec = spec[:, :, :, inds]
    dirs = dirs[inds]
    return time, lon, lat, spec, freq, dirs
