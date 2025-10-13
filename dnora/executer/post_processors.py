from abc import ABC, abstractmethod
from dnora.modelrun import ModelRun
from dnora.file_module import FileNames
import scipy
import numpy as np
import xarray as xr
import geo_parameters as gp
from geo_skeletons import GriddedSkeleton
from dnora import msg
import pathlib
from datetime import datetime

SWAN_VARS = {
    "Hsig": gp.wave.Hs,
    "Tm01": gp.wave.Tm01,
    "Tm02": gp.wave.Tm02,
    "Tm_10": gp.wave.Tm_10,
    "Depth": gp.ocean.WaterDepth,
    "Windv_x": gp.wind.XWind,
    "Windv_y": gp.wind.YWind,
    "Dspr": gp.wave.Spr,
    "Dir": gp.wave.Dirm,
    "TPsmoo": gp.wave.Tp("tps"),
    "RTpeak": gp.wave.Tp,
    "PkDir": gp.wave.Dirp,
    "Vel_x": gp.ocean.XCurrent,
    "Vel_y": gp.ocean.YCurrent,
}


class PostProcessor(ABC):
    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> None:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass


def read_swan_mat_to_ds(filename: str, lon: np.ndarray, lat=np.ndarray) -> xr.Dataset:
    """Reads a SWAN mat-file and creates an xarray Dataset with standard names"""
    # This creates a dict with keys
    # 'Hsig_20200215_000000', 'RTpeak_20200215_000000', 'TPsmoo_20200215_000000' ...
    mat = scipy.io.loadmat(filename)
    # Create class
    wave_grid = GriddedSkeleton.add_time()
    for key, value in SWAN_VARS.items():
        wave_grid = wave_grid.add_datavar(value(key))

    key = list(SWAN_VARS.keys())[0]
    keys = [a for a in mat.keys() if key == a.split("_")[0]]
    keys.sort()
    start_time = "T".join(keys[0].split("_")[1:])
    end_time = "T".join(keys[-1].split("_")[1:])
    if not start_time and not end_time:  # Stationary run
        start_time = f"{datetime.today():%Y-%m-%d %H:00}"
        end_time = f"{datetime.today():%Y-%m-%d %H:00}"

    obj = wave_grid(lon=lon, lat=lat, time=(start_time, end_time))

    for swan_var, param in SWAN_VARS.items():

        if swan_var in mat.keys():  # Stationary run
            keys = [swan_var]
        else:
            keys = [a for a in mat.keys() if swan_var == "_".join(a.split("_")[0:-2])]

        keys.sort()
        if keys:
            msg.plain(f"Decoding variable {swan_var} as {param}...")
        for n, kk in enumerate(keys):
            obj.ind_insert(
                swan_var,
                np.flipud(np.where(np.isnan(mat.get(kk)), np.nan, mat.get(kk))),
                time=n,
            )

    return obj.ds()


class SwanMatToNc(PostProcessor):
    def __init__(self, file_names_to_convert: list[str]):
        self._file_names_to_convert = file_names_to_convert

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> None:

        lon = model.grid().lon()
        lat = model.grid().lat()

        for filename in self._file_names_to_convert:
            filename = pathlib.Path(filename.replace("'", "")).stem

            input_file = f"{file_object.get_folder()}/{filename}.mat"
            output_file = f"{file_object.get_folder()}/{filename}.nc"

            msg.from_file(input_file)
            msg.blank()
            ds = read_swan_mat_to_ds(input_file, lon=lon, lat=lat)

            msg.blank()
            msg.to_file(output_file)
            ds.to_netcdf(output_file)

    def __str__(self):
        return "Converting SWAN output mat-file to NetCDF..."


class SwashMatToNc(PostProcessor):
    def __str__(self):
        return "Converting SWASH output mat-file to NetCDF..."

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> None:
        """
        Function to convert the swash mat-file (with only:'Botlev' & 'Wavelev')
        to a netcdf file
        Parameters:
        input_file = file.mat
        output_file = filename.nc
        ----------
        eta  : 3D free surface elevation [time,x, y]
        lon  : [lon_min, lon_max]
        lat  : [lat_min, lat_max]
        dt  : timestep of simulation in seconds"""

        input_file = f"{file_object.get_folder()}/{model.grid().name}.mat"
        output_file = f"{file_object.get_folder()}/{model.grid().name}.nc"
        dt = 1

        lon = model.grid().edges("lon", native=True)
        lat = model.grid().edges("lat", native=True)

        mat = scipy.io.loadmat(input_file)
        x = mat[list(mat.keys())[5]].shape[0]
        y = mat[list(mat.keys())[5]].shape[1]
        t = len(list(mat.keys())) - 4
        eta = np.zeros((t, x, y))
        depth = mat[list(mat.keys())[4]]  # mat['Botlev']
        lon = np.linspace(lon[0], lon[1], num=y)
        lat = np.linspace(lat[0], lat[1], num=x)
        time = np.linspace(0, t * dt, num=t)
        for i in range(5, len(list(mat.keys())), 1):
            eta[i - 5, :, :] = mat[list(mat.keys())[i]]

        # create xarray
        ds = xr.Dataset(
            {
                "eta": xr.DataArray(
                    eta,
                    coords={"time": time, "lat": lat, "lon": lon},
                    dims=["time", "lat", "lon"],
                    attrs={"units": "m", "long_name": "surface elevation"},
                ),
                "depth": xr.DataArray(
                    depth,
                    coords={"lat": lat, "lon": lon},
                    dims=["lat", "lon"],
                    attrs={"units": "m", "long_name": "depth"},
                ),
            }
        )

        ds.time.attrs["units"] = "s"
        # save xarray to netcdf
        ds.to_netcdf(output_file)


class HosOceanToNc(PostProcessor):
    def __str__(self):
        return "Converting HOS-Ocean output file to NetCDF..."

    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> None:
        """
            Function to convert the HOS-ocean file e.g., 3d.dat
            to a netcdf file
            Parameters:
            input_file = HOS-ocean output of eta, e.g., 3d.dat
            output_file = filename.nc
            ----------
            eta  : 3D free surface elevation [time,y, x]
        Returns
        -------
            ds : xr.Dataset
                eta: 3D surface elevation [time,y, x]
        """

        input_file = ""
        output_file = ""
        with open(input_file) as f:
            lines = f.readlines()

        # remove lines with # comments
        lines = [x for x in lines if not x.startswith("#")]
        lines = [s.replace("\n", "") for s in lines]

        I = int(lines[2].split(",")[1].split("=")[1])
        J = int(lines[2].split()[-1])

        title = lines[0].split('"')[1]

        # remove lines with titles
        lines = [x for x in lines if not x.startswith("T")]
        lines = [x for x in lines if not x.startswith("V")]
        lines = [x for x in lines if not x.startswith("Z")]
        lines = [s.split() for s in lines]

        timestep = int(len(lines) / (I * J))

        x = np.zeros(I * J)
        y = np.zeros(I * J)
        eta = np.zeros(len(lines))

        for i in range(len(lines)):
            if i < I * J:
                x[i] = float(lines[i][0])
                y[i] = float(lines[i][1])
                eta[i] = float(lines[i][2])
            else:
                eta[i] = float(lines[i][0])

        eta_3d = np.zeros((timestep, J, I))
        eta_3d[0, :, :] = eta[0 : I * J].reshape((-J, I))

        # fill x, y
        x = x[0:I]
        y = y[0::I]

        # fill eta_3d
        for t in range(timestep):
            eta_3d[t, :, :] = eta[t * (I * J) : (t + 1) * (I * J)].reshape((-J, I))

        # create xarray
        ds = xr.Dataset(
            {
                "eta": xr.DataArray(
                    eta_3d,
                    coords={"time": np.arange(timestep), "y": y, "x": x},
                    dims=["time", "y", "x"],
                    attrs={"units": "m", "long_name": title},
                )
            }
        )
        # save xarray ro netcdf
        ds.to_netcdf(output_file)
