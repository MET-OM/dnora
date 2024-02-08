from abc import ABC, abstractmethod
from dnora.modelrun import ModelRun
from dnora.file_module import FileNames
import scipy
import numpy as np
import xarray as xr


class PostProcessor(ABC):
    @abstractmethod
    def __call__(self, model: ModelRun, file_object: FileNames, **kwargs) -> None:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass


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

        input_file = f"{file_object.get_folder()}/{self.grid().name}.mat"
        output_file = f"{file_object.get_folder()}/{self.grid().name}.nc"
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
