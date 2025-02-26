from __future__ import annotations  # For TYPE_CHECKING

import numpy as np
import netCDF4

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dnora.modelrun.modelrun import ModelRun
    from dnora.file_module import FileNames

from dnora.type_manager.dnora_types import DnoraDataType
from dnora import msg
from dnora import file_module
from .data_writers import DataWriter
from dnora.type_manager.spectral_conventions import SpectralConvention
import warnings


class SpectraWriter(DataWriter):
    """Writes the boundary spectra to a certain file format.

    This object is provided to the .export_spectra() method.
    """

    _convention = None

    def convention(self) -> SpectralConvention:
        """Defines in which format the incoming spectra should be.

        The conventions to choose from are predetermined:

        OCEAN:    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        MATHVEC:  Mathematical convention in vector
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        WW3:      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        if isinstance(self._convention, str):
            self._convention = SpectralConvention[self._convention.upper()]
        return self._convention


class WW3(SpectraWriter):
    def __init__(
        self, convention=SpectralConvention.WW3, one_file: bool = None
    ) -> None:
        if one_file is not None:
            warnings.warn(
                "Set 'one_file' as a keyword in the 'export_spectra'-method instead!",
                DeprecationWarning,
                stacklevel=2,
            )
        self.one_file = one_file
        self._convention = convention

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        var_names: dict = None,
        one_file: bool = True,
        squeeze_lonlat: bool = False,
        **kwargs,
    ) -> list[str]:
        """You can change the variable names by giving a dictionary e.g. var_names = {'efth': 'SPEC'}
        To write one file for each lon/lat points, set one_file = False
        To write longitude/latitude withouth a time dimension, set squeeze_lonlat = True
        """

        msg.info("Writing WAVEWATCH-III netcdf-output")
        boundary = model.spectra()

        if self.one_file is not None:
            one_file = self.one_file

        if one_file:
            # clean = False keeps possible #LON/#LAT tags
            filename = file_object.get_filepath(clean=False)
            if len(boundary.inds()) == 1:
                filename = file_object.replace_placeholders(
                    filename, lon=boundary.lon()[0], lat=boundary.lat()[0]
                )

            output_files = file_module.clean_filename(filename)
            self.write_netcdf(
                boundary,
                output_files,
                n=None,
                var_names=var_names,
                one_file=one_file,
                squeeze_lonlat=squeeze_lonlat,
            )

        else:
            output_files = []
            for n in boundary.inds():
                output_file = file_object.get_filepath(
                    lon=boundary.lon()[n], lat=boundary.lat()[n]
                )
                output_files.append(output_file)

                msg.plain(f"Point {n} >> {output_file}")
                self.write_netcdf(
                    boundary,
                    output_file,
                    n=n,
                    var_names=var_names,
                    one_file=one_file,
                    squeeze_lonlat=squeeze_lonlat,
                )

        return output_files

    def write_netcdf(
        self,
        boundary,
        output_file: str,
        n: int,
        var_names: dict,
        one_file: bool,
        squeeze_lonlat: bool,
    ) -> None:
        """Writes WW3 compatible netcdf spectral output from a list containing xarray datasets."""
        var_names = var_names or {}
        time_var = var_names.get("time", "time")
        inds_var = var_names.get("inds", "station")
        freq_var = var_names.get("freq", "frequency")
        dirs_var = var_names.get("dirs", "direction")
        efth_var = var_names.get("efth", "efth")
        lon_var = var_names.get("lon", "longitude")
        lat_var = var_names.get("lat", "latitude")

        msg.plain("Using the following names in the netcdf-file:")
        msg.plain(f"'time' >> '{time_var}'")
        msg.plain(f"'inds' >> '{inds_var}'")
        msg.plain(f"'freq' >> '{freq_var}'")
        msg.plain(f"'dirs' >> '{dirs_var}'")
        msg.plain(f"'efth' >> '{efth_var}'")
        msg.plain(f"'lon' >> '{lon_var}'")
        msg.plain(f"'lat' >> '{lat_var}'")

        root_grp = netCDF4.Dataset(output_file, "w", format="NETCDF4")
        #################### dimensions
        root_grp.createDimension(time_var, None)
        if one_file:
            root_grp.createDimension(inds_var, len(boundary.inds()))
        else:
            root_grp.createDimension(inds_var, 1)
        root_grp.createDimension("string16", 16)
        root_grp.createDimension(freq_var, len(boundary.freq()))
        root_grp.createDimension(dirs_var, len(boundary.dirs()))

        #######################################################
        ####################### variables
        time = root_grp.createVariable(time_var, np.float64, (time_var,))
        station = root_grp.createVariable(inds_var, np.int32, (inds_var,))
        frequency = root_grp.createVariable(freq_var, np.float32, (freq_var,))
        direction = root_grp.createVariable(dirs_var, np.float32, (dirs_var,))
        efth = root_grp.createVariable(
            efth_var,
            np.float32,
            (
                time_var,
                inds_var,
                freq_var,
                dirs_var,
            ),
        )

        if squeeze_lonlat:
            lonlat_variables = (inds_var,)
        else:
            lonlat_variables = (time_var, inds_var)

        latitude = root_grp.createVariable(
            lat_var,
            np.float32,
            lonlat_variables,
        )
        longitude = root_grp.createVariable(
            lon_var,
            np.float32,
            lonlat_variables,
        )
        station_name = root_grp.createVariable(
            "station_name",
            "S1",
            (
                inds_var,
                "string16",
            ),
        )
        string16 = root_grp.createVariable("string16", np.int32, ("string16",))

        ########################## Attributes
        time.units = "seconds since 1970-01-01 00:00:00 UTC"
        time.calendar = "standard"
        time.standard_name = time_var
        time.axis = "T"

        station.long_name = "station id"
        station.axis = "X"

        frequency.units = "s-1"
        frequency.long_name = "frequency of center band"
        frequency.standard_name = "sea_surface_wave_frequency"
        frequency.globwave_name = freq_var
        frequency.valid_min = 0
        frequency.valid_max = 10
        frequency.axis = "Y"

        direction.units = "degree"
        direction.long_name = "sea surface wave to direction"
        direction.standard_name = "sea_surface_wave_to_direction"
        direction.globwave_name = "direction"
        direction.valid_min = 0
        direction.valid_max = 360
        direction.axis = "Z"

        longitude.units = "degree_east"
        longitude.long_name = "longitude"
        longitude.standard_name = "longitude"
        longitude.valid_min = -180
        longitude.valid_max = 180
        # longitude:_FillValue = 9.96921e+36f ;
        longitude.content = "TX"
        longitude.associates = "time station"

        latitude.units = "degree_north"
        latitude.long_name = "latitude"
        latitude.standard_name = "latitude"
        latitude.valid_min = -90
        latitude.valid_max = 90
        # latitude:_FillValue = 9.96921e+36f ;
        latitude.content = "TX"
        latitude.associates = "time station"

        station_name.long_name = "station name"
        station_name.content = "XW"
        station_name.associates = "station string16"

        station.long_name = "station id"
        station.axis = "X"

        string16.long_name = "station_name number of characters"
        string16.axis = "W"

        efth.long_name = "sea surface wave directional variance spectral density"
        efth.standard_name = "sea_surface_wave_directional_variance_spectral_density"
        efth.globwave_name = "directional_variance_spectral_density"
        efth.units = "m2 s rad-1"
        efth.scale_factor = 1
        efth.add_offset = 0
        efth.valid_min = 0
        # efth.valid_max = 1.0e+20
        # efth._FillValue = 9.96921e+36
        efth.content = "TXYZ"
        efth.associates = "time station frequency direction"
        #######################################################
        ############## Pass data
        time[:] = boundary.time().values.astype("datetime64[s]").astype("float64")
        frequency[:] = boundary.freq()
        direction[:] = boundary.dirs()

        if one_file:
            station[:] = boundary.inds()
            efth[:] = boundary.spec()
            if squeeze_lonlat:
                lonlat_shape = (len(boundary.lon()),)
            else:
                lonlat_shape = (len(boundary.time()), len(boundary.lon()))
            longitude[:] = np.full(lonlat_shape, boundary.lon(), dtype=float)
            latitude[:] = np.full(lonlat_shape, boundary.lat(), dtype=float)
        else:
            station[:] = 1
            efth[:] = boundary.spec(inds=[n])
            if squeeze_lonlat:
                lonlat_shape = (1,)
            else:
                lonlat_shape = (len(boundary.time()), 1)
            longitude[:] = np.full(lonlat_shape, boundary.lon()[n], dtype=float)
            latitude[:] = np.full(lonlat_shape, boundary.lat()[n], dtype=float)
        # longitude[:] = bnd_out.longitude.values
        # latitude[:] = bnd_out.latitude.values
        station_name[:] = 1

        root_grp.close()
        return


class SWAN(SpectraWriter):
    def convention(self) -> SpectralConvention:
        """Convention of spectra"""
        return SpectralConvention.MET

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        factor: float = 1e-4,
        **kwargs,
    ) -> str:
        filename = file_object.get_filepath()
        boundary = model.spectra()
        days = boundary.days(datetime=False)

        with open(filename, "w") as file_out:
            file_out.write("SWAN   1\n")
            file_out.write("$ Data produced by " + boundary.name + "\n")
            file_out.write("TIME\n")
            file_out.write("          1\n")
            file_out.write("LONLAT\n")
            file_out.write("          " + format(len(boundary.inds())) + "\n")
            for k in range(len(boundary.inds())):
                file_out.write(
                    "   "
                    + format(boundary.lon()[k], ".4f")
                    + "  "
                    + format(boundary.lat()[k], ".4f")
                    + "\n"
                )
            file_out.write("AFREQ\n")
            file_out.write("          " + str(len(boundary.freq())) + "\n")
            for l in range(len(boundary.freq())):
                file_out.write("   " + format(boundary.freq()[l], ".4f") + "\n")
            file_out.write("NDIR\n")
            file_out.write("          " + format(len(boundary.dirs())) + "\n")
            for m in range(len(boundary.dirs())):
                file_out.write("   " + format(boundary.dirs()[m], ".1f") + "\n")
            file_out.write("QUANT\n")
            file_out.write("          1\n")
            file_out.write("VaDens\n")
            file_out.write("m2/Hz/degr \n")
            file_out.write("-32767\n")
            # first day
            # msg.to_file(f"{output_path}")

            for day in days:
                msg.plain(day)
                times = boundary.time(time=slice(day, day))
                for tim in times:
                    time_stamp = (
                        str(tim).split("-")[0]
                        + str(tim).split("-")[1]
                        + str(tim).split("-")[2][:2]
                        + "."
                        + str(tim).split("-")[2][3:5]
                        + "0000\n"
                    )
                    file_out.write(time_stamp)
                    for n in range(len(boundary.x())):
                        file_out.write("FACTOR\n")
                        file_out.write(format(factor, "1.0E") + "\n")
                        S = boundary.spec(
                            time=slice(str(tim), str(tim)), inds=n
                        ).squeeze()

                        # SWAN uses m*m/Hz/deg normalization
                        np.savetxt(file_out, S * np.pi / (180 * factor), fmt="%-10.0f")
        return filename
