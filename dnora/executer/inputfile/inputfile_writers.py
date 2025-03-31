from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
import re
import os
from pathlib import Path

# Import objects
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from dnora.modelrun.modelrun import ModelRun
    from dnora.file_module import FileNames

from dnora.utils.grid import create_swan_segment_coords
from dnora import msg, file_module
from dnora.type_manager.dnora_types import DnoraFileType
from .ww3_functions import (
    ww3_grid,
    ww3_prnc,
    ww3_specfile_list,
    ww3_bounc,
    ww3_shel,
    ww3_spectral_output_list,
)


class InputFileWriter(ABC):
    def file_type(self) -> DnoraFileType:
        return DnoraFileType.INPUT

    @abstractmethod
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        **kwargs,
    ) -> Union[str, list[str]]:
        pass


class Null(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        **kwargs,
    ) -> str:
        return ""


class SWAN(InputFileWriter):
    def __init__(
        self,
        timestep=10,
        f_low=0.04,
        f_high=1.0,
        n_freq=31,
        n_dir=36,
        output_var="HSIGN RTP TPS PDIR TM01 TMM10 DIR DSPR DEP",
    ):
        self.default_calibrations = {
            "wind": 1,
            "wcap": 0.5000e-04,
            "waterlevel": 1,
            "current": 1,
            "ice": 1,
        }
        self.swan_timestep = timestep
        self.f_low = f_low
        self.f_high = f_high
        self.n_freq = n_freq
        self.n_dir = n_dir
        self.output_var = output_var
        return

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        calibrate: dict[str, float] = None,
        use_wind: bool = True,
        use_waterlevel: bool = True,
        use_spectra: bool = True,
        use_current: bool = True,
        use_ice: bool = True,
    ):
        forcing = model.wind()
        boundary = model.spectra()
        waterlevel = model.waterlevel()
        oceancurrent = model.current()
        ice = model.ice()
        grid = model.grid()
        grid_path = exported_files["grid"][-1]
        forcing_path = exported_files["wind"][-1]
        boundary_path = exported_files["spectra"][-1]
        oceancurrent_path = exported_files["current"][-1]
        waterlevel_path = exported_files["waterlevel"][-1]
        ice_path = exported_files["ice"]

        filename = file_object.get_filepath()

        if calibrate is None:
            calibrate = {}

        if forcing is None and use_wind == True:
            msg.info(
                "No wind object provided. Wind information will NOT be written to SWAN input file!"
            )
            use_wind = False

        if waterlevel is None and use_waterlevel == True:
            msg.info(
                "No waterlevel object provided. Waterlevel information will NOT be written to SWAN input file!"
            )
            use_waterlevel = False

        if oceancurrent is None and use_current == True:
            msg.info(
                "No current object provided. OceanCurrent information will NOT be written to SWAN input file!"
            )
            use_current = False

        if boundary is None and use_spectra == True:
            msg.info(
                "No spectra object provided. Spectra information will NOT be written to SWAN input file!"
            )
            use_spectra = False

        if ice is None and use_ice == True:
            msg.info(
                "No ice object provided. Ice information will NOT be written to SWAN input file!"
            )
            use_ice = False

        # Define start and end times of model run
        DATE_START = model.time(crop_with="all")[0]
        DATE_END = model.time(crop_with="all")[-1]
        STR_START = DATE_START.strftime("%Y%m%d.%H%M%S")
        STR_END = DATE_END.strftime("%Y%m%d.%H%M%S")
        HOTSTART_FILE = (
            "hotstart_"
            + grid.name
            + "_"
            + (pd.Timestamp(DATE_START) - pd.Timedelta(hours=1)).strftime("%Y%m%d%H%M")
        )

        spec_lon, spec_lat = grid.output_points()

        # spec_points = [(x, y) for x, y in zip(spec_lon, spec_lat)]
        # STR_FORCING_START = STR_START
        # STR_FORCING_END = STR_END

        # # For wind forcing, range of all time steps in .asc file needs to be specified
        # if forcing is not None:
        #     if (forcing.name() == "ERA5"):  # bugfix: for era5 complete 24 hour timesteps are written to asc file.
        #         STR_FORCING_START = pd.Timestamp(DATE_START).strftime('%Y%m%d') + '.000000'
        #         STR_FORCING_END = pd.Timestamp(DATE_END).strftime('%Y%m%d') + '.230000'

        delta_X = np.round(np.diff(grid.edges("lon")), 5)[0]
        delta_Y = np.round(np.diff(grid.edges("lat")), 5)[0]

        factor = {}
        for calib_type in ["wind", "waterlevel", "current", "ice"]:
            factor[calib_type] = (
                calibrate.get(calib_type) or self.default_calibrations.get(calib_type)
            ) * 0.001
        factor["wcap"] = calibrate.get("wcap") or self.default_calibrations.get("wcap")

        with open(filename, "w") as file_out:
            file_out.write("$************************HEADING************************\n")
            file_out.write("$ \n")
            file_out.write(" PROJ '" + grid.name + "' 'T24' \n")
            file_out.write("$ \n")
            file_out.write("$*******************MODEL INPUT*************************\n")
            file_out.write("$ \n")
            file_out.write("SET NAUT \n")
            file_out.write("$ \n")
            file_out.write("MODE NONSTATIONARY TWOD \n")
            file_out.write("COORD SPHE CCM \n")
            file_out.write(
                "CGRID "
                + str(grid.lon()[0])
                + " "
                + str(grid.lat()[0])
                + " 0. "
                + str(delta_X)
                + " "
                + str(delta_Y)
                + " "
                + str(grid.nx() - 1)
                + " "
                + str(grid.ny() - 1)
                + " CIRCLE %d %f %f %d \n"
                % (self.n_dir, self.f_low, self.f_high, self.n_freq)
            )
            file_out.write("$ \n")

            file_out.write(
                "INPGRID BOTTOM "
                + str(grid.lon()[0])
                + " "
                + str(grid.lat()[0])
                + " 0. "
                + str(grid.nx() - 1)
                + " "
                + str(grid.ny() - 1)
                + " "
                + str((delta_X / (grid.nx() - 1)).round(8))
                + " "
                + str((delta_Y / (grid.ny() - 1)).round(8))
                + "\n"
            )
            file_out.write(
                "READINP BOTTOM 1 '" + grid_path.split("/")[-1] + "' 3 0 FREE \n"
            )
            file_out.write("$ \n")
            if use_spectra:
                lons, lats = create_swan_segment_coords(
                    grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
                )

                bound_string = "BOUNDSPEC SEGMENT XY"

                for lon, lat in zip(lons, lats):
                    bound_string += f" {lon:.4f} {lat:.4f}"
                bound_string += " VARIABLE FILE 0 "
                bound_string += f"'{boundary_path.split('/')[-1]}'\n"
                file_out.write(bound_string)

                # file_out.write('BOU NEST \''+boundary_path.split('/')[-1]+'\' OPEN \n')
                file_out.write("$ \n")

            if use_wind:
                self.output_var = self.output_var + " WIND"
                delta_Xf = np.round(np.diff(forcing.edges("lon")), 5)[0]
                delta_Yf = np.round(np.diff(forcing.edges("lat")), 5)[0]

                file_out.write(
                    "INPGRID WIND "
                    + str(forcing.lon()[0].round(3))
                    + " "
                    + str(forcing.lat()[0].round(3))
                    + " 0. "
                    + str(forcing.nx() - 1)
                    + " "
                    + str(forcing.ny() - 1)
                    + " "
                    + str((delta_Xf / (forcing.nx() - 1)).round(6))
                    + " "
                    + str((delta_Yf / (forcing.ny() - 1)).round(6))
                    + " NONSTATIONARY "
                    + STR_START
                    + f" {forcing.dt():.0f} HR "
                    + STR_END
                    + "\n"
                )

                file_out.write(
                    "READINP WIND "
                    + str(factor["wind"])
                    + "  '"
                    + forcing_path.split("/")[-1]
                    + "' 3 0 0 1 FREE \n"
                )
                file_out.write("$ \n")
            else:
                file_out.write("WIND 0 0 \n")  # no wind forcing

            if use_waterlevel:
                self.output_var = self.output_var + " WATLEV"
                delta_Xf = np.round(np.diff(waterlevel.edges("lon")), 5)[0]
                delta_Yf = np.round(np.diff(waterlevel.edges("lat")), 5)[0]

                file_out.write(
                    "INPGRID WLEV "
                    + str(waterlevel.lon()[0])
                    + " "
                    + str(waterlevel.lat()[0])
                    + " 0. "
                    + str(waterlevel.nx() - 1)
                    + " "
                    + str(waterlevel.ny() - 1)
                    + " "
                    + str((delta_Xf / (waterlevel.nx() - 1)).round(6))
                    + " "
                    + str((delta_Yf / (waterlevel.ny() - 1)).round(6))
                    + " NONSTATIONARY "
                    + STR_START
                    + f" {waterlevel.dt():.0f} HR "
                    + STR_END
                    + "\n"
                )

                file_out.write(
                    "READINP WLEV "
                    + str(factor["waterlevel"])
                    + "  '"
                    + waterlevel_path.split("/")[-1]
                    + "' 3 0 1 FREE \n"
                )
                file_out.write("$ \n")
            else:
                pass

            if use_current:
                self.output_var = self.output_var + " VEL"
                # delta_Xf = np.round(np.abs(oceancurrent.lon()[-1] - oceancurrent.lon()[0]), 5)
                # delta_Yf = np.round(np.abs(oceancurrent.lat()[-1] - oceancurrent.lat()[0]), 5)
                delta_Xf = np.round(np.diff(oceancurrent.edges("lon")), 5)[0]
                delta_Yf = np.round(np.diff(oceancurrent.edges("lat")), 5)[0]

                file_out.write(
                    "INPGRID CUR "
                    + str(oceancurrent.lon()[0].round(3))
                    + " "
                    + str(oceancurrent.lat()[0].round(3))
                    + " 0. "
                    + str(oceancurrent.nx() - 1)
                    + " "
                    + str(oceancurrent.ny() - 1)
                    + " "
                    + str((delta_Xf / (oceancurrent.nx() - 1)).round(6))
                    + " "
                    + str((delta_Yf / (oceancurrent.ny() - 1)).round(6))
                    + " EXC 32767 NONSTATIONARY "
                    + STR_START
                    + f" {oceancurrent.dt():.0f} HR "
                    + STR_END
                    + " \n"
                )

                file_out.write(
                    "READINP CUR "
                    + str(factor["current"])
                    + "  '"
                    + oceancurrent_path.split("/")[-1]
                    + "' 3 0 0 1 FREE \n"
                )
                file_out.write("$ \n")
            else:
                pass

            if use_ice:
                delta_Xf = np.round(np.diff(ice.edges("lon")), 5)[0]
                delta_Yf = np.round(np.diff(ice.edges("lat")), 5)[0]
                for i in range(len(ice_path)):
                    if ice_path[i].split("/")[-1].startswith("sic") or ice_path[
                        i
                    ].split("/")[-1].startswith("ice"):
                        ICE_NAME = "AICE"
                    elif ice_path[i].split("/")[-1].startswith("sit"):
                        ICE_NAME = "HICE"
                    # self.output_var = self.output_var + " " + ICE_NAME # comment due to AICE/HICE not available as output in nc-format in SWAN
                    file_out.write(
                        "INPGRID "
                        + ICE_NAME
                        + " "
                        + str(ice.lon()[0].round(3))
                        + " "
                        + str(ice.lat()[0].round(3))
                        + " 0. "
                        + str(ice.nx() - 1)
                        + " "
                        + str(ice.ny() - 1)
                        + " "
                        + str((delta_Xf / (ice.nx() - 1)).round(6))
                        + " "
                        + str((delta_Yf / (ice.ny() - 1)).round(6))
                        + " NONSTATIONARY "
                        + STR_START
                        + f" {ice.dt():.0f} HR "
                        + STR_END
                        + " \n"
                    )
                    file_out.write(
                        "READINP "
                        + ICE_NAME
                        + " "
                        + str(factor["ice"])
                        + "  '"
                        + ice_path[i].split("/")[-1]
                        + "' \n"
                    )
                    file_out.write("$ \n")
            else:
                pass

            if os.path.isfile(grid_path.split("/")[0] + "/" + HOTSTART_FILE) is True:
                file_out.write("INITIAL HOTSTART '" + HOTSTART_FILE + "'" "\n")

            file_out.write("GEN3 WESTH cds2=" + str(factor["wcap"]) + " AGROW" + "\n")
            file_out.write("FRICTION JON 0.067 \n")
            file_out.write("PROP BSBT \n")
            file_out.write("NUM ACCUR NONST 1 \n")
            file_out.write("$ \n")
            file_out.write("$*******************************************************\n")

            file_out.write("$ Generate block-output \n")
            temp_list = forcing_path.split("/")
            forcing_folder = "/".join(temp_list[0:-1])
            file_out.write(
                "BLOCK 'COMPGRID' HEAD '"
                + grid.name
                + "_"
                + STR_START.split(".")[0]
                + ".nc"
                + "' & \n"
            )
            file_out.write(
                "LAY 1 " + self.output_var + " OUTPUT " + STR_START + " 1 HR \n"
            )
            file_out.write("$ \n")
            if len(spec_lon) > 0:
                file_out.write("POINTS 'pkt' &\n")
                for slon, slat in zip(spec_lon, spec_lat):
                    file_out.write(str(slon) + " " + str(slat) + " &\n")
                file_out.write(
                    "SPECOUT 'pkt' SPEC2D ABS '"
                    + grid.name
                    + "_"
                    + STR_START.split(".")[0]
                    + "_spec.nc"
                    + "' & \n"
                )
                file_out.write("OUTPUT " + STR_START + " 1 HR \n")
            else:
                pass
            file_out.write(
                "COMPUTE "
                + STR_START
                + " %d MIN " % self.swan_timestep
                + STR_END
                + "\n"
            )
            file_out.write(
                "HOTFILE 'hotstart_"
                + grid.name
                + "_"
                + STR_END.replace(".", "")[:-2]
                + "'"
                + " FREE \n"
            )
            file_out.write("STOP \n")

        return filename


# class SWANold(InputFileWriter):
#     def __call__(
#         self,
#         model: ModelRun,
#         file_object: FileNames,
#         exported_files: dict[str, list[str]],
#         calib_wind: float = 1.0,
#         calib_wcap: float = 0.5000e-04,
#         use_wind: bool = True,
#         use_waterlevel: bool = True,
#         **kwargs,
#     ) -> str:
#         forcing = model.wind()
#         waterlevel = model.waterlevel()
#         grid = model.grid()
#         grid_path = exported_files["grid"][-1]
#         forcing_path = exported_files["wind"][-1]
#         boundary_path = exported_files["spectra"][-1]
#         waterlevel_path = exported_files["waterlevel"][-1]

#         filename = file_object.get_filepath()

#         spec_lon, spec_lat = grid.output_points()
#         spec_points = [(x, y) for x, y in zip(spec_lon, spec_lat)]

#         if forcing is None and use_wind:
#             msg.info(
#                 "No wind object provided. Wind information will NOT be written to SWAN input file!"
#             )
#             use_wind = False

#         if waterlevel is None and use_waterlevel:
#             msg.info(
#                 "No waterlevel object provided. Waterlevel information will NOT be written to SWAN input file!"
#             )
#             use_waterlevel = False

#         # Define start and end times of model run
#         DATE_START = model.time(crop_with="all")[0]
#         DATE_END = model.time(crop_with="all")[-1]
#         STR_START = DATE_START.strftime("%Y%m%d.%H%M%S")
#         STR_END = DATE_END.strftime("%Y%m%d.%H%M%S")

#         delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 5)
#         delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 5)

#         delta_Xf = np.round(np.abs(forcing.lon()[-1] - forcing.lon()[0]), 5)
#         delta_Yf = np.round(np.abs(forcing.lat()[-1] - forcing.lat()[0]), 5)

#         factor_wind = calib_wind * 0.001

#         delta_Xwl = np.round(np.abs(waterlevel.lon()[-1] - waterlevel.lon()[0]), 5)
#         delta_Ywl = np.round(np.abs(waterlevel.lat()[-1] - waterlevel.lat()[0]), 5)

#         with open(filename, "w") as file_out:
#             file_out.write("$************************HEADING************************\n")
#             file_out.write("$ \n")
#             file_out.write(" PROJ '" + grid.name + "' 'T24' \n")
#             file_out.write("$ \n")
#             file_out.write("$*******************MODEL INPUT*************************\n")
#             file_out.write("$ \n")
#             file_out.write("SET NAUT \n")
#             file_out.write("$ \n")
#             file_out.write("MODE NONSTATIONARY TWOD \n")
#             file_out.write("COORD SPHE CCM \n")
#             file_out.write(
#                 "CGRID "
#                 + str(grid.lon()[0])
#                 + " "
#                 + str(grid.lat()[0])
#                 + " 0. "
#                 + str(delta_X)
#                 + " "
#                 + str(delta_Y)
#                 + " "
#                 + str(grid.nx() - 1)
#                 + " "
#                 + str(grid.ny() - 1)
#                 + " CIRCLE 36 0.04 1.0 31 \n"
#             )
#             file_out.write("$ \n")

#             file_out.write(
#                 "INPGRID BOTTOM "
#                 + str(grid.lon()[0])
#                 + " "
#                 + str(grid.lat()[0])
#                 + " 0. "
#                 + str(grid.nx() - 1)
#                 + " "
#                 + str(grid.ny() - 1)
#                 + " "
#                 + str((delta_X / (grid.nx() - 1)).round(6))
#                 + " "
#                 + str((delta_Y / (grid.ny() - 1)).round(6))
#                 + "\n"
#             )
#             file_out.write(
#                 "READINP BOTTOM 1 '" + grid_path.split("/")[-1] + "' 3 0 FREE \n"
#             )
#             file_out.write("$ \n")

#             lons, lats = create_swan_segment_coords(
#                 grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
#             )
#             bound_string = "BOUNDSPEC SEGMENT XY"
#             for lon, lat in zip(lons, lats):
#                 bound_string += f" {lon:.2f} {lat:.2f}"
#             bound_string += " VARIABLE FILE 0 "
#             bound_string += f"'{boundary_path.split('/')[-1]}'\n"
#             file_out.write(bound_string)
#             # file_out.write('BOU NEST \''+boundary_path.split('/')[-1]+'\' OPEN \n')
#             file_out.write("$ \n")
#             dt = forcing.time().to_series().diff().dt.total_seconds().values[-1] / 3600
#             dt_wl = waterlevel.time().to_series().diff().dt.total_seconds().values[-1] / 3600
#             if use_wind:
#                 file_out.write(
#                     "INPGRID WIND "
#                     + str(forcing.lon()[0])
#                     + " "
#                     + str(forcing.lat()[0])
#                     + " 0. "
#                     + str(forcing.nx() - 1)
#                     + " "
#                     + str(forcing.ny() - 1)
#                     + " "
#                     + str((delta_Xf / (forcing.nx() - 1)).round(6))
#                     + " "
#                     + str((delta_Yf / (forcing.ny() - 1)).round(6))
#                     + " NONSTATIONARY "
#                     + STR_START
#                     + f" {dt:.0f} HR "
#                     + STR_END
#                     + "\n"
#                 )
#                 file_out.write(
#                     "READINP WIND "
#                     + str(factor_wind)
#                     + "  '"
#                     + forcing_path.split("/")[-1]
#                     + "' 3 0 0 1 FREE \n"
#                 )
#                 file_out.write("$ \n")
#             else:
#                 file_out.write("OFF QUAD \n")
#             if use_waterlevel:
#                 file_out.write(
#                     "INPGRID WLEVEL "
#                     + str(waterlevel.lon()[0])
#                     + " "
#                     + str(waterlevel.lat()[0])
#                     + " 0. "
#                     + str(waterlevel.nx() - 1)
#                     + " "
#                     + str(waterlevel.ny() - 1)
#                     + " "
#                     + str((delta_Xwl / (waterlevel.nx() - 1)).round(6))
#                     + " "
#                     + str((delta_Ywl / (waterlevel.ny() - 1)).round(6))
#                     + " NONSTATIONARY "
#                     + STR_START
#                     + f" {dt_wl:.0f} HR "
#                     + STR_END
#                     + "\n"
#                 )
#                 file_out.write(
#                     "READINP WLEVEL"
#                     + "  '"
#                     + waterlevel_path.split("/")[-1]
#                     + "' 3 0 0 1 FREE \n"
#                 )
#                 file_out.write("$ \n")
#             file_out.write("GEN3 WESTH cds2=" + str(calib_wcap) + "\n")
#             file_out.write("FRICTION JON 0.067 \n")
#             file_out.write("PROP BSBT \n")
#             file_out.write("NUM ACCUR NONST 1 \n")
#             file_out.write("$ \n")
#             file_out.write("$*******************************************************\n")
#             file_out.write("$ Generate block-output \n")
#             temp_list = forcing_path.split("/")
#             forcing_folder = "/".join(temp_list[0:-1])
#             file_out.write(
#                 "BLOCK 'COMPGRID' HEAD '"
#                 + grid.name
#                 + "_"
#                 + STR_START.split(".")[0]
#                 + ".nc"
#                 + "' & \n"
#             )
#             file_out.write(
#                 "LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT "
#                 + STR_START
#                 + " 1 HR \n"
#             )
#             file_out.write("$ \n")
#             if spec_points:
#                 file_out.write("POINTS 'pkt' &\n")
#                 for lon, lat in spec_points:
#                     file_out.write(str(lon) + " " + str(lat) + " &\n")
#                 file_out.write(
#                     "SPECOUT 'pkt' SPEC2D ABS '"
#                     + grid.name
#                     + "_"
#                     + STR_START.split(".")[0]
#                     + "_spec.nc"
#                     + "' & \n"
#                 )
#                 file_out.write("OUTPUT " + STR_START + " 1 HR \n")
#             else:
#                 pass
#             file_out.write("COMPUTE " + STR_START + " 10 MIN " + STR_END + "\n")
#             file_out.write("STOP \n")

#         return filename


class SWASH(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        boundary: str,
        bound_side_command: str = "BOU SIDE #BOUNDARY CCW CON REG #HS #TP #DIRP ",
        dt: float = 0.001,  # [s]
        **kwargs,
    ) -> str:
        grid = model.grid()
        filename = file_object.get_filepath()
        grid_path = exported_files["grid"][-1]

        DATE_START = model.time()[0]
        DATE_END = model.time()[-1]
        STR_START = pd.Timestamp(DATE_START).strftime("%H%M%S")
        STR_END = pd.Timestamp(DATE_END).strftime("%H%M%S")

        delta_X = np.round(np.abs(grid.edges("lon")[-1] - grid.edges("lon")[0]), 8)
        delta_Y = np.round(np.abs(grid.edges("lat")[-1] - grid.edges("lat")[0]), 8)

        waveseries = model.get("waveseries")

        hs = kwargs.get("hs")
        tp = kwargs.get("tp")
        dirp = kwargs.get("dirp")

        if waveseries is not None:
            if hs is None:
                msg.info(f"Using waveseries object {waveseries.name} for Hs!")
                hs = waveseries.hs(inds=0)[0]
            if tp is None:
                msg.info(f"Using waveseries object {waveseries.name} for Tp!")
                tp = waveseries.tp(inds=0)[0]
            if dirp is None:
                msg.info(f"Using waveseries object {waveseries.name} for Dirp!")
                dirp = waveseries.dirp(inds=0, dir_type="from")[0]
        else:
            hs = kwargs.get("hs") or 0.5
            tp = kwargs.get("tp") or 20
            dirp = kwargs.get("dirp") or 0

        msg.plain(f"Using hs = {hs:.2f}")
        msg.plain(f"Using tp = {tp:.1f}")
        msg.plain(f"Using dirp = {hs:.0f}")

        bound_side_command = re.sub("#HS", f"{hs:.2f}", bound_side_command)
        bound_side_command = re.sub("#TP", f"{tp:.1f}", bound_side_command)
        bound_side_command = re.sub("#DIRP", f"{dirp:.0f}", bound_side_command)
        bound_side_command = re.sub("#BOUNDARY", boundary, bound_side_command)

        msg.plain(f"Using bound_side_command: {bound_side_command}")

        with open(filename, "w") as file_out:
            file_out.write("$************************HEADING************************\n")
            file_out.write("$ \n")
            file_out.write(" PROJ '" + grid.name + "' 'T24' \n")
            file_out.write("$ \n")
            file_out.write("$*******************MODEL INPUT*************************\n")
            file_out.write("$ \n")
            file_out.write("SET NAUT \n")
            file_out.write("$ \n")
            file_out.write("MODE NONSTATIONARY TWOD \n")
            file_out.write("COORD SPHE CCM \n")
            file_out.write(
                "CGRID REG "
                + str(grid.edges("lon")[0])
                + " "
                + str(grid.edges("lat")[0])
                + " 0. "
                + str(delta_X)
                + " "
                + str(delta_Y)
                + " "
                + str(grid.nx() - 1)
                + " "
                + str(grid.ny() - 1)
                + " \n"
            )
            file_out.write("$ \n")
            file_out.write("VERT 1 \n")
            file_out.write("$ \n")
            file_out.write(
                "INPGRID BOTTOM "
                + str(grid.lon()[0])
                + " "
                + str(grid.lat()[0])
                + " 0. "
                + str(grid.nx() - 1)
                + " "
                + str(grid.ny() - 1)
                + " "
                + str((delta_X / (grid.nx() - 1)).round(8))
                + " "
                + str((delta_Y / (grid.ny() - 1)).round(8))
                + " EXC -999 \n"
            )
            file_out.write(
                "READINP BOTTOM 1 '" + grid_path.split("/")[-1] + "' 3 0 FREE \n"
            )
            file_out.write("$ \n")
            file_out.write(bound_side_command + " \n")
            # file_out.write('BOU NEST \''+add_folder_to_filename(self.bnd_filename, self.bnd_folder)+'\' OPEN \n')
            file_out.write("$ \n")
            file_out.write("$*******************************************************\n")
            file_out.write("$ OUTPUT REQUESTS \n")
            temp_list = grid_path.split("/")
            forcing_folder = "/".join(temp_list[0:-1])
            file_out.write("BLOCK 'COMPGRID' NOHEAD '" + grid.name + ".mat" + "' & \n")
            file_out.write("LAY 3 WATL BOTL OUTPUT " + STR_START + " 5 SEC \n")
            file_out.write("$ \n")
            file_out.write("COMPUTE " + STR_START + f" {dt:.3f} SEC " + STR_END + "\n")
            file_out.write("STOP \n")

        return filename


class REEF3D(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        option: str = "REEF3D",
        edges: list[str] = ["W"],
        nproc: int = 1,
        rot_angles: int = 0,
        wave_input: str = "SPEC1D",
        **kwargs,
    ) -> str:
        grid = model.grid()
        waveseries = model.waveseries()
        filename = file_object.get_filepath()
        grid_path = exported_files["grid"][-1]

        geodat = pd.read_csv(grid_path, sep=" ")  # read geo.dat
        geodat.columns = ["x", "y", "z"]

        if option == "DiveMESH":
            filename = "/".join(filename.split("/")[:-1]) + "/control.txt"
            with open(filename, "w") as file_out:
                if "W" in self.edges:
                    file_out.write("C 11 6 // West side: wave generation" "\n")
                    file_out.write("C 12 7 // side: numerical beach" "\n")
                    file_out.write("C 13 7 // side: numerical beach" "\n")
                    file_out.write("C 14 7 // side: numerical beach" "\n")
                elif "N" in self.edges:
                    file_out.write("C 12 6 // North side: wave generation" "\n")
                    file_out.write("C 13 3 // side: numerical beach" "\n")
                    file_out.write("C 14 3 // side: numerical beach" "\n")
                    file_out.write("C 11 7 // side: numerical beach" "\n")
                elif "E" in self.edges:
                    file_out.write("C 14 6 // East side: wave generation" "\n")
                    file_out.write("C 11 7 // side: numerical beach" "\n")
                    file_out.write("C 12 7 // side: numerical beach" "\n")
                    file_out.write("C 13 7 // side: numerical beach" "\n")
                elif "S" in self.edges:
                    file_out.write("C 13 6 // South side: wave generation" "\n")
                    file_out.write("C 11 7 // side: numerical beach" "\n")
                    file_out.write("C 12 7 // side: numerical beach" "\n")
                    file_out.write("C 14 7 // side: numerical beach" "\n")

                file_out.write("C 15 21 // bottom: wall boundary" "\n")
                file_out.write("C 16 3 // top: symmetry plane" "\n")
                file_out.write(" \n")

                file_out.write(
                    "B 1 "
                    + str(int(grid.dx().round(0)))
                    + " // horizontal mesh size dx"
                    "\n"
                )
                file_out.write(
                    "B 2 "
                    + str(int(geodat.x.max() / int(grid.dx().round(0))))
                    + " "
                    + str(int(geodat.y.max() / int(grid.dx().round(0))))
                    + " 10 // number of cells in x, y and z directions"
                    "\n"
                )
                file_out.write(
                    "B 10 0.0 "
                    + str(int(geodat.x.max()))
                    + " 0.0 "
                    + str(int(geodat.y.max()))
                    + " 0.0 1.0 // rectangular domain size"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write("B 103 5 // vertical grid clustering" "\n")
                file_out.write(
                    "B 113 2.5 // the stretching factor for the vertical grid clustering"
                    "\n"
                )
                file_out.write(
                    "B 116 1.0 // the focal point for the vertical grid clustering, which is water depth here"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write("G 10 1 // turn geodat on/off" "\n")
                file_out.write(
                    "G 13 "
                    + str(self.rot_angle)
                    + " // rotation angle of geo coordinates around vertical axis"
                    "\n"
                )
                file_out.write("G 15 2 // local inverse distance interpolation" "\n")
                file_out.write("G 20 0 // use automatic grid size off" "\n")
                file_out.write("G 31 14 // number of smoothing iterations" "\n")
                # file_out.write('G 41 1' '\n')
                file_out.write(" \n")

                file_out.write(
                    "M 10 " + str(self.nproc) + " // number of processors" "\n"
                )
                file_out.write("M 20 2 // decomposition method 2" "\n")
        elif option == "REEF3D":
            with open(filename, "w") as file_out:
                file_out.write("A 10 3  // choose the model reef::fnpf" "\n")
                file_out.write(
                    "A 310 3 // 3rd-order runge-kutta for fsfbc time treatment" "\n"
                )
                file_out.write(
                    "A 311 5 // 5th-order weno for fsfbc spatial treatment including wetting-drying"
                    "\n"
                )
                file_out.write("A 320 1 // 2nd-order laplace" "\n")
                file_out.write(" \n")

                file_out.write(
                    "A 341 2.0 // size of coastal relaxation zone by a factor of the horizontal cell size"
                    "\n"
                )
                file_out.write("A 343 1   // turn on wetting-drying" "\n")
                file_out.write(
                    "A 345 0.001 // wetting-drying water depth threshold" "\n"
                )
                file_out.write(
                    "A 346 2.1   // added viscosity within the coastal relaxation zone"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write(
                    "A 350 1 // viscosity damping wave breaking algorithm" "\n"
                )
                file_out.write(
                    "A 351 3 // breaking wave detection for both deep and shallow water"
                    "\n"
                )
                file_out.write(
                    "A 352 3 // additional filtering for viscosity based breaking for both deep and shallow water"
                    "\n"
                )
                file_out.write("A 361 5 // filtering outer iterations" "\n")
                file_out.write("A 362 2 // filtering inner iterations" "\n")
                file_out.write(
                    "A 365 1.86 // artificial viscosity for breaking wave energy dissipation"
                    "\n"
                )
                file_out.write(" \n")
                if wave_input == "SPEC1D":
                    file_out.write("B 85 10 // spectrum file" "\n")
                    file_out.write("B 90 1 // wave input" "\n")
                    file_out.write("B 92 31 // 1st-order irregular wave" "\n")

                elif wave_input == "JONSWAP":
                    file_out.write("B 85 2 // jonswap" "\n")
                    file_out.write("B 90 1 // wave input" "\n")
                    file_out.write("B 92 31 // 1st-order irregular wave" "\n")
                    file_out.write(
                        "B 93 "
                        + str(waveseries.hs()[0])
                        + " "
                        + str(waveseries.tp()[0])
                        + " // wave height, wave period"
                        "\n"
                    )
                    file_out.write(
                        "B 134 "
                        + str(waveseries.sprm()[0])
                        + " // spreading parameter for the directional spreading functions"
                        "\n"
                    )

                file_out.write(" \n")
                file_out.write(
                    "B 96 200.0 400.0 // wave generation zone length and numerical beach length"
                    "\n"
                )
                # file_out.write('B 107 0.0 '+str(int(geodat.x.max()))+' 0.0 0.0 200.0 // wave generation zone length and numerical beach length' '\n')
                # file_out.write('B 107 0.0 '+str(int(geodat.x.max()))+' '+str(int(geodat.y.max()))+' '+str(int(geodat.x.max()))+' 200.0 // customised numerical beach at the side walls' '\n')
                # file_out.write('B 107 25000.0 12000.0 0.0 16000.0 200.0 // customised numerical beach at the end of the tank' '\n')
                # file_out.write('B 107 0.0 0.0 2900.0 3500.0 200.0 // customised numerical beach at the side walls' '\n')
                # file_out.write('B 108 0.0 0.0 0.0 '+str(int(geodat.y.max()))+' 200.0 // customised wave generation zone' '\n')
                file_out.write("B 98 2 // relaxation method 2 for wave generation" "\n")
                file_out.write("B 99 2 // relaxation method 2 for numerical beach" "\n")
                file_out.write(" \n")

                file_out.write(
                    "F 60 " + str(grid.topo().max().round(1)) + " // still water depth"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write("G 50 1 // read in geo bathymetry" "\n")
                file_out.write(" \n")

                file_out.write(
                    "I 30 0 // turn off full tank initialisation, one can turn it on for a quick check of the setup"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write("N 41 1800.0 // simulation time" "\n")
                file_out.write("N 47 1.0 // cfl number" "\n")
                file_out.write(" \n")

                file_out.write("M 10 " + str(nproc) + " // number of processors" "\n")
                file_out.write(" \n")

                file_out.write("P 180 1 // turn on .vtp free surface printout" "\n")
                file_out.write(
                    "P 185 0.0 1800.0 0.5 // print out .vtp files interval based on simulation time window"
                    "\n"
                )
                file_out.write(" \n")

                file_out.write("W 22 -9.81 // gravity" "\n")
                file_out.write(" \n")

        return filename


class HOS_ocean(InputFileWriter):
    def __init__(
        self,
        n1=256,
        n2=64,
        xlen=None,
        ylen=80,
        T_stop=100,
        f_out=1,
        toler=1.0e-7,
        n=4,
        Ta=0,
        depth=100,
        Tp_real=10,
        Hs_real=4.5,
        gamma=3.3,
        beta=0.78,
        random_phases=1,
        tecplot=11,
        i_out_dim=1,
        i_3d=1,
        i_a_3d=0,
        i_2d=0,
        i_prob=0,
        i_sw=0,
    ):
        self.n1 = n1  # default is 256 at HOS-ocean-1.5/sources/HOS/variables_3D.f90
        self.n2 = n2  # default is 256 at HOS-ocean-1.5/sources/HOS/variables_3D.f90
        self.T_stop = T_stop
        self.f_out = f_out
        self.toler = toler
        self.n = n
        self.Ta = Ta
        self.depth = depth  # np.mean(self.topo()[self.land_sea_mask()])
        self.Tp_real = Tp_real
        self.Hs_real = Hs_real
        self.gamma = gamma
        self.beta = beta
        self.random_phases = random_phases
        self.tecplot = tecplot
        self.i_out_dim = i_out_dim
        self.i_3d = i_3d
        self.i_a_3d = i_a_3d
        self.i_2d = i_2d
        self.i_prob = i_prob
        self.i_sw = i_sw
        if xlen is None:
            self.xlen = (self.n1 * 9.81 * self.Tp_real**2) / (
                8 * 2 * np.pi
            )  # according to n1*2*np.pi/xlen = 5 k_p
            self.ylen = (self.n2 * 9.81 * self.Tp_real**2) / (8 * 2 * np.pi)
        else:
            self.xlen = xlen
            self.ylen = ylen
        return

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        **kwargs,
    ) -> str:
        # Create input file name
        filename = file_object.get_filepath()
        __, folder = file_module.split_filepath(filename)

        file_object.create_folder(folder=folder + "/Results")

        with open(filename, "w") as file_out:
            file_out.write("Restart previous computation :: i_restart        :: 0\n")
            file_out.write("Choice of computed case      :: i_case           :: 3\n")

            file_out.write("--- Geometry of the horizontal domain\n")
            file_out.write(
                "Length in x-direction        :: xlen             :: "
                + format(self.xlen, ".1f")
                + "\n"
            )
            file_out.write(
                "Length in y-direction        :: ylen             :: "
                + format(self.ylen, ".1f")
                + "\n"
            )

            file_out.write("--- Time stuff \n")
            file_out.write(
                "Duration of the simulation   :: T_stop           :: "
                + format(self.T_stop, ".1f")
                + "\n"
            )
            file_out.write(
                "Sampling frequency (output)  :: f_out            :: "
                + format(self.f_out, ".1f")
                + "\n"
            )
            file_out.write(
                "Tolerance of RK scheme       :: toler            :: "
                + format(self.toler, ".2e")
                + "\n"
            )
            file_out.write(
                "Dommermuth initialisation    :: n                :: "
                + format(self.n, ".0f")
                + "\n"
            )
            file_out.write(
                "Dommermuth initialisation    :: Ta               :: "
                + format(self.Ta, ".1f")
                + "\n"
            )

            file_out.write("--- Physical dimensional parameters \n")
            file_out.write("Gravity                      :: grav             :: 9.81\n")
            file_out.write(
                "Water depth                  :: depth            :: "
                + format(self.depth, ".1f")
                + "\n"
            )

            file_out.write("--- Irregular waves (i_case=3) \n")
            file_out.write(
                "Peak period in s             :: Tp_real          :: "
                + format(self.Tp_real, ".1f")
                + "\n"
            )
            file_out.write(
                "Significant wave height in m :: Hs_real          :: "
                + format(self.Hs_real, ".1f")
                + "\n"
            )
            file_out.write(
                "JONSWAP Spectrum             :: gamma            :: "
                + format(self.gamma, ".1f")
                + "\n"
            )
            file_out.write(
                "Directionality (Dysthe)      :: beta             :: "
                + format(self.beta, ".5f")
                + "\n"
            )
            file_out.write(
                "Random phases generation     :: random_phases    :: "
                + format(self.random_phases, ".0f")
                + "\n"
            )

            file_out.write("--- Output files \n")
            file_out.write(
                "Tecplot version              :: tecplot          :: "
                + format(self.tecplot, ".0f")
                + "\n"
            )
            file_out.write(
                "Output: 1-dim. ; 0-nondim.   :: i_out_dim        :: "
                + format(self.i_out_dim, ".0f")
                + "\n"
            )
            file_out.write(
                "3d free surface quantities   :: i_3d             :: "
                + format(self.i_3d, ".0f")
                + "\n"
            )
            file_out.write(
                "3d modes                     :: i_a_3d           :: "
                + format(self.i_a_3d, ".0f")
                + "\n"
            )
            file_out.write(
                "2d free surface, center line :: i_2d             :: "
                + format(self.i_2d, ".0f")
                + "\n"
            )
            file_out.write(
                "Wave probes in domain        :: i_prob           :: "
                + format(self.i_prob, ".0f")
                + "\n"
            )
            file_out.write(
                'Swense output 1="yes",0="no" :: i_sw             :: '
                + format(self.i_sw, ".0f")
                + "\n"
            )
        return filename


def apply_folder_on_server(
    exported_files: list[str], folder_on_server: str
) -> list[str]:
    """Replaces true export folder with a given folder on server"""
    obj_exported_to = []
    for grid_filename in exported_files:
        fn = Path(grid_filename).name
        obj_exported_to.append(str(Path(folder_on_server, fn)))

    return obj_exported_to


class WW3Grid(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        folder_on_server: str = "",
        **kwargs,
    ) -> str:
        grid = model.grid()
        spectral_grid = model.get("spectralgrid")
        if file_object.get_filename() == "":
            filename = file_object.get_folder() + "/ww3_grid.nml"
        else:
            filename = file_object.get_filepath()

        if spectral_grid is not None:
            freq1 = spectral_grid.freq()[0]
            nth = len(spectral_grid.dirs())
            nk = len(spectral_grid.freq())
            dirshift = spectral_grid.dirs()[0]
        else:
            freq1 = 0.04118
            nth = 36
            nk = 32
            dirshift = 0
        if folder_on_server:
            grid_exported_to = apply_folder_on_server(
                exported_files["grid"], folder_on_server
            )
        else:
            grid_exported_to = exported_files["grid"]

        ww3_grid(
            grid,
            filename,
            grid_exported_to,
            freq1,
            nth,
            nk,
            dirshift,
        )

        return filename


class WW3Wind(InputFileWriter):
    def file_type(self) -> DnoraFileType:
        return DnoraFileType.WIND

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        folder_on_server: str = "",
        **kwargs,
    ) -> str:
        if file_object.get_filename() == "":
            filename = file_object.get_folder() + "/ww3_prnc.nml"
        else:
            filename = file_object.get_filepath()
        if folder_on_server:
            wind_exported_to = apply_folder_on_server(
                exported_files["wind"], folder_on_server
            )
        else:
            wind_exported_to = exported_files["wind"]
        ww3_prnc(filename, wind_exported_to)

        return filename


class WW3Spectra(InputFileWriter):

    def file_type(self) -> DnoraFileType:
        return DnoraFileType.SPECTRA

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        method: str = "nearest",
        verbose_level: int = 2,
        folder_on_server: str = "",
        **kwargs,
    ) -> str:
        msg.to_file(file_object.get_folder() + "/spectral_boundary_files.list")
        spectra_exported_to = apply_folder_on_server(
            exported_files["spectra"], folder_on_server
        )
        ww3_specfile_list(
            file_object.get_folder() + "/spectral_boundary_files.list",
            spectra_exported_to,
        )

        if file_object.get_filename() == "":
            filename = file_object.get_folder() + "/ww3_bounc.nml"
        else:
            filename = file_object.get_filepath()
        if method == "nearest":
            method = 1
        if method == "linear":
            method = 2
        if method not in [1, 2]:
            msg.warning(
                f"Cannot undestand method={method} (not nearest/linear). Uning nearest neighbour."
            )
            method = 1

        ww3_bounc(filename, method, verbose_level)

        return filename


class WW3(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        homog: dict[tuple[float, float]] = None,
        folder_on_server: str = "/server/namelists/",
        **kwargs,
    ) -> str:
        """To use homogeneous input, set all the variables in order as: homog = {'wind': [1,4]}"""
        if homog is None:
            homog = {}
        lons, lats = model.grid().output_points()
        msg.to_file(file_object.get_folder() + "/spectral_points.list")
        ww3_spectral_output_list(
            file_object.get_folder() + "/spectral_points.list", lons, lats
        )
        start_time = model.start_time(crop_with="all").strftime("%Y%m%d %H0000")
        end_time = model.end_time(crop_with="all").strftime("%Y%m%d %H0000")
        if file_object.get_filename() == "":
            filename = file_object.get_folder() + "/ww3_shel.nml"
        else:
            filename = file_object.get_filepath()

        forcing = {}
        forcing["wind"] = model.wind() is not None
        forcing["waterlevel"] = model.waterlevel() is not None
        forcing["current"] = model.current() is not None

        ww3_shel(filename, folder_on_server, start_time, end_time, forcing, homog)

        return filename
