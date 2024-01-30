from __future__ import annotations

from abc import ABC, abstractmethod
from copy import copy
import os
import numpy as np
import pandas as pd

# Import objects
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    # from ...wnd.wnd_mod import Forcing
    # from ...grd.grd_mod import Grid
    # from ...bnd.bnd_mod import Boundary
    # from ...spc.spc_mod import Spectra
    # from ...wsr.wsr_mod import WaveSeries
    from ...modelrun.modelrun import ModelRun
    from ...file_module import FileNames

from ...aux_funcs import create_swan_segment_coords
from ... import msg
from ... import file_module

from .ww3_functions import (
    ww3_grid,
    ww3_prnc,
    ww3_specfile_list,
    ww3_bounc,
    ww3_shel,
    ww3_spectral_output_list,
)


class InputFileWriter(ABC):
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
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        calib_wind: float = 1.0,
        calib_wcap: float = 0.5000e-04,
        use_wind: bool = True,
        **kwargs,
    ) -> str:
        forcing = model.forcing()
        grid = model.grid()
        grid_path = exported_files["grid"][-1]
        forcing_path = exported_files["forcing"][-1]
        boundary_path = exported_files["boundary"][-1]

        filename = file_object.get_filepath()

        spec_lon, spec_lat = grid.output_points()
        spec_points = [(x, y) for x, y in zip(spec_lon, spec_lat)]

        if forcing is None and use_wind == True:
            msg.info(
                "No forcing object provided. Wind information will NOT be written to SWAN input file!"
            )
            use_wind = False

        # Define start and end times of model run
        DATE_START = model.time(crop=True)[0]
        DATE_END = model.time(crop=True)[-1]
        STR_START = DATE_START.strftime("%Y%m%d.%H%M%S")
        STR_END = DATE_END.strftime("%Y%m%d.%H%M%S")

        delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 5)
        delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 5)

        delta_Xf = np.round(np.abs(forcing.lon()[-1] - forcing.lon()[0]), 5)
        delta_Yf = np.round(np.abs(forcing.lat()[-1] - forcing.lat()[0]), 5)

        factor_wind = calib_wind * 0.001

        with open(filename, "w") as file_out:
            file_out.write("$************************HEADING************************\n")
            file_out.write("$ \n")
            file_out.write(" PROJ '" + grid.name() + "' 'T24' \n")
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
                + " CIRCLE 36 0.04 1.0 31 \n"
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
                + str((delta_X / (grid.nx() - 1)).round(6))
                + " "
                + str((delta_Y / (grid.ny() - 1)).round(6))
                + "\n"
            )
            file_out.write(
                "READINP BOTTOM 1 '" + grid_path.split("/")[-1] + "' 3 0 FREE \n"
            )
            file_out.write("$ \n")

            lons, lats = create_swan_segment_coords(
                grid.boundary_mask(), grid.lon_edges(), grid.lat_edges()
            )
            bound_string = "BOUNDSPEC SEGMENT XY"
            for lon, lat in zip(lons, lats):
                bound_string += f" {lon:.2f} {lat:.2f}"
            bound_string += " VARIABLE FILE 0 "
            bound_string += f"'{boundary_path.split('/')[-1]}'\n"
            file_out.write(bound_string)
            # file_out.write('BOU NEST \''+boundary_path.split('/')[-1]+'\' OPEN \n')
            file_out.write("$ \n")

            if use_wind:
                file_out.write(
                    "INPGRID WIND "
                    + str(forcing.lon()[0])
                    + " "
                    + str(forcing.lat()[0])
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
                    + str(factor_wind)
                    + "  '"
                    + forcing_path.split("/")[-1]
                    + "' 3 0 0 1 FREE \n"
                )
                file_out.write("$ \n")
            else:
                file_out.write("OFF QUAD \n")
            file_out.write("GEN3 WESTH cds2=" + str(calib_wcap) + "\n")
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
                + grid.name()
                + "_"
                + STR_START.split(".")[0]
                + ".nc"
                + "' & \n"
            )
            file_out.write(
                "LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT "
                + STR_START
                + " 1 HR \n"
            )
            file_out.write("$ \n")
            if spec_points:
                file_out.write("POINTS 'pkt' &\n")
                for lon, lat in spec_points:
                    file_out.write(str(lon) + " " + str(lat) + " &\n")
                file_out.write(
                    "SPECOUT 'pkt' SPEC2D ABS '"
                    + grid.name()
                    + "_"
                    + STR_START.split(".")[0]
                    + "_spec.nc"
                    + "' & \n"
                )
                file_out.write("OUTPUT " + STR_START + " 1 HR \n")
            else:
                pass
            file_out.write("COMPUTE " + STR_START + " 10 MIN " + STR_END + "\n")
            file_out.write("STOP \n")

        return filename


class SWASH(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        bound_side_command: str = "BOU SIDE W CCW CON REG 0.5 14 270 ",
        **kwargs,
    ) -> str:
        grid = model.grid()
        filename = file_object.get_filepath()
        grid_path = exported_files["grid"][-1]
        DATE_START = model.time(crop=True)[0]
        DATE_END = model.time(crop=True)[-1]
        STR_START = pd.Timestamp(DATE_START).strftime("%H%M%S")
        STR_END = pd.Timestamp(DATE_END).strftime("%H%M%S")

        delta_X = np.round(np.abs(grid.edges("lon")[-1] - grid.edges("lon")[0]), 8)
        delta_Y = np.round(np.abs(grid.edges("lat")[-1] - grid.edges("lat")[0]), 8)

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
            file_out.write(self.bound_side_command + " \n")
            # file_out.write('BOU NEST \''+add_folder_to_filename(self.bnd_filename, self.bnd_folder)+'\' OPEN \n')
            file_out.write("$ \n")
            file_out.write("$*******************************************************\n")
            file_out.write("$ OUTPUT REQUESTS \n")
            temp_list = grid_path.split("/")
            forcing_folder = "/".join(temp_list[0:-1])
            file_out.write("BLOCK 'COMPGRID' NOHEAD '" + grid.name + ".mat" + "' & \n")
            file_out.write("LAY 3 WATL BOTL OUTPUT " + STR_START + " 5 SEC \n")
            file_out.write("$ \n")
            file_out.write("COMPUTE " + STR_START + " 0.001 SEC " + STR_END + "\n")
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


class WW3Grid(InputFileWriter):
    # def __init__(self):
    #     self.scaling = 10**6
    #     return

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        folder_on_server: str = "/server/gridfiles/",
        **kwargs,
    ) -> str:
        grid = model.grid()
        spectral_grid = model.spectral_grid()
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
        ww3_grid(
            grid,
            filename,
            exported_files["grid"],
            folder_on_server,
            freq1,
            nth,
            nk,
            dirshift,
        )

        return filename


class WW3Forcing(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        folder_on_server: str = "/server/windfiles/",
        **kwargs,
    ) -> str:
        if file_object.get_filename() == "":
            filename = file_object.get_folder() + "/ww3_prnc.nml"
        else:
            filename = file_object.get_filepath()
        ww3_prnc(filename, exported_files["forcing"], folder_on_server)

        return filename


class WW3Boundary(InputFileWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        exported_files: dict[str, list[str]],
        method: str = "nearest",
        verbose_level: int = 2,
        folder_on_server: str = "/server/boundaryfiles/",
        **kwargs,
    ) -> str:
        ww3_specfile_list(
            file_object.get_folder() + "/spectral_boundary_files.list",
            exported_files["boundary"],
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
        if homog is None:
            homog = {}
        lons, lats = model.grid().output_points()
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
        forcing["wnd"] = model.forcing() is not None
        forcing["wlv"] = model.waterlevel() is not None
        forcing["ocr"] = model.oceancurrent() is not None

        ww3_shel(filename, folder_on_server, start_time, end_time, forcing, homog)

        return filename
