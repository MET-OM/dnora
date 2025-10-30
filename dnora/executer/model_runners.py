from subprocess import Popen
from abc import ABC, abstractmethod
from dnora.file_module import FileNames
from dnora.type_manager.model_formats import ModelFormat
from dnora.type_manager.dnora_types import DnoraFileType
from .post_processors import PostProcessor, SwashMatToNc, HosOceanToNc, SwanMatToNc
from dnora import msg
import shutil
import os
import json
import xarray as xr
import numpy as np
import pandas as pd

from dnora.type_manager.dnora_objects import dnora_objects
from dnora.type_manager.dnora_types import data_type_from_string


class ModelRunner(ABC):
    """Runs the model."""

    def __init__(self):
        pass

    @abstractmethod
    def preferred_format(self) -> str:
        """For the file format using defauts.py, e.g. ModelFormat.SWAN"""
        return

    def post_processors(self) -> list[PostProcessor]:
        return []

    @abstractmethod
    def __call__(
        self,
        input_file: str = None,
        model_folder: str = None,
    ) -> None:
        """Runs the model executable"""

        return


class SWAN(ModelRunner):
    def __init__(self):
        self._post_processors = []

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.SWAN

    def post_processors(self) -> list[PostProcessor]:
        return self._post_processors

    def __call__(self, file_object, model_folder, nproc=4, **kwargs) -> None:
        if model_folder:
            msg.plain(f"{model_folder}/swanrun >>> {file_object.get_folder()}/swanrun")
            shutil.copy2(
                f"{model_folder}/swanrun", f"{file_object.get_folder()}/swanrun"
            )
            msg.plain(
                f"{model_folder}/swan.exe >>> {file_object.get_folder()}/swan.exe"
            )
            shutil.copy2(
                f"{model_folder}/swan.exe", f"{file_object.get_folder()}/swan.exe"
            )

        # Set post processor if mat-output is requested
        with open(f"{file_object.get_folder()}/{file_object.get_filename()}", "r") as f:
            line = True
            file_names_to_convert = []
            while line:
                line = f.readline()
                if ".mat" in line:
                    mat_file = [f for f in line.split(" ") if ".mat" in f][0]
                    file_names_to_convert.append(mat_file)

            if file_names_to_convert:
                msg.info(
                    "Mat-file output requested. Setting up post-processor to convert to Netcdf..."
                )
                self._post_processors = [SwanMatToNc(file_names_to_convert)]

        try:
            msg.plain(f"swanrun -input {file_object.get_filename()} -omp {nproc}")
            p = Popen(
                ["swanrun", "-input", file_object.get_filename(), "-omp", f"{nproc}"],
                cwd=file_object.get_folder(),
            )
            p.wait()
        except FileNotFoundError:
            msg.advice(
                "swanrun not found! 1) Set the environmental variable DNORA_SWAN_PATH=/path/to/swanfolder, 2) provide keyword model_folder=/path/to/swanfolder or 3) add the SWAN directiory to your environmental PATH"
            )

        return


class SWASH(ModelRunner):
    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.SWASH

    def post_processors(self) -> list[PostProcessor]:
        return [SwashMatToNc()]

    def __call__(self, file_object: FileNames, model_folder: str) -> None:
        try:
            p = Popen(
                [f"{model_folder}/swashrun", "-input", file_object.get_filename()],
                cwd=file_object.get_folder(),
            )
        except FileNotFoundError:
            msg.advice(
                "swashrun not found! 1) Set the environmental variable DNORA_SWASH_PATH=/path/to/swashfolder, 2) provide keyword model_folder=/path/to/swashfolder or 3) add the SWASH directiory to your environmental PATH"
            )

        p.wait()


WW3_DEFAULT_INPUTFILE_NAMES = {
    DnoraFileType.SPECTRA: "ww3_bounc.nml",
    DnoraFileType.GRID: "ww3_grid.nml",
    DnoraFileType.WIND: "ww3_prnc.nml",
    DnoraFileType.INPUT: "ww3_shel.nml",
    DnoraFileType.CURRENT: "ww3_prnc.nml",
    DnoraFileType.WATERLEVEL: "ww3_prnc.nml",
    DnoraFileType.ICE: "ww3_prnc.nml",
}


class WW3(ModelRunner):
    def __init__(self, program: str):
        """E.g. program = 'grid' to run ww3_grid etc."""
        self.program = program
        if program == "shel":
            self._post_processors = [WW3("ounf"), WW3("ounp")]
        else:
            self._post_processors = []
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.WW3

    def post_processors(self) -> list[PostProcessor]:
        return self._post_processors

    def __call__(self, file_object, model_folder, nproc=4, **kwargs) -> None:

        # Copy model executables if needed
        if model_folder:
            from_file = f"{model_folder}/ww3_{self.program}"
            to_file = file_object.get_folder()
            msg.copy_file(from_file, to_file)
            shutil.copy(from_file, to_file)

        # Target input file, e.g. ww3_prnc.nml
        to_file = WW3_DEFAULT_INPUTFILE_NAMES.get(file_object.obj_type)

        # Written input files, e.g. 'ww3_prcn_wind.nml', or ['ww3_prnc_ice.nml.sic', 'ww3_prnc_ice.nml.sit']
        if self.program in ["shel", "grid", "bounc", "prnc"]:
            from_files = [file_object.get_filename()]
            out_files = [f"{file_object.get_folder()}/ww3_{self.program}"]
            if file_object.obj_type == DnoraFileType.ICE:
                from_files = [f"{from_files[0]}.sic", f"{from_files[0]}.sit"]
                out_files = [f"{out_files[0]}_sic", f"{out_files[0]}_sit"]
            out_files = [f"{of}.out" for of in out_files]

        for from_file, out_file in zip(from_files, out_files):
            from_path = f"{file_object.get_folder()}/{from_file}"
            if not os.path.exists(from_path):
                msg.plain(f"{from_path} not found. Skipping...")
                continue
            if from_file != to_file:

                msg.copy_file(from_file, to_file)
                shutil.copy(from_path, f"{file_object.get_folder()}/{to_file}")

            msg.info(f"Running ww3_{self.program}...")
            msg.to_file(out_file)
            with open(out_file, "w") as outfile:
                try:
                    p = Popen(
                        [f"ww3_{self.program}"],
                        cwd=file_object.get_folder(),
                        stdout=outfile,
                    )
                    p.wait()
                except FileNotFoundError:
                    msg.advice(
                        f"ww3_{self.program} not found! 1) Set the environmental variable DNORA_WW3_PATH=/path/to/swashfolder or 2) provide keyword model_folder=/path/to/ww3folder"
                    )

        return


class HOS_ocean(ModelRunner):
    def __init__(self):
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.HOS_OCEAN

    def post_processors(self) -> list[PostProcessor]:
        return [HosOceanToNc()]

    def __call__(self, input_file: str, model_folder: str) -> None:
        print("Running HOS_ocean------------------->>>>>>>>>>>>>>>>>>>>>>>>>>")
        p = Popen(["HOS-ocean"], cwd=model_folder)
        p.wait()


class REEF3D(ModelRunner):
    def __init__(self, nproc=1):
        self.nproc = nproc
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.REEF3D

    def __call__(self, input_file: str, model_folder: str) -> None:
        p = Popen(["DiveMESH"], cwd=model_folder)
        p.wait()

        if self.nproc == 1:
            p = Popen(["REEF3D"], cwd=model_folder)
        else:
            p = Popen(
                ["/usr/bin/mpirun -n " + str(self.nproc) + " --oversubscribe REEF3D"],
                cwd=model_folder,
            )
        p.wait()


class VesselIcingPreProcessor(ModelRunner):
    def __init__(self, data_type: str) -> None:
        self._data_type = data_type

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.VESSEL_ICING

    def __call__(self, file_object: FileNames, model_folder: str) -> None:
        with open(f"{file_object.get_folder()}/mi-fieldcalc.json", "r") as f:
            config = json.load(f)

        obj_type_common = data_type_from_string("grid")
        obj_common = dnora_objects.get(obj_type_common)
        common_grid = obj_common.from_ds(xr.open_dataset(config.get("grid")))

        if self._data_type == "grid":
            out_data = common_grid
        else:
            obj_type = data_type_from_string(self._data_type)
            obj = dnora_objects.get(obj_type)
            raw_data = obj.from_ds(xr.open_dataset(config.get(self._data_type)))
            out_data = raw_data.resample.grid(common_grid, verbose=True)

        if self._data_type in ["wavegrid", "waveseries"]:
            outfile = f"{file_object.get_folder()}/mi-fieldcalc_wave.nc"
        else:
            outfile = f"{file_object.get_folder()}/mi-fieldcalc_{self._data_type}.nc"
        msg.to_file(outfile)
        out_data.ds().to_netcdf(outfile)


class VesselIcing(ModelRunner):
    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.VESSEL_ICING

    def __call__(self, file_object: FileNames, model_folder: str) -> None:
        try:
            import mi_fieldcalc as mifc
        except ImportError:
            raise ImportError(
                "Icing model is not installed! Install it with 'conda install mi-fieldcalc'"
            )
        with open(file_object.get_filepath(), "r") as f:
            config = json.load(f)

        data = {key: xr.open_dataset(value) for key, value in config.items()}
        """V vesselIcingMincog(const V sal, const V wave, const V x_wind, const V y_wind, const V airtemp, const V rh,
                        const V sst, const V p, const V Pw, /*const V aice,*/ const V depth, const V vs,
                        const V alpha, const V zmin, const V zmax, const int alt)"""

        msg.plain(f"Calculating icing:")
        icing = []
        for time in data.get("wind").time.values:
            hs = data.get("wavegrid").sel(time=time).hs.values
            sss = data.get("ocean").sel(time=time).sss.values
            uwnd = data.get("wind").sel(time=time).u.values
            vwnd = data.get("wind").sel(time=time).v.values
            t2m = np.full(data.get("ocean").sel(time=time).sss.shape, -20.0)
            relh = np.full(data.get("ocean").sel(time=time).sss.shape, 80.0)
            sst = data.get("ocean").sel(time=time).sst.values
            eta = data.get("waterlevel").sel(time=time).eta.values
            tm = data.get("wavegrid").sel(time=time).tm01.values
            sic = data.get("ice").sel(time=time).sic.values
            topo = data.get("grid").topo.values
            msg.plain(f"{pd.to_datetime(time):%Y-%m-%d %H:%M}")
            msg.plain(
                f"\t<< Mean over grid: wind=({np.mean(uwnd):.2f}, {np.mean(vwnd):.2f}), t2m={np.mean(t2m):.2f}, r={np.mean(relh):.2f}, sst={np.mean(sst):.2f}, sss={np.mean(sss):.2f}, hs={np.mean(hs):.2f}, tm={np.mean(tm):.2f}, sic={np.mean(sic):.2f}, depth={np.mean(topo):.2f}, waterlevel={np.mean(eta):.2f}"
            )
            icing.append(
                mifc.vesselIcingMincog(
                    sss,
                    hs,
                    uwnd,
                    vwnd,
                    t2m,  # t2m
                    relh,  # r
                    sst,
                    eta,
                    tm,
                    sic,
                    topo,
                    5,  # vs
                    3.14,  # alpha
                    4,  # zmin
                    4,  # zmax
                    1,  # alt
                    -1e20,  # undef
                )
            )
            msg.plain(f"\t>> Mean icing over grid: {np.mean(icing):.4f}")

        # Add icing variable to Dataset
        coords = {
            "lon": data.get("wind").lon,
            "lat": data.get("wind").lat,
            "time": data.get("wind").time,
        }
        da = xr.DataArray(
            np.stack(icing, axis=0),
            dims=("time", "lat", "lon"),
            coords=coords,
            name="icing",
        )
        da.attrs["grid_mapping"] = "crs"
        da.attrs["long_name"] = "MINCOG icing"
        da.attrs["units"] = "mm/h"
        outfile = f"{file_object.get_folder()}/mi-fieldcalc.nc"
        msg.to_file(outfile)
        da.to_netcdf(outfile)
