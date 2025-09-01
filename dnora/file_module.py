from __future__ import annotations
from pathlib import Path
from .defaults import read_defaults
from dataclasses import dataclass
import pandas as pd
import re
from dnora import msg
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from dnora.modelrun import ModelRun
    from dnora.type_manager.dnora_objects import DnoraObject
from dnora.type_manager.dnora_types import DnoraDataType, DnoraFileType
from dnora.type_manager.model_formats import ModelFormat


@dataclass
class FileNames:
    model: ModelRun
    obj_type: Union[DnoraDataType, DnoraFileType] = None
    obj_name: str = None
    format: ModelFormat = ModelFormat.MODELRUN
    filename: str = None
    folder: str = None
    dateformat: str = None
    dateformat_folder: str = None
    edge_object: DnoraDataType = None
    time_object: DnoraDataType = None
    start_time: str = None
    end_time: str = None

    def __post_init__(self):
        self._defaults = read_defaults("export_defaults.yml", from_module=True)

        self.fallback = self._defaults[ModelFormat.MODELRUN.name]
        self.primary = self._defaults[self.format.name]

        if self.edge_object is None:
            if (
                isinstance(self.obj_type, DnoraDataType)
                and self.model.get(self.obj_type) is not None
            ):
                self.edge_object = self.obj_type
            else:
                self.edge_object = DnoraDataType.GRID  # Always present in ModelRun

        if self.obj_name is None and self.obj_type is not None:
            if self.model.get(self.obj_type) is not None:
                self.obj_name = self.model[self.obj_type].name

    def get_start_time(self):
        if self.start_time is not None:
            return self.start_time
        time_object = self.time_object or "all"
        return self.model.start_time(crop_with=time_object)

    def get_end_time(self):
        if self.end_time is not None:
            return self.end_time
        time_object = self.time_object or "all"
        return self.model.end_time(crop_with=time_object)

    def get_dateformat(self) -> str:
        return self.dateformat or get_default_value(
            "dateformat", self.obj_type, self.primary, self.fallback
        )

    def get_dateformat_folder(self) -> str:
        return self.dateformat_folder or get_default_value(
            "dateformat_folder", self.obj_type, self.primary, self.fallback
        )

    def get_filename(
        self,
        extension: str = None,
        start_time: str = None,
        end_time: str = None,
        lon: tuple[float, float] = None,
        lat: tuple[float, float] = None,
        x: tuple[float, float] = None,
        y: tuple[float, float] = None,
        clean: bool = True,
    ) -> str:
        if self.filename is not None:
            filename = self.filename
        else:
            filename = get_default_value(
                "filename", self.obj_type, self.primary, self.fallback
            )
        start_time = start_time or self.get_start_time()
        end_time = end_time or self.get_end_time()
        lon = lon or self.model[self.edge_object].edges("lon", strict=True)
        lat = lat or self.model[self.edge_object].edges("lat", strict=True)
        x = x or self.model[self.edge_object].edges("x", strict=True)
        y = y or self.model[self.edge_object].edges("y", strict=True)

        filename = self.replace_placeholders(
            filename, start_time, end_time, lon, lat, x, y
        )
        filename_suffix = Path(filename).suffix[1:]
        # pytest_ww3_spectra_00.0000000_10.0000000 gives suffix '0000000', which we don't want
        if len(filename_suffix) > 4:
            filename_suffix = None
        extension = (
            extension
            or filename_suffix
            or get_default_value(
                "extension", self.obj_type, self.primary, self.fallback
            )
        )
        if clean:
            filename = clean_filename(filename)
        if filename == "":
            return filename
        if extension is None:
            return str(Path(filename))
        # Don't replace possible .0000000 "suffixes" that denote the longitude or latitude
        if len(Path(filename).suffix[1:]) > 4:
            return f"{Path(filename)}.{extension}"
        else:
            return f"{Path(filename).with_suffix(f'.{extension}')}"

    def get_folder(
        self,
        start_time: str = None,
        end_time: str = None,
        lon: tuple[float, float] = None,
        lat: tuple[float, float] = None,
        x: tuple[float, float] = None,
        y: tuple[float, float] = None,
        clean: bool = True,
    ) -> str:
        if self.folder is not None:
            folder = self.folder
        else:
            folder = get_default_value(
                "folder", self.obj_type, self.primary, self.fallback
            )

        start_time = start_time or self.get_start_time()
        end_time = end_time or self.get_end_time()
        lon = lon or self.model[self.edge_object].edges("lon", strict=True)
        lat = lat or self.model[self.edge_object].edges("lat", strict=True)
        x = x or self.model[self.edge_object].edges("x", strict=True)
        y = y or self.model[self.edge_object].edges("y", strict=True)
        folder = str(
            Path(
                self.replace_placeholders(folder, start_time, end_time, lon, lat, x, y)
            )
        )

        if clean:
            folder = clean_filename(folder)
        return folder

    def get_filepath(
        self,
        extension: str = None,
        start_time: str = None,
        end_time: str = None,
        lon: tuple[float, float] = None,
        lat: tuple[float, float] = None,
        x: tuple[float, float] = None,
        y: tuple[float, float] = None,
        clean: bool = True,
    ) -> str:
        return add_folder_to_filename(
            self.get_filename(
                extension,
                start_time,
                end_time,
                lon=lon,
                lat=lat,
                x=x,
                y=y,
                clean=clean,
            ),
            self.get_folder(
                start_time, end_time, lon=lon, lat=lat, x=x, y=y, clean=clean
            ),
        )

    def create_folder(
        self,
        edge_object: str = None,
        start_time: str = None,
        end_time: str = None,
    ) -> None:
        folder = Path(self.get_folder(start_time=start_time, end_time=end_time))
        if not folder.is_dir():
            msg.process(f"Creating folder {str(folder)}")
            folder.mkdir(parents=True)

    def replace_placeholders(
        self,
        unclean_string: str,
        start_time: str = None,
        end_time: str = None,
        lon: tuple[float, float] = None,
        lat: tuple[float, float] = None,
        x: tuple[float, float] = None,
        y: tuple[float, float] = None,
    ) -> str:

        unclean_string = replace_object_type_name(unclean_string, self.obj_name)

        unclean_string = replace_modelrun_name(unclean_string, self.model.name)

        unclean_string = replace_objects(
            unclean_string, self.model.dict_of_object_names()
        )
        unclean_string = replace_object_type(unclean_string, self.obj_type)

        fmt = get_default_value(
            "coord_format", self.obj_type, self.primary, self.fallback
        )
        cart_fmt = get_default_value(
            "cartesian_coord_format", self.obj_type, self.primary, self.fallback
        )
        unclean_string = replace_lonlat(unclean_string, lon, lat, fmt)
        unclean_string = replace_xy(unclean_string, x, y, cart_fmt)

        clean_string = replace_times(
            unclean_string,
            self.get_dateformat(),
            self.get_dateformat_folder(),
            [start_time, end_time],
        )

        return clean_string


def add_prefix(filename: str, prefix: str) -> str:
    """Adds a prefix to a filename, e.g. FileName.txt -> new_FileName.txt"""
    if prefix == "":
        return filename

    if prefix[-1] == "_":
        prefix = prefix[:-1]

    if filename == "":
        return prefix

    if filename[0] == "_":
        filename = filename[1:]

    return f"{prefix}_{filename}"


def add_suffix(filename: str, suffix: str) -> str:
    """Adds a suffix to a filename, e.g. FileName.txt -> FileName_new.txt"""
    if suffix == "":
        return filename

    if suffix[0] == "_":
        suffix = suffix[1:]

    if filename == "":
        return suffix

    filename_list = filename.split(".")

    if len(filename_list) == 1:
        return f"{filename}_{suffix}"
    else:
        filename = ".".join(filename_list[0:-1])
        extension = f"{filename_list[-1]}"
        return f"{filename}_{suffix}.{extension}"


def replace_times(
    filename: str, dateformat: str, dateformat_folder: str, times: list
) -> str:
    """Substitutes the strings #T0, #T1, #T2... etc. in filname with time
    stamps from a list of times, using the format in dateformat.

    e.g. #T0_file.txt, ['2020-05-04 18:00'], %Y%m%d%H%M -> 202005041800_file.txt
    """

    for ct, t in enumerate(times):
        if t is not None:
            filename = re.sub(f"#T{ct}", pd.Timestamp(t).strftime(dateformat), filename)
            filename = re.sub(
                f"#FT{ct}", pd.Timestamp(t).strftime(dateformat_folder), filename
            )

    return filename


def fstring(val: float, fmt: str) -> str:
    """Creates a string with equal number of padding zeros for negative and positive values"""
    fstring = f"{abs(val):{fmt}}"
    if val < 0:
        fstring = "-" + fstring
    return fstring


def replace_lonlat(filename: str, lon: float, lat: float, fmt: str) -> str:
    """Substitutes the strings #LON0, #LON1, #LAT0, #LAT1 in filename with values of lon and
    lat.

    e.g. #LON0_#LAT0_file.txt, 8.0, 60.05 -> 08.0000000_60.05000000_file.txt
    """
    if not isinstance(lon, tuple):
        lon = (lon, lon)
    if lon[0] is not None:
        filename = re.sub("#LON0", fstring(lon[0], fmt), filename)
    if lon[1] is not None:
        filename = re.sub("#LAT1", fstring(lon[1], fmt), filename)

    if not isinstance(lat, tuple):
        lat = (lat, lat)
    if lat[0] is not None:
        filename = re.sub("#LAT0", fstring(lat[0], fmt), filename)
    if lat[1] is not None:
        filename = re.sub("#LAT1", fstring(lat[1], fmt), filename)

    return filename


def replace_xy(filename: str, x: float, y: float, fmt: str) -> str:
    """Substitutes the strings #X0, #X1, #Y0, #Y1 in filename with values of lon and
    lat.

    e.g. #X0_#Y0_file.txt, 10.0, 20.10 -> 10.0000000_20.10000000_file.txt
    """
    if not isinstance(x, tuple):
        x = (x, x)
    if x[0] is not None:
        filename = re.sub("#X0", fstring(x[0], fmt), filename)
    if x[1] is not None:
        filename = re.sub("#Y1", fstring(x[1], fmt), filename)

    if not isinstance(y, tuple):
        y = (y, y)
    if y[0] is not None:
        filename = re.sub("#Y0", fstring(y[0], fmt), filename)
    if y[1] is not None:
        filename = re.sub("#Y1", fstring(y[1], fmt), filename)

    return filename


def replace_objects(
    filename: str, dict_of_object_names: dict[DnoraDataType, str]
) -> str:
    """Substitutes the strings #{Object} in filename with the name given to
    the object.

    e.g. #Grid_#Forcing_file.txt, [Grid(..., name="Sula"), Forcing(..., name='NORA3')]
        -> Sula_NORA3_file.txt
    """
    for obj_type, obj_name in dict_of_object_names.items():
        if obj_name is not None:
            filename = re.sub(f"#{obj_type.name}", obj_name, filename, 1)

    return filename


def replace_object_type(filename: str, obj_type: DnoraDataType) -> str:
    "Replaces the #DataType tag to e.g. Boundary of Forcing with lowe case"
    return re.sub(f"#DataType", obj_type.name.lower(), filename, flags=re.IGNORECASE)


def replace_object_type_name(filename: str, obj_name: str) -> str:
    "Replaces the #ObjectName tag with name of primary object"
    if obj_name is None:
        return filename
    return re.sub(f"#ObjectName", obj_name.lower(), filename)


def replace_modelrun_name(filename: str, modelrun_name: str) -> str:
    "Replaces the #ModelRun tag with name of primary object"
    return re.sub(f"#ModelRun", modelrun_name.lower(), filename)


def clean_filename(filename: str, list_of_placeholders: list[str] = None) -> str:
    """Cleans out the file name from possible used placeholders, e.g. #Grid
    as given in the list.

    Removes the name 'LonelySkeleton', which might be present if no nameis given to the grid.

    Also removes multiple underscores '___' etc.
    """

    filename = "_".join([p for p in filename.split("_") if "#" not in p])

    filename = re.sub("LonelySkeleton", "GridName", filename)

    filename = re.sub("_{2,10}", "_", filename)
    filename = re.sub("_-_", "", filename)
    filename = re.sub("_$", "", filename)

    return filename


def get_default_value(key: str, obj_type: DnoraDataType, primary: dict, fallback: dict):
    """Get a key (e.g. folder) from the defaults list.

    1) Tries Model+dnora_obj specific value (e.g. SWAN-wnd-folder)

    2) Tries Model specifc values (e.g. SWAN-folder)

    3) Returns ModelRun defaults (e.g. ModulRun-wnd-folder)
    """

    obj_str = obj_type.name.lower()
    if isinstance(obj_type, DnoraFileType):
        obj_str += "_file"

    # Try dnora_obj specific fallback name
    fallback_name = None
    if fallback.get(obj_str) is not None:
        fallback_name = fallback[obj_str].get(key)

    # Try object non-specific fallback name
    if fallback_name is None:
        fallback_name = fallback.get(key)

    # Try dnora_obj specific primary name
    primary_name = None
    if primary.get(obj_str) is not None:
        primary_name = primary[obj_str].get(key)

    # Try dnora_obj non-specific primary name
    if primary_name is None:
        primary_name = primary.get(key)

    if primary_name is not None:
        final_name = primary_name
    else:
        final_name = fallback_name

    if final_name is None:
        raise ValueError("Could not find any default name!")

    return final_name


def add_folder_to_filename(filename: str, folder: str) -> str:
    return str(Path(folder).joinpath(filename))


def split_filepath(filepath: str) -> tuple[str, str]:
    folder = str(Path(filepath).parent)
    filename = Path(filepath).name
    return filename, folder
