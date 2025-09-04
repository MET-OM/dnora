import unittest
import sys
import numpy as np
import os
import yaml
from dnora import file_module
from pathlib import Path
from dnora.type_manager.dnora_types import DnoraDataType, DnoraFileType
from dnora.defaults import read_defaults

defaults = read_defaults("export_defaults.yml", from_module=True)

# with open("data/defaults.yml", "r") as file:
#    defaults = yaml.safe_load(file)


class GetDefaultValues(unittest.TestCase):
    def test_filename(self):
        filename = file_module.get_default_value(
            key="filename",
            obj_type=DnoraDataType.SPECTRA,
            primary=defaults["SWAN"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "spec#SPECTRA#GRID#T0_#T1")

        filename = file_module.get_default_value(
            key="filename",
            obj_type=DnoraDataType.SPECTRA1D,
            primary=defaults["SWAN"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "spec1d_#SPECTRA1D_#GRID_#T0_#T1")

    def test_folder(self):
        filename = file_module.get_default_value(
            key="folder",
            obj_type=DnoraDataType.GRID,
            primary=defaults["REEF3D"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "#GRID_REEF3D")

        filename = file_module.get_default_value(
            key="folder",
            obj_type=DnoraFileType.INPUT,
            primary=defaults["HOS_OCEAN"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "MyHOS_oceanFolder")

    def test_dateformat(self):
        filename = file_module.get_default_value(
            key="dateformat",
            obj_type=DnoraDataType.WIND,
            primary=defaults["SWAN"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "%Y%m%d")

        filename = file_module.get_default_value(
            key="dateformat",
            obj_type=DnoraDataType.WIND,
            primary=defaults["WW3"],
            fallback=defaults["MODELRUN"],
        )
        self.assertEqual(filename, "%Y%m%dT%H%M")


class AddFolder(unittest.TestCase):
    def test_basic(self):
        path = file_module.add_folder_to_filename(filename="filename", folder="folder")
        self.assertEqual(Path(path).as_posix(), "folder/filename")

    def test_empty_file(self):
        path = file_module.add_folder_to_filename(filename="", folder="folder")
        self.assertEqual(path, "folder")

    def test_empty_folder(self):
        path = file_module.add_folder_to_filename(filename="filename", folder="")
        self.assertEqual(path, "filename")

    def test_with_extension(self):
        path = file_module.add_folder_to_filename(
            filename="filename.asc", folder="folder"
        )
        self.assertEqual(Path(path).as_posix(), "folder/filename.asc")


class SplitPath(unittest.TestCase):
    def test_basic(self):
        filename, folder = file_module.split_filepath(filepath="folder/filename")
        self.assertEqual(filename, "filename")
        self.assertEqual(folder, "folder")

    def test_empty_file(self):
        filename, folder = file_module.split_filepath(filepath="folder/")
        self.assertEqual(folder, ".")
        self.assertEqual(filename, "folder")

    def test_empty_folder(self):
        filename, folder = file_module.split_filepath(filepath="/filename")
        self.assertEqual(Path(folder).as_posix(), "/")
        self.assertEqual(filename, "filename")

    def test_with_extension(self):
        filename, folder = file_module.split_filepath(filepath="folder/filename.asc")
        self.assertEqual(filename, "filename.asc")
        self.assertEqual(folder, "folder")

    def test_multiple_folders(self):
        filename, folder = file_module.split_filepath(
            filepath="folder1/folder2/filename.asc"
        )
        self.assertEqual(filename, "filename.asc")
        self.assertEqual(Path(folder).as_posix(), "folder1/folder2")


class Clean(unittest.TestCase):
    def test_trivial(self):
        filename = file_module.clean_filename(
            filename="filename", list_of_placeholders=defaults["list_of_placeholders"]
        )
        self.assertEqual(filename, "filename")

    def test_empty(self):
        filename = file_module.clean_filename(
            filename="", list_of_placeholders=defaults["list_of_placeholders"]
        )
        self.assertEqual(filename, "")

    def test_grid(self):
        filename = file_module.clean_filename(
            filename="filename_#GRID",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_boundary(self):
        filename = file_module.clean_filename(
            filename="filename_#SPECTRA",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_spectra(self):
        filename = file_module.clean_filename(
            filename="filename_#SPECTRA1D",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_forcing(self):
        filename = file_module.clean_filename(
            filename="filename_#WIND",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_time(self):
        filename = file_module.clean_filename(
            filename="filename_#T0",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_#T1",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_#T0_#T1",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_#T0_#T1",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_lonlat(self):
        filename = file_module.clean_filename(
            filename="filename_#LON",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_#LAT",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_E#LON",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_N#LAT",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_#LON_#LAT",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

        filename = file_module.clean_filename(
            filename="filename_E#LON_N#LAT",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename")

    def test_combinations(self):
        filename = file_module.clean_filename(
            filename="filename_#LON_#GRID#SPECTRA_SaveThis_#WIND__",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename_SaveThis")

        filename = file_module.clean_filename(
            filename="filename_#LON_#T0_#T1#GRID#SPECTRA_SaveThis_#WIND__",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename_SaveThis")

        filename = file_module.clean_filename(
            filename="filename_#LON_#T0_#T1#GRID_#SPECTRA_SaveThis_E#LON#WIND__",
            list_of_placeholders=defaults["list_of_placeholders"],
        )
        self.assertEqual(filename, "filename_SaveThis")


class ReplaceObjects(unittest.TestCase):
    def test_all(self):
        dict_of_object_names = {
            DnoraDataType.GRID: "Sulafjorden",
            DnoraDataType.SPECTRA: "NORA3",
            DnoraDataType.WIND: "MEPS",
            DnoraDataType.SPECTRA1D: "FromBoundary",
        }
        filename = file_module.replace_objects(
            "file_#GIRD_#GRID_#SPECTRA_#WIND_#SPECTRA1D",
            dict_of_object_names,
        )
        self.assertEqual(filename, "file_#GIRD_Sulafjorden_NORA3_MEPS_FromBoundary")


class ReplaceLonLat(unittest.TestCase):
    def test_all(self):
        filename = file_module.replace_lonlat(
            "file_#LON0_#ModelRun_#LAT0", lon=5.3, lat=60.2222, fmt="010.7f"
        )
        self.assertEqual(filename, "file_05.3000000_#ModelRun_60.2222000")


class ReplaceTime(unittest.TestCase):
    def test_empty(self):
        filename = file_module.replace_times(
            "file_#T0", dateformat="%Y", dateformat_folder="", times=[]
        )
        self.assertEqual(filename, "file_#T0")

    def test_T0(self):
        filename = file_module.replace_times(
            "file_#T1",
            dateformat="%Y",
            dateformat_folder="",
            times=["2020-06-05 05:00"],
        )
        self.assertEqual(filename, "file_#T1")

        filename = file_module.replace_times(
            "file_#T0",
            dateformat="%Y",
            dateformat_folder="",
            times=["2020-06-05 05:00"],
        )
        self.assertEqual(filename, "file_2020")

        filename = file_module.replace_times(
            "file_#T0_#GRID_#T1",
            dateformat="%Y-%m-%d",
            dateformat_folder="",
            times=["2020-06-05 05:00"],
        )
        self.assertEqual(filename, "file_2020-06-05_#GRID_#T1")

    def test_both(self):
        filename = file_module.replace_times(
            "file_#T1",
            dateformat="%Y-%m-%d %H:%M",
            dateformat_folder="",
            times=["2020-06-05 05:01", "2020-07-05 05:01"],
        )
        self.assertEqual(filename, "file_2020-07-05 05:01")

        filename = file_module.replace_times(
            "file_#T0",
            dateformat="%Y",
            dateformat_folder="",
            times=["2020-06-05 05:00", "2020-07-05 05:00"],
        )
        self.assertEqual(filename, "file_2020")

        filename = file_module.replace_times(
            "file_#T0_#GRID_#T1",
            dateformat="%Y-%m-%d",
            dateformat_folder="",
            times=["2020-06-05 05:00", "2020-07-05 05:00"],
        )
        self.assertEqual(filename, "file_2020-06-05_#GRID_2020-07-05")

    def test_folder(self):
        filename = file_module.replace_times(
            "#FT0/file_#T1",
            dateformat="%Y-%m-%d %H:%M",
            dateformat_folder="%Y",
            times=["2020-06-05 05:01", "2020-07-05 05:01"],
        )
        self.assertEqual(filename, "2020/file_2020-07-05 05:01")

        filename = file_module.replace_times(
            "#FT0/file_#T0",
            dateformat="%Y-%m-%d",
            dateformat_folder="%Y/%m",
            times=["2020-06-05 05:00", "2020-07-05 05:00"],
        )
        self.assertEqual(filename, "2020/06/file_2020-06-05")

        filename = file_module.replace_times(
            "#FT1/file_#T0_#GRID_#T1",
            dateformat="%Y-%m-%d",
            dateformat_folder="%Y-%m",
            times=["2020-06-05 05:00", "2020-07-05 05:00"],
        )
        self.assertEqual(filename, "2020-07/file_2020-06-05_#GRID_2020-07-05")


class AddPrefix(unittest.TestCase):
    def test_all(self):
        filename = file_module.add_prefix("file_#T0", "ddd")
        self.assertEqual(filename, "ddd_file_#T0")

        filename = file_module.add_prefix("", "ddd")
        self.assertEqual(filename, "ddd")

        filename = file_module.add_prefix("_asd", "ddd")
        self.assertEqual(filename, "ddd_asd")

        filename = file_module.add_prefix("asd", "ddd_")
        self.assertEqual(filename, "ddd_asd")

        filename = file_module.add_prefix("_asd", "ddd_")
        self.assertEqual(filename, "ddd_asd")


class AddSuffix(unittest.TestCase):
    def test_all(self):
        filename = file_module.add_suffix("file_#T0", "ddd")
        self.assertEqual(filename, "file_#T0_ddd")

        filename = file_module.add_suffix("", "ddd")
        self.assertEqual(filename, "ddd")

        filename = file_module.add_suffix("asd.txt", "ddd")
        self.assertEqual(filename, "asd_ddd.txt")

        filename = file_module.add_suffix("asd", "_ddd")
        self.assertEqual(filename, "asd_ddd")

        filename = file_module.add_suffix("asd.txt", "_ddd")
        self.assertEqual(filename, "asd_ddd.txt")


if __name__ == "__main__":
    unittest.main()
