import unittest
import sys
import numpy as np
import os
import yaml
from dnora import file_module

with open('data/defaults.yml', 'r') as file:
  defaults = yaml.safe_load(file)

class GetDefaultValues(unittest.TestCase):
	def test_filename(self):
		filename = file_module.get_default_value(key='filename', dnora_obj='boundary', primary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'spec#Boundary#Grid#T0#T1')

		filename = file_module.get_default_value(key='filename', dnora_obj='spectra', primary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'spec1d_#Spectra_#Grid_#T0-#T1')

	def test_folder(self):
		filename = file_module.get_default_value(key='folder', dnora_obj='grid', primary=defaults['REEF3D'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'output')

		filename = file_module.get_default_value(key='folder', dnora_obj='input_file', primary=defaults['REEF3D'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'MyHOS_oceanFolder')

	def test_dateformat(self):
		filename = file_module.get_default_value(key='dateformat', dnora_obj='forcing', primary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%d')

		filename = file_module.get_default_value(key='dateformat', dnora_obj='forcing', primary=defaults['WW3'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%dT%H%M')

class AddFolder(unittest.TestCase):
    def test_basic(self):
        path = file_module.add_folder_to_filename(filename='filename', folder='folder')
        self.assertEqual(path, 'folder/filename')

    def test_empty_file(self):
        path = file_module.add_folder_to_filename(filename='', folder='folder')
        self.assertEqual(path, 'folder')

    def test_empty_folder(self):
        path = file_module.add_folder_to_filename(filename='filename', folder='')
        self.assertEqual(path, 'filename')

    def test_with_extension(self):
        path = file_module.add_folder_to_filename(filename='filename.asc', folder='folder')
        self.assertEqual(path, 'folder/filename.asc')

class SplitPath(unittest.TestCase):
    def test_basic(self):
        filename, folder = file_module.split_filepath(filepath='folder/filename')
        self.assertEqual(filename, 'filename')
        self.assertEqual(folder, 'folder')

    def test_empty_file(self):
        filename, folder = file_module.split_filepath(filepath='folder/')
        self.assertEqual(folder, '.')
        self.assertEqual(filename, 'folder')

    def test_empty_folder(self):
        filename, folder = file_module.split_filepath(filepath='/filename')
        self.assertEqual(folder, '/')
        self.assertEqual(filename, 'filename')

    def test_with_extension(self):
        filename, folder = file_module.split_filepath(filepath='folder/filename.asc')
        self.assertEqual(filename, 'filename.asc')
        self.assertEqual(folder, 'folder')

    def test_multiple_folders(self):
        filename, folder = file_module.split_filepath(filepath='folder1/folder2/filename.asc')
        self.assertEqual(filename, 'filename.asc')
        self.assertEqual(folder, 'folder1/folder2')

class Clean(unittest.TestCase):
    def test_trivial(self):
        filename = file_module.clean(filename='filename', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_empty(self):
        filename = file_module.clean(filename='', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, '')

    def test_grid(self):
        filename = file_module.clean(filename='filename_#Grid', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_boundary(self):
        filename = file_module.clean(filename='filename_#Boundary', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_spectra(self):
        filename = file_module.clean(filename='filename_#Spectra', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_forcing(self):
        filename = file_module.clean(filename='filename_#Forcing', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_time(self):
        filename = file_module.clean(filename='filename_#T0', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_#T1', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_#T0-#T1', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_#T0_#T1', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_lonlat(self):
        filename = file_module.clean(filename='filename_#Lon', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_#Lat', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_E#Lon', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_N#Lat', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_#Lon_#Lat', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

        filename = file_module.clean(filename='filename_E#Lon_N#Lat', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename')

    def test_combinations(self):
        filename = file_module.clean(filename='filename_#Lon_#Grid#BoundarySaveThis#Forcing__', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename_SaveThis')

        filename = file_module.clean(filename='filename_#Lon_#T0-#T1#Grid#BoundarySaveThis#Forcing__', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename_SaveThis')

        filename = file_module.clean(filename='filename_#Lon_#T0-#T1#Grid#BoundarySaveThisE#Lon#Forcing__', list_of_placeholders=defaults['list_of_placeholders'])
        self.assertEqual(filename, 'filename_SaveThis')

class ReplaceObjects(unittest.TestCase):
    def test_all(self):
        dict_of_object_names={'ModelRun': 'TestModel', 'Grid': 'Sulafjorden', 'Boundary': 'NORA3', 'Forcing': 'MEPS', 'Spectra': 'FromBoundary'}
        filename = file_module.replace_objects('file_#Gird_#ModelRun_#Grid_#Boundary_#Forcing_#Spectra', dict_of_object_names)
        self.assertEqual(filename, 'file_#Gird_TestModel_Sulafjorden_NORA3_MEPS_FromBoundary')

class ReplaceLonLat(unittest.TestCase):
    def test_all(self):
        filename = file_module.replace_lonlat('file_#Lon_#ModelRun_#Lat', lon=5.3, lat=60.2222)
        self.assertEqual(filename, 'file_05.3000000_#ModelRun_60.2222000')

class ReplaceTime(unittest.TestCase):
    def test_empty(self):
        filename = file_module.replace_times('file_#T0', '%Y', [])
        self.assertEqual(filename, 'file_#T0')

    def test_T0(self):
        filename = file_module.replace_times('file_#T1', '%Y', ['2020-06-05 05:00'])
        self.assertEqual(filename, 'file_#T1')

        filename = file_module.replace_times('file_#T0', '%Y', ['2020-06-05 05:00'])
        self.assertEqual(filename, 'file_2020')

        filename = file_module.replace_times('file_#T0_#Grid_#T1', '%Y-%m-%d', ['2020-06-05 05:00'])
        self.assertEqual(filename, 'file_2020-06-05_#Grid_#T1')

    def test_both(self):
        filename = file_module.replace_times('file_#T1', '%Y-%m-%d %H:%M', ['2020-06-05 05:01', '2020-07-05 05:01'])
        self.assertEqual(filename, 'file_2020-07-05 05:01')

        filename = file_module.replace_times('file_#T0', '%Y', ['2020-06-05 05:00', '2020-07-05 05:00'])
        self.assertEqual(filename, 'file_2020')

        filename = file_module.replace_times('file_#T0_#Grid_#T1', '%Y-%m-%d', ['2020-06-05 05:00', '2020-07-05 05:00'])
        self.assertEqual(filename, 'file_2020-06-05_#Grid_2020-07-05')

class AddPrefix(unittest.TestCase):
    def test_all(self):
        filename = file_module.add_prefix('file_#T0', 'ddd')
        self.assertEqual(filename, 'ddd_file_#T0')

        filename = file_module.add_prefix('', 'ddd')
        self.assertEqual(filename, 'ddd')

        filename = file_module.add_prefix('_asd', 'ddd')
        self.assertEqual(filename, 'ddd_asd')

        filename = file_module.add_prefix('asd', 'ddd_')
        self.assertEqual(filename, 'ddd_asd')

        filename = file_module.add_prefix('_asd', 'ddd_')
        self.assertEqual(filename, 'ddd_asd')

class AddSuffix(unittest.TestCase):
    def test_all(self):
        filename = file_module.add_suffix('file_#T0', 'ddd')
        self.assertEqual(filename, 'file_#T0_ddd')

        filename = file_module.add_suffix('', 'ddd')
        self.assertEqual(filename, 'ddd')

        filename = file_module.add_suffix('asd.txt', 'ddd')
        self.assertEqual(filename, 'asd_ddd.txt')

        filename = file_module.add_suffix('asd', '_ddd')
        self.assertEqual(filename, 'asd_ddd')

        filename = file_module.add_suffix('asd.txt', '_ddd')
        self.assertEqual(filename, 'asd_ddd.txt')


if __name__ == '__main__':
	unittest.main()
