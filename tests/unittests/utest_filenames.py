import unittest
import sys
sys.path.insert(0, "../../")
from dnora.grd import Grid
from dnora.wnd import Forcing
from dnora.bnd import Boundary
from dnora.mdl import ModelRun
from dnora.aux import clean_filename
import numpy as np
import os

class FileNames(unittest.TestCase):
	def test_name_grid(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)
		filestring = 'test_#Grid_#Forcing_#Boundary_#T0_#T1_test2'
		datestring = '%Y%m%dT%H%MZ'

		filename = model.filename(filestring=filestring, datestring=datestring, extension='test')
		exp_filename = 'test_TestGrid_#Forcing_#Boundary_20200101T0000Z_20200102T0600Z_test2.test'

		self.assertEqual(filename, exp_filename)

		exp_clean_filename = 'test_TestGrid_20200101T0000Z_20200102T0600Z_test2.test'
		self.assertEqual(clean_filename(filename, ['#Forcing', '#Boundary']), exp_clean_filename)

	def test_name_grid_forcing(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)
		model._forcing = Forcing(grid, name='TestForcing')
		filestring = 'test_#Grid_#Forcing_#Boundary_#T0_#T1_test2'
		datestring = '%Y%m%dT%H%MZ'

		filename = model.filename(filestring=filestring, datestring=datestring, extension='test')
		exp_filename = 'test_TestGrid_TestForcing_#Boundary_20200101T0000Z_20200102T0600Z_test2.test'

		self.assertEqual(filename, exp_filename)

		exp_clean_filename = 'test_TestGrid_TestForcing_20200101T0000Z_20200102T0600Z_test2.test'
		self.assertEqual(clean_filename(filename, ['#Forcing', '#Boundary']), exp_clean_filename)

	def test_name_grid_forcing_boundary(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)
		model._forcing = Forcing(grid, name='TestForcing')
		model._boundary = Boundary(grid, name='TestBoundary')

		filestring = 'test_#Grid_#Forcing_#Boundary_#T0_#T1_test2'
		datestring = '%Y%m%dT%H%MZ'

		filename = model.filename(filestring=filestring, datestring=datestring, extension='test')
		exp_filename = 'test_TestGrid_TestForcing_TestBoundary_20200101T0000Z_20200102T0600Z_test2.test'

		self.assertEqual(filename, exp_filename)

		exp_clean_filename = 'test_TestGrid_TestForcing_TestBoundary_20200101T0000Z_20200102T0600Z_test2.test'
		self.assertEqual(clean_filename(filename, ['#Forcing', '#Boundary']), exp_clean_filename)

	def test_grid_exported(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		# Should produce default name
		filename = model.grid_exported_as(out_format='SWAN')
		exp_filename = 'TestGrid_SWAN.asc'
		self.assertEqual(filename, exp_filename)

		folder = model.grid_exported_to(out_format='SWAN')
		exp_folder = f"{os.getcwd()}/output"
		self.assertEqual(folder, exp_folder)

		path = model.grid_exported_path(out_format='SWAN')
		exp_path = f"{exp_folder}/{exp_filename}"
		self.assertEqual(path, exp_path)

	def test_empty_forcing_exported(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		# No forcing exists, so should give empyt name
		filename = model.forcing_exported_as(out_format='SWAN')
		exp_filename = ''
		self.assertEqual(filename, exp_filename)

		folder = model.forcing_exported_to(out_format='SWAN')
		exp_folder = ''
		self.assertEqual(folder, exp_folder)


		path = model.forcing_exported_path(out_format='SWAN')
		exp_path = ''
		self.assertEqual(path, exp_path)

	def test_forcing_exported(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		model._forcing = Forcing(grid, name='TestForcing')

		filename = model.forcing_exported_as(out_format='SWAN')
		exp_filename = 'windTestForcingTestGrid20200101_20200102.asc'
		self.assertEqual(filename, exp_filename)

		folder = model.forcing_exported_to(out_format='SWAN')
		exp_folder = f"{os.getcwd()}/output"
		self.assertEqual(folder, exp_folder)


		path = model.forcing_exported_path(out_format='SWAN')
		exp_path = f"{exp_folder}/{exp_filename}"
		self.assertEqual(path, exp_path)

	def test_empty_boundary_exported(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		# No boundary exists, so should give empyt name
		filename = model.boundary_exported_as(out_format='SWAN')
		exp_filename = ''
		self.assertEqual(filename, exp_filename)

		folder = model.boundary_exported_to(out_format='SWAN')
		exp_folder = ''
		self.assertEqual(folder, exp_folder)


		path = model.boundary_exported_path(out_format='SWAN')
		exp_path = ''
		self.assertEqual(path, exp_path)

	def test_boundary_exported(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		model._boundary = Boundary(grid, name='TestBoundary')

		filename = model.boundary_exported_as(out_format='SWAN')
		exp_filename = 'specTestBoundaryTestGrid20200101_20200102.asc'
		self.assertEqual(filename, exp_filename)

		folder = model.boundary_exported_to(out_format='SWAN')
		exp_folder = f"{os.getcwd()}/output"
		self.assertEqual(folder, exp_folder)


		path = model.boundary_exported_path(out_format='SWAN')
		exp_path = f"{exp_folder}/{exp_filename}"
		self.assertEqual(path, exp_path)

	def test_input_file_written(self):

		grid = Grid(lon=(0,1), lat=(0,1), name='TestGrid')
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		model = ModelRun(grid, start_time, end_time)

		filename = model.input_file_written_as(out_format='SWAN')
		exp_filename = 'input_20200101_TestGrid.swn'
		self.assertEqual(filename, exp_filename)

		folder = model.input_file_written_to(out_format='SWAN')
		exp_folder = 'MySWANFolder'
		self.assertEqual(folder, exp_folder)


		path = model.input_file_written_path(out_format='SWAN')
		exp_path = f"{exp_folder}/{exp_filename}"
		self.assertEqual(path, exp_path)

if __name__ == '__main__':
	unittest.main()
