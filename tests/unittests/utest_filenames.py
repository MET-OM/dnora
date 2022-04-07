import unittest
import sys
import numpy as np
import os
import yaml
from dnora.aux import get_default_value

with open('data/defaults.yml', 'r') as file:
  defaults = yaml.safe_load(file)

class FileNames(unittest.TestCase):
	def test_filename(self):
		filename = get_default_value(key='filename', module='bnd', primary='primary_name', secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'primary_name')
	
		filename = get_default_value(key='filename', module='bnd', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'spec#Boundary#Grid#T0#T1')

		filename = get_default_value(key='filename', module='spc', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'spec_#Boundary_#Grid_#T0-#T1')

	def test_folder(self):
		filename = get_default_value(key='folder', module='wnd', primary='primary_folder', secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'primary_folder')

		filename = get_default_value(key='folder', module='wnd', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'wndoutput')

		filename = get_default_value(key='folder', module='bnd', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'swanoutput')

		filename = get_default_value(key='folder', module='spc', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'swanoutput')

		filename = get_default_value(key='folder', module='wnd', primary=None, secondary=defaults['WW3'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'wndoutput')

		filename = get_default_value(key='folder', module='bnd', primary=None, secondary=defaults['WW3'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'output')

		filename = get_default_value(key='folder', module='bnd', primary=None, secondary=defaults['SWASH'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'output')

		filename = get_default_value(key='folder', module='spc', primary=None, secondary=defaults['SWASH'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'output')

	def test_dateformat(self):
		filename = get_default_value(key='dateformat', module='wnd', primary='primary_dateformat', secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, 'primary_dateformat')

		filename = get_default_value(key='dateformat', module='wnd', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%d')

		filename = get_default_value(key='dateformat', module='spc', primary=None, secondary=defaults['SWAN'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%d')

		filename = get_default_value(key='dateformat', module='grd', primary=None, secondary=defaults['SWASH'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%d')

		filename = get_default_value(key='dateformat', module='wnd', primary=None, secondary=defaults['SWASH'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%dT%H%M')

		filename = get_default_value(key='dateformat', module='spc', primary=None, secondary=defaults['SWASH'], fallback=defaults['ModelRun'])
		self.assertEqual(filename, '%Y%m%dT%H%M')

if __name__ == '__main__':
	unittest.main()
