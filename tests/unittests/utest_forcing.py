import unittest
import sys
sys.path.insert(0, "../../")
from dnora.wnd import Forcing
from dnora.grd import Grid
import numpy as np

class ForcingMethodsInit(unittest.TestCase):
	def test_init(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		forcing = Forcing(grid, name='TestForcing')

		self.assertEqual(forcing.nt(), 0)
		self.assertEqual(forcing.ny(), 0)
		self.assertEqual(forcing.ny(), 0)
		self.assertEqual(forcing.name(), 'TestForcing')
		self.assertIsNone(forcing.time())
		self.assertEqual(forcing.size(), (1,1,0))
		self.assertIsNone(np.testing.assert_almost_equal(np.array([[[]]]),forcing.u()))
		self.assertIsNone(np.testing.assert_almost_equal(np.array([[[]]]),forcing.v()))

if __name__ == '__main__':
	unittest.main()
