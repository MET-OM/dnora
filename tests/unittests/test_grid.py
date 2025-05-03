import unittest
import sys
sys.path.insert(0, "../../")
from dnora.grd import Grid
import numpy as np

class GridMethodsInit(unittest.TestCase):
	def test_init(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		self.assertEqual(np.array_equal(grid.lon(),np.array([1,2])), True)
		self.assertEqual(np.array_equal(grid.lat(),np.array([0,2])), True)
		self.assertEqual(grid.nx(), 2)
		self.assertEqual(grid.ny(), 2)
		self.assertEqual(grid.size(), (2,2))

	def test_init_one_point_in_lat(self):
		grid = Grid(lon=(1,2), lat=(0,0))
		self.assertEqual(np.array_equal(grid.lon(),np.array([1, 2])), True)
		self.assertEqual(np.array_equal(grid.lat(),np.array([0])), True)
		self.assertEqual(grid.nx(), 2)
		self.assertEqual(grid.ny(), 1)
		self.assertEqual(grid.size(), (1,2))

	def test_init_one_point_in_lon(self):
		grid = Grid(lon=(1,1), lat=(0,2))
		self.assertEqual(np.array_equal(grid.lon(),np.array([1])), True)
		self.assertEqual(np.array_equal(grid.lat(),np.array([0,2])), True)
		self.assertEqual(grid.nx(), 1)
		self.assertEqual(grid.ny(), 2)
		self.assertEqual(grid.size(), (2,1))

	def test_init_one_point(self):
		grid = Grid(lon=(1,1), lat=(0,0))
		self.assertEqual(np.array_equal(grid.lon(),np.array([1])), True)
		self.assertEqual(np.array_equal(grid.lat(),np.array([0])), True)
		self.assertEqual(grid.nx(), 1)
		self.assertEqual(grid.ny(), 1)
		self.assertEqual(grid.size(), (1,1))


class GridMethodsNxSpacing(unittest.TestCase):
	def test_nx_spacing_y_trivial(self):
		grid = Grid(lon=(1,2), lat=(0,0))
		grid.set_spacing(nx=11, ny=11)
		self.assertEqual(grid.nx(), 11)
		self.assertEqual(grid.ny(), 1)
		self.assertEqual(grid.size(), (1,11))
		self.assertAlmostEqual(grid.dlon(), 0.1)
		self.assertAlmostEqual(grid.dlat(), 0.)
		self.assertAlmostEqual(grid.dx(), 10108.629694959884)
		self.assertAlmostEqual(grid.dy(), 0.)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.001, 1/10),grid.lon()))
		self.assertAlmostEqual(np.array([0.]), grid.lat())

	def test_nx_spacing_x_trivial(self):
		grid = Grid(lon=(1,1), lat=(0,2))
		grid.set_spacing(nx=11, ny=11)
		self.assertEqual(grid.nx(), 1)
		self.assertEqual(grid.ny(), 11)
		self.assertEqual(grid.size(), (11,1))
		self.assertAlmostEqual(grid.dlon(), 0.)
		self.assertAlmostEqual(grid.dlat(), 0.2)
		self.assertAlmostEqual(grid.dx(), 0.)
		self.assertAlmostEqual(grid.dy(), 20217.259389919767)

		self.assertAlmostEqual(np.array([1.]), grid.lon())
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.001, 2/10),grid.lat()))

	def test_nx_spacing(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		grid.set_spacing(nx=11, ny=11)
		self.assertEqual(grid.nx(), 11)
		self.assertEqual(grid.ny(), 11)
		self.assertEqual(grid.size(), (11,11))
		self.assertAlmostEqual(grid.dlon(), 0.1)
		self.assertAlmostEqual(grid.dlat(), 0.2)
		self.assertAlmostEqual(grid.dx(), 10107.09006262059)
		self.assertAlmostEqual(grid.dy(), 20217.259389919767)

	def test_dm_spacing_y_trivial(self):
		grid = Grid(lon=(1,2), lat=(0,0))
		grid.set_spacing(dm=1000)
		self.assertEqual(grid.nx(), 112)
		self.assertEqual(grid.ny(), 1)
		self.assertEqual(grid.size(), (1,112))
		self.assertAlmostEqual(grid.dlon(), 1/111)
		self.assertAlmostEqual(grid.dlat(), 0.)
		self.assertAlmostEqual(grid.dx(), 992.8118450407029)
		self.assertAlmostEqual(grid.dy(), 0.)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.001, 1/111),grid.lon()))
		self.assertAlmostEqual(np.array([0.]), grid.lat())

	def test_dm_spacing_x_trivial(self):
		grid = Grid(lon=(1,1), lat=(0,2))
		grid.set_spacing(dm=1000)
		self.assertEqual(grid.nx(), 1)
		self.assertEqual(grid.ny(), 223)
		self.assertEqual(grid.size(), (223,1))
		self.assertAlmostEqual(grid.dlon(), 0.)
		self.assertAlmostEqual(grid.dlat(), 2/222)
		self.assertAlmostEqual(grid.dx(), 0.)
		self.assertAlmostEqual(grid.dy(), 997.2639160946972)

		self.assertAlmostEqual(np.array([1.]), grid.lon())
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.001, 2/222),grid.lat()))

	def test_dm_spacing(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		grid.set_spacing(dm=1000)
		self.assertEqual(grid.nx(), 112)
		self.assertEqual(grid.ny(), 223)
		self.assertEqual(grid.size(), (223,112))
		self.assertAlmostEqual(grid.dlon(), 1/111)
		self.assertAlmostEqual(grid.dlat(), 2/222)
		self.assertAlmostEqual(grid.dx(), 992.6606311502364)
		self.assertAlmostEqual(grid.dy(), 997.2639160946972)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.001, 1/111),grid.lon()))
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.001, 2/222),grid.lat()))

class GridMethodsDlonSpacing(unittest.TestCase):
	def test_dlon_spacing_y_trivial(self):
		grid = Grid(lon=(1,2), lat=(0,0))
		grid.set_spacing(dlon=0.1, dlat=0.1)
		self.assertEqual(grid.nx(), 11)
		self.assertEqual(grid.ny(), 1)
		self.assertEqual(grid.size(), (1,11))
		self.assertAlmostEqual(grid.dlon(), 0.1)
		self.assertAlmostEqual(grid.dlat(), 0.)
		self.assertAlmostEqual(grid.dx(), 10108.629694959884)
		self.assertAlmostEqual(grid.dy(), 0.)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.01, 0.1),grid.lon()))
		self.assertAlmostEqual(np.array([0.]), grid.lat())


	def test_dlon_spacing_x_trivial(self):
		grid = Grid(lon=(1,1), lat=(0,2))
		grid.set_spacing(dlon=0.1, dlat=0.1)
		self.assertEqual(grid.nx(), 1)
		self.assertEqual(grid.ny(), 21)
		self.assertEqual(grid.size(), (21,1))
		self.assertAlmostEqual(grid.dlon(), 0.)
		self.assertAlmostEqual(grid.dlat(), 0.1)
		self.assertAlmostEqual(grid.dx(), 0.)
		self.assertAlmostEqual(grid.dy(), 10589.993013767498)

		self.assertAlmostEqual(np.array([1.]), grid.lon())
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.01, 0.1),grid.lat()))

	def test_dlon_spacing(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		grid.set_spacing(dlon=0.1, dlat=0.1)
		self.assertEqual(grid.nx(), 11)
		self.assertEqual(grid.ny(), 21)
		self.assertEqual(grid.size(), (21,11))
		self.assertAlmostEqual(grid.dlon(), 0.1)
		self.assertAlmostEqual(grid.dlat(), 0.1)
		self.assertAlmostEqual(grid.dx(), 10107.09006262059)
		self.assertAlmostEqual(grid.dy(), 10589.993013767498)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.01, 0.1),grid.lon()))
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.01, 0.1),grid.lat()))

	def test_dlon_spacing_rounding(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		grid.set_spacing(dlon=0.11, dlat=0.09)
		self.assertEqual(grid.nx(), 10)
		self.assertEqual(grid.ny(), 23)
		self.assertEqual(grid.size(), (23,10))
		self.assertAlmostEqual(grid.dlon(), 0.1111111111111111)
		self.assertAlmostEqual(grid.dlat(), 0.09090909090909091)
		self.assertAlmostEqual(grid.dx(), 11117.799068882648)
		self.assertAlmostEqual(grid.dy(), 9669.124056048586)

		self.assertAlmostEqual(grid.lon()[0], 1.)
		self.assertAlmostEqual(grid.lon()[-1], 2.)
		self.assertAlmostEqual(grid.lat()[0], 0.)
		self.assertAlmostEqual(grid.lat()[-1], 2.)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.01, 0.1111111111111111),grid.lon()))
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.01, 0.09090909090909091),grid.lat()))

	def test_dlon_floating_edge(self):
		grid = Grid(lon=(1,2), lat=(0,2))
		grid.set_spacing(dlon=0.11, dlat=0.09, floating_edge=True)
		self.assertEqual(grid.nx(), 10)
		self.assertEqual(grid.ny(), 23)
		self.assertEqual(grid.size(), (23,10))
		self.assertAlmostEqual(grid.dlon(), 0.11)
		self.assertAlmostEqual(grid.dlat(), 0.09)
		self.assertAlmostEqual(grid.dx(), 11006.65444371987)
		self.assertAlmostEqual(grid.dy(), 9572.432815488102)

		self.assertAlmostEqual(grid.lon()[0], 1.)
		self.assertAlmostEqual(grid.lon()[-1], 1.99)
		self.assertAlmostEqual(grid.lat()[0], 0.)
		self.assertAlmostEqual(grid.lat()[-1], 1.98)

		self.assertIsNone(np.testing.assert_almost_equal(np.arange(1, 2.0, 0.11),grid.lon()))
		self.assertIsNone(np.testing.assert_almost_equal(np.arange(0, 2.0, 0.09),grid.lat()))

if __name__ == '__main__':
	unittest.main()
