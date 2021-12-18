import unittest
import sys
sys.path.insert(0, "../../")
from dnora.aux import create_time_stamps, u_v_from_dir
from dnora.grd import Grid
import numpy as np

class TimeStamps(unittest.TestCase):
	def test_stride(self):
		"""Test that stride works"""
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 06:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12)

		self.assertEqual(len(st),3)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')
		self.assertEqual(st[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 23:00')
		self.assertEqual(et[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 06:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')
		self.assertEqual(ft[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

	def test_lead_time(self):
		"""Test that leadtime works"""
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 00:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12, lead_time=24)

		self.assertEqual(len(st),3)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')
		self.assertEqual(st[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 23:00')
		self.assertEqual(et[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2019-12-31 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2019-12-31 12:00')
		self.assertEqual(ft[2].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')

	def test_negative_lead_time(self):
		"""Test that negative lead_time works"""
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 00:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12, lead_time=-24)

		self.assertEqual(len(st),3)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')
		self.assertEqual(st[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 23:00')
		self.assertEqual(et[2].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-02 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2020-01-02 12:00')
		self.assertEqual(ft[2].strftime('%Y-%m-%d %H:%M'), '2020-01-03 00:00')

	def test_last_file(self):
		"""Test that last_file works"""
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-02 00:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12, last_file='2020-01-01 12:00')

		self.assertEqual(len(st),2)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 23:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

	def test_hours_per_file(self):
		"""Test that last_file works"""
		start_time = '2020-01-01 00:00'
		end_time = '2020-01-10 00:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12, last_file='2020-01-01 12:00', hours_per_file=48)

		self.assertEqual(len(st),2)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-03 11:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

	def test_between_files(self):
		"""Test that it works to start and stop between files"""
		start_time = '2020-01-01 02:00'
		end_time = '2020-01-01 17:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12)

		self.assertEqual(len(st),2)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 02:00')
		self.assertEqual(st[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 17:00')

		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')
		self.assertEqual(ft[1].strftime('%Y-%m-%d %H:%M'), '2020-01-01 12:00')

	def test_start_after_last_file(self):
		"""Test that it works to start and stop between files"""
		start_time = '2020-01-01 11:00'
		end_time = '2020-01-05 17:00'

		st, et, ft = create_time_stamps(start_time, end_time, stride=12, last_file='2020-01-01 00:00', hours_per_file=72)

		self.assertEqual(len(st),1)

		self.assertEqual(st[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 11:00')
		self.assertEqual(et[0].strftime('%Y-%m-%d %H:%M'), '2020-01-03 23:00')
		self.assertEqual(ft[0].strftime('%Y-%m-%d %H:%M'), '2020-01-01 00:00')

class WindDir(unittest.TestCase):
	def test_main_directions(self):
		ws = 1
		u,v = u_v_from_dir(ws,wdir=0)
		self.assertAlmostEqual(u,0)
		self.assertAlmostEqual(v,-ws)

		u,v = u_v_from_dir(ws,wdir=90)
		self.assertAlmostEqual(u,-ws)
		self.assertAlmostEqual(v,0)

		ws = 2
		u,v = u_v_from_dir(ws,wdir=180)
		self.assertAlmostEqual(u,0)
		self.assertAlmostEqual(v,ws)

		u,v = u_v_from_dir(ws,wdir=270)
		self.assertAlmostEqual(u,ws)
		self.assertAlmostEqual(v,0)


	def test_intercardinal_directions(self):
		wss = [16, 4, 9, 25]
		wds = [45, 135, 225, 315]

		for n in range(len(wss)):
			ws = wss[n]
			wd = wds[n]
			u,v = u_v_from_dir(ws,wd)
			self.assertAlmostEqual(u,-np.sin(wd*np.pi/180)*ws)
			self.assertAlmostEqual(v,-np.cos(wd*np.pi/180)*ws)


	def test_rabdom(self):
		wss = [17.21051191389413, 28.14738017880853, 4.5045008368608075, 8.214697597599898, 18.236349202702314, 10.332676364125264, 16.18891903347092, 6.821218174835286, 3.6122019719718423, 27.07374666146082, 5.164601452349762, 22.598053105326606, 28.460914144820652, 29.91596032156754, 18.499791087402073]
		wds = [28.45239639403414, 331.601264970793, 290.6987083729437, 234.89049223957838, 351.3836667740746, 312.44884615774475, 81.00506442163805, 184.71490864230987, 115.3222319980693, 85.4316085482833, 66.41164501632166, 93.19206404543112, 12.30946847939622, 143.29788976308757, 255.9299740907829]


		for n in range(len(wss)):
			ws = wss[n]
			wd = wds[n]
			u,v = u_v_from_dir(ws,wd)
			self.assertAlmostEqual(u,-np.sin(wd*np.pi/180)*ws)
			self.assertAlmostEqual(v,-np.cos(wd*np.pi/180)*ws)

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
