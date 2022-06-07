import unittest
import sys
sys.path.insert(0, "../../")
from dnora.aux_funcs import create_time_stamps, u_v_from_dir
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

if __name__ == '__main__':
	unittest.main()
