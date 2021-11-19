import unittest
import sys
sys.path.insert(0, "../../")
from dnora.aux import create_time_stamps, u_v_from_dir
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
		u,v = u_v_from_dir(ws=1,wdir=0)
		self.assertAlmostEqual(u,0)
		self.assertAlmostEqual(v,-1)

		u,v = u_v_from_dir(ws=1,wdir=90)
		self.assertAlmostEqual(u,-1)
		self.assertAlmostEqual(v,0)

		u,v = u_v_from_dir(ws=2,wdir=180)
		self.assertAlmostEqual(u,0)
		self.assertAlmostEqual(v,2)

		u,v = u_v_from_dir(ws=2,wdir=270)
		self.assertAlmostEqual(u,2)
		self.assertAlmostEqual(v,0)


	def test_intercardinal_directions(self):
		u,v = u_v_from_dir(ws=16,wdir=45)
		self.assertAlmostEqual(u,-np.cos(45*np.pi/180)*16)
		self.assertAlmostEqual(v,-np.sin(45*np.pi/180)*16)

		u,v = u_v_from_dir(ws=4,wdir=135)
		self.assertAlmostEqual(u,np.cos(135*np.pi/180)*4)
		self.assertAlmostEqual(v,np.sin(135*np.pi/180)*4)

		u,v = u_v_from_dir(ws=9,wdir=225)
		self.assertAlmostEqual(u,-np.cos(225*np.pi/180)*9)
		self.assertAlmostEqual(v,-np.sin(225*np.pi/180)*9)

		u,v = u_v_from_dir(ws=25,wdir=315)
		self.assertAlmostEqual(u,np.cos(315*np.pi/180)*25)
		self.assertAlmostEqual(v,np.sin(315*np.pi/180)*25)


if __name__ == '__main__':
	unittest.main()
		
