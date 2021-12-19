import unittest
import sys
sys.path.insert(0, "../../")
from dnora import bnd
import numpy as np

class BoundarySpectralConventions(unittest.TestCase):
	def test_ocean_to_met(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='Met')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(18,36+18,1.), 36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_ocean(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='Ocean')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(18,36+18,1.), 36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ocean_to_ww3(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='WW3')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ocean_to_math(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='Math')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

if __name__ == '__main__':
	unittest.main()
