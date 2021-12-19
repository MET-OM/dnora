import unittest
import sys
sys.path.insert(0, "../../")
from dnora import bnd
import numpy as np
from copy import copy

class BoundarySpectralConventionsOcean(unittest.TestCase):
	def test_ocean_to_met(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='Met')
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

	def test_ocean_to_mathvec(self):
		D=np.arange(0,360,10.)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='MathVec')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

class BoundarySpectralConventionsWW3(unittest.TestCase):
	def test_ww3_to_ocean(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Ocean')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,1.), np.arange(0,36,1.)])))

	def test_ww3_to_met(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Met')
		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		stemp = np.mod(np.arange(18,36+18,1.), 36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))


	def test_ww3_to_math(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Math')
		if not isinstance(bnd_processor, list):
			bnd_processor = [bnd_processor]

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ww3_to_mathvec(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='MathVec')
		if not isinstance(bnd_processor, list):
			bnd_processor = [bnd_processor]

		Snew = copy(S)
		Dnew = copy(D)

		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,1.), np.arange(0,36,1.)])))

class BoundarySpectralConventionsMet(unittest.TestCase):

	def test_met_to_ocean(self):
		D=np.arange(0,360,10.)
		S=np.array([ np.mod(np.arange(18,36+18,1.), 36),  np.mod(np.arange(18,36+18,1.), 36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='Ocean')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.arange(0,36,1.)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_ww3(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(18,36+18,1.), 36), np.mod(np.arange(18,36+18,1.), 36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='WW3')
		Snew = copy(S)
		Dnew = copy(D)

		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)


		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_math(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(18,36+18,1.), 36), np.mod(np.arange(18,36+18,1.), 36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='Math')
		Snew = copy(S)
		Dnew = copy(D)

		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_mathvec(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(18,36+18,1.), 36), np.mod(np.arange(18,36+18,1.), 36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='MathVec')
		Snew = copy(S)
		Dnew = copy(D)

		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,1.), np.arange(0,36,1.)])))

class BoundarySpectralConventionsMath(unittest.TestCase):
	def test_math_to_ocean(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='Ocean')

		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.arange(0,36,1.)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_math_to_ww3(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='WW3')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_math_to_met(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='Met')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(18,36+18,1.), 36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_math_to_mathvec(self):
		D=np.arange(0,360,10.)
		S=np.array([np.mod(np.arange(0,-36,-1.)+9,36), np.mod(np.arange(0,-36,-1.)+9,36)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='MathVec')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-10.)+90,360)))
		stemp = np.arange(0,36,1.)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

class BoundarySpectralConventionsMathVec(unittest.TestCase):
	def test_mathvec_to_ocean(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Ocean')
		Snew, Dnew=bnd_processor(spec=S,dirs=D)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_mathvec_to_met(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Met')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		stemp = np.mod(np.arange(18,36+18,1.), 36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_mathvec_to_ww3(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='WW3')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_mathvec_to_math(self):
		D=np.mod(np.arange(0,-360,-10.)+90,360)
		S=np.array([np.arange(0,36,1.), np.arange(0,36,1.)])

		bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Math')

		Snew = copy(S)
		Dnew = copy(D)
		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)

		self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,10.)))
		stemp = np.mod(np.arange(0,-36,-1.)+9,36)
		self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

if __name__ == '__main__':
	unittest.main()
