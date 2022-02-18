import unittest
import sys
sys.path.insert(0, "../../")
from dnora import bnd
import numpy as np
from copy import copy

def load_test_spec(shifted=False, math=False):
	f = np.loadtxt('data/freq.test')
	if shifted:
		if math:
			D = np.loadtxt('data/dir_math_shifted.test')
		else:
			D = np.loadtxt('data/dir_shifted.test')
	else:
		if math:
			D = np.loadtxt('data/dir_math.test')
		else:
			D = np.loadtxt('data/dir.test')

	S=np.ones((2,2,len(f),len(D)), float)
	S[0,0,:,:]=np.loadtxt('data/spec1.test')
	S[0,1,:,:]=np.loadtxt('data/spec2.test')
	S[1,0,:,:]=np.loadtxt('data/spec3.test')
	S[1,1,:,:]=np.loadtxt('data/spec4.test')

	return S, D

def loop_conventions(list_of_conventions, S, D):
	Snew = copy(S)
	Dnew = copy(D)
	for n in range(len(list_of_conventions)-1):
		cur_c = list_of_conventions[n]
		wan_c = list_of_conventions[n+1]
		bnd_processor = bnd.process.processor_for_convention_change(current_convention=cur_c, wanted_convention=wan_c)
		if not isinstance(bnd_processor, list):
			bnd_processor = [bnd_processor]

		for processor in bnd_processor:
			Snew, Dnew=processor(spec=Snew,dirs=Dnew)
	return Snew, Dnew

class BoundarySpectralConventionsOcean(unittest.TestCase):
	def test_ocean_to_met(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='Met')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.mod(np.arange(0,36,dD/10)+18,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ocean_to_ww3(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='WW3')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ocean_to_math(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='Math')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ocean_to_mathvec(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Ocean', wanted_convention='MathVec')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

class BoundarySpectralConventionsWW3(unittest.TestCase):
	def test_ww3_to_ocean(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Ocean')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])))

	def test_ww3_to_met(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Met')
			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			stemp = np.mod(np.arange(0,36,dD/10)+18,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))


	def test_ww3_to_math(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='Math')
			if not isinstance(bnd_processor, list):
				bnd_processor = [bnd_processor]

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_ww3_to_mathvec(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='WW3', wanted_convention='MathVec')
			if not isinstance(bnd_processor, list):
				bnd_processor = [bnd_processor]

			Snew = copy(S)
			Dnew = copy(D)

			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])))

class BoundarySpectralConventionsMet(unittest.TestCase):

	def test_met_to_ocean(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([ np.mod(np.arange(0,36,dD/10)+18,36),  np.mod(np.arange(0,36,dD/10)+18,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='Ocean')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.arange(0,36,dD/10)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_ww3(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,36,dD/10)+18,36), np.mod(np.arange(0,36,dD/10)+18,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='WW3')
			Snew = copy(S)
			Dnew = copy(D)

			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_math(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,36,dD/10)+18,36), np.mod(np.arange(0,36,dD/10)+18,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='Math')
			Snew = copy(S)
			Dnew = copy(D)

			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([stemp, stemp])))

	def test_met_to_mathvec(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,36,dD/10)+18,36), np.mod(np.arange(0,36,dD/10)+18,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Met', wanted_convention='MathVec')
			Snew = copy(S)
			Dnew = copy(D)

			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])))

class BoundarySpectralConventionsMath(unittest.TestCase):
	def test_math_to_ocean(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='Ocean')

			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.arange(0,36,dD/10)
			mask=stemp>35.999
			stemp[mask]=0
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_math_to_ww3(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='WW3')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_math_to_met(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='Met')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.mod(np.arange(0,36,dD/10)+18,36)
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_math_to_mathvec(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.arange(0,360,dD)
			S=np.array([np.mod(np.arange(0,-36,-dD/10)+9,36), np.mod(np.arange(0,-36,-dD/10)+9,36)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='Math', wanted_convention='MathVec')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.mod(np.arange(0,-360,-dD)+90,360)))
			stemp = np.arange(0,36,dD/10)
			mask=stemp>35.999
			stemp[mask]=0
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

class BoundarySpectralConventionsMathVec(unittest.TestCase):
	def test_mathvec_to_ocean(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Ocean')
			Snew, Dnew=bnd_processor(spec=S,dirs=D)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_mathvec_to_met(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Met')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			stemp = np.mod(np.arange(0,36,dD/10)+18,36)
			mask=stemp>35.999
			stemp[mask]=0
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_mathvec_to_ww3(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='WW3')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			mask=stemp>35.999
			stemp[mask]=0
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

	def test_mathvec_to_math(self):
		for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
			D=np.mod(np.arange(0,-360,-dD)+90,360)
			S=np.array([np.arange(0,36,dD/10), np.arange(0,36,dD/10)])

			bnd_processor = bnd.process.processor_for_convention_change(current_convention='MathVec', wanted_convention='Math')

			Snew = copy(S)
			Dnew = copy(D)
			for processor in bnd_processor:
				Snew, Dnew=processor(spec=Snew,dirs=Dnew)

			self.assertIsNone(np.testing.assert_almost_equal(Dnew,np.arange(0,360,dD)))
			stemp = np.mod(np.arange(0,-36,-dD/10)+9,36)
			mask=stemp>35.999
			stemp[mask]=0
			self.assertIsNone(np.testing.assert_almost_equal(Snew,[stemp, stemp]))

class BoundarySpectralConventionsCircular(unittest.TestCase):
	def test_ocean_circular(self):
		S, D = load_test_spec()

		list_of_conventions = ['Ocean', 'Met', 'WW3', 'Math', 'MathVec', 'Ocean']

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		S, D = load_test_spec(shifted=True) # Spectra starting from 7.5 deg

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_met_circular(self):
		S, D = load_test_spec()

		list_of_conventions = ['Met', 'WW3', 'MathVec', 'Ocean', 'Math', 'Met']

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		S, D = load_test_spec(shifted=True) # Spectra starting from 7.5 deg

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_ww3_circular(self):
		S, D = load_test_spec(math=True)

		list_of_conventions = ['WW3', 'Ocean', 'Math', 'MathVec', 'Met', 'WW3']

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		S, D = load_test_spec(shifted=True, math=True) # Spectra starting from 7.5 deg

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_math_circular(self):
		S, D = load_test_spec()

		list_of_conventions = ['Math', 'MathVec', 'WW3', 'Ocean', 'Met', 'WW3', 'Math']

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		S, D = load_test_spec(shifted=True) # Spectra starting from 7.5 deg

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

	def test_mathvec_circular(self):
		S, D = load_test_spec(math=True)

		list_of_conventions = ['MathVec', 'WW3', 'Ocean', 'Met', 'Math',  'WW3', 'MathVec']

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		S, D = load_test_spec(shifted=True, math=True) # Spectra starting from 7.5 deg

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))

		list_of_conventions.reverse()

		Snew, Dnew = loop_conventions(list_of_conventions, S, D)
		self.assertIsNone(np.testing.assert_almost_equal(Dnew,D))
		self.assertIsNone(np.testing.assert_almost_equal(Snew,S))


if __name__ == '__main__':
	unittest.main()
