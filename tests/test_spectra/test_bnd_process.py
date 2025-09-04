import unittest
import dnora
import numpy as np
from copy import copy
from dnora.type_manager.spectral_conventions import SpectralConvention
import os


def load_test_spec(shifted=False, math=False):
    f = np.linspace(0.0345, 0.5476, 30)  # np.loadtxt('data/freq.test')
    if shifted:
        if math:
            D = np.mod(
                np.linspace(82.5, -262.5, 24), 360
            )  # np.loadtxt('data/dir_math_shifted.test')
        else:
            D = np.linspace(7.5, 352.5, 24)  # np.loadtxt('data/dir_shifted.test')
    else:
        if math:
            D = np.mod(
                np.linspace(90.0, -255.0, 24), 360
            )  # np.loadtxt('data/dir_math.test')
        else:
            D = np.linspace(0.0, 345.0, 24)  # np.loadtxt('data/dir.test')

    S = np.ones((2, 2, len(f), len(D)), float)

    # Get the current file's directory (the test file's directory)
    current_dir = os.path.dirname(__file__)

    # Construct the absolute path to the data file
    data_file1 = os.path.join(current_dir, "data", "spec1.test")
    data_file2 = os.path.join(current_dir, "data", "spec2.test")
    data_file3 = os.path.join(current_dir, "data", "spec3.test")
    data_file4 = os.path.join(current_dir, "data", "spec4.test")
    S[0, 0, :, :] = np.loadtxt(data_file1)
    S[0, 1, :, :] = np.loadtxt(data_file2)
    S[1, 0, :, :] = np.loadtxt(data_file3)
    S[1, 1, :, :] = np.loadtxt(data_file4)

    return S, D


def loop_conventions(list_of_conventions, S, D):
    Snew = copy(S)
    inds = np.array(range(S.shape[0]))
    Dnew = copy(D)
    freq = np.linspace(0.0345, 0.5476, 30)
    for n in range(len(list_of_conventions) - 1):
        cur_c = list_of_conventions[n]
        wan_c = list_of_conventions[n + 1]
        bnd_processor = dnora.process.spectra.spectral_processor_for_convention_change(
            current_convention=cur_c, wanted_convention=wan_c
        )
        if not isinstance(bnd_processor, list):
            bnd_processor = [bnd_processor]

        for processor in bnd_processor:
            Snew, Dnew, _, _, _ = processor(
                spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
            )
    return Snew, Dnew


class BoundarySpectralConventionsOcean(unittest.TestCase):
    def test_ocean_to_met(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            inds = np.array(range(S.shape[0]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.OCEAN,
                    wanted_convention=SpectralConvention.MET,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.mod(np.arange(0, 36, dD / 10) + 18, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_ocean_to_ww3(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.OCEAN,
                    wanted_convention=SpectralConvention.WW3,
                )
            )

            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_ocean_to_math(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.OCEAN,
                    wanted_convention=SpectralConvention.MATH,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_ocean_to_mathvec(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.OCEAN,
                    wanted_convention=SpectralConvention.MATHVEC,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )

            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            self.assertIsNone(np.testing.assert_almost_equal(Snew, S))


class BoundarySpectralConventionsWW3(unittest.TestCase):
    def test_ww3_to_ocean(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.WW3,
                    wanted_convention=SpectralConvention.OCEAN,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Snew,
                    np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)]),
                )
            )

    def test_ww3_to_met(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.WW3,
                    wanted_convention=SpectralConvention.MET,
                )
            )
            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            stemp = np.mod(np.arange(0, 36, dD / 10) + 18, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_ww3_to_math(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.WW3,
                    wanted_convention=SpectralConvention.MATH,
                )
            )
            if not isinstance(bnd_processor, list):
                bnd_processor = [bnd_processor]

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_ww3_to_mathvec(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.WW3,
                    wanted_convention=SpectralConvention.MATHVEC,
                )
            )
            if not isinstance(bnd_processor, list):
                bnd_processor = [bnd_processor]

            Snew = copy(S)
            Dnew = copy(D)

            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Snew,
                    np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)]),
                )
            )


class BoundarySpectralConventionsMet(unittest.TestCase):
    def test_met_to_ocean(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MET,
                    wanted_convention=SpectralConvention.OCEAN,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.arange(0, 36, dD / 10)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_met_to_ww3(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MET,
                    wanted_convention=SpectralConvention.WW3,
                )
            )
            Snew = copy(S)
            Dnew = copy(D)

            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_met_to_math(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MET,
                    wanted_convention=SpectralConvention.MATH,
                )
            )
            Snew = copy(S)
            Dnew = copy(D)

            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            self.assertIsNone(
                np.testing.assert_almost_equal(Snew, np.array([stemp, stemp]))
            )

    def test_met_to_mathvec(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                    np.mod(np.arange(0, 36, dD / 10) + 18, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MET,
                    wanted_convention=SpectralConvention.MATHVEC,
                )
            )
            Snew = copy(S)
            Dnew = copy(D)

            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Snew,
                    np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)]),
                )
            )


class BoundarySpectralConventionsMath(unittest.TestCase):
    def test_math_to_ocean(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATH,
                    wanted_convention=SpectralConvention.OCEAN,
                )
            )

            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.arange(0, 36, dD / 10)
            mask = stemp > 35.999
            stemp[mask] = 0
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))

    def test_math_to_ww3(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATH,
                    wanted_convention=SpectralConvention.WW3,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_math_to_met(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATH,
                    wanted_convention=SpectralConvention.MET,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.mod(np.arange(0, 36, dD / 10) + 18, 36)
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))

    def test_math_to_mathvec(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.arange(0, 360, dD)
            S = np.array(
                [
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                    np.mod(np.arange(0, -36, -dD / 10) + 9, 36),
                ]
            )
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATH,
                    wanted_convention=SpectralConvention.MATHVEC,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(
                    Dnew, np.mod(np.arange(0, -360, -dD) + 90, 360)
                )
            )
            stemp = np.arange(0, 36, dD / 10)
            mask = stemp > 35.999
            stemp[mask] = 0
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))


class BoundarySpectralConventionsMathVec(unittest.TestCase):
    def test_mathvec_to_ocean(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATHVEC,
                    wanted_convention=SpectralConvention.OCEAN,
                )
            )
            Snew, Dnew, _, _, _ = bnd_processor(
                spec=S, dirs=D, freq=freq, inds=inds, times=None
            )

            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_mathvec_to_met(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATHVEC,
                    wanted_convention=SpectralConvention.MET,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            stemp = np.mod(np.arange(0, 36, dD / 10) + 18, 36)
            mask = stemp > 35.999
            stemp[mask] = 0
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))

    def test_mathvec_to_ww3(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATHVEC,
                    wanted_convention=SpectralConvention.WW3,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            mask = stemp > 35.999
            stemp[mask] = 0
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))

    def test_mathvec_to_math(self):
        for dD in [1.0, 2.0, 3.0, 5.0, 6.0, 9.0, 10.0, 15.0, 18.0, 22.5, 30.0, 45.0]:
            D = np.mod(np.arange(0, -360, -dD) + 90, 360)
            S = np.array([np.arange(0, 36, dD / 10), np.arange(0, 36, dD / 10)])
            inds = np.array(range(S.shape[0]))
            freq = np.array([[0.1]])
            S = np.reshape(S, (S.shape[0], len(freq), S.shape[-1]))
            bnd_processor = (
                dnora.process.spectra.spectral_processor_for_convention_change(
                    current_convention=SpectralConvention.MATHVEC,
                    wanted_convention=SpectralConvention.MATH,
                )
            )

            Snew = copy(S)
            Dnew = copy(D)
            for processor in bnd_processor:
                Snew, Dnew, _, _, _ = processor(
                    spec=Snew, dirs=Dnew, freq=freq, inds=inds, times=None
                )
            Snew = np.squeeze(Snew)
            self.assertIsNone(
                np.testing.assert_almost_equal(Dnew, np.arange(0, 360, dD))
            )
            stemp = np.mod(np.arange(0, -36, -dD / 10) + 9, 36)
            mask = stemp > 35.999
            stemp[mask] = 0
            self.assertIsNone(np.testing.assert_almost_equal(Snew, [stemp, stemp]))


class BoundarySpectralConventionsCircular(unittest.TestCase):
    def test_ocean_circular(self):
        S, D = load_test_spec()

        list_of_conventions = [
            SpectralConvention.OCEAN,
            SpectralConvention.MET,
            SpectralConvention.WW3,
            SpectralConvention.MATH,
            SpectralConvention.MATHVEC,
            SpectralConvention.OCEAN,
        ]

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        S, D = load_test_spec(shifted=True)  # Spectra starting from 7.5 deg

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_met_circular(self):
        S, D = load_test_spec()

        list_of_conventions = [
            SpectralConvention.MET,
            SpectralConvention.WW3,
            SpectralConvention.MATHVEC,
            SpectralConvention.OCEAN,
            SpectralConvention.MATH,
            SpectralConvention.MET,
        ]

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        S, D = load_test_spec(shifted=True)  # Spectra starting from 7.5 deg

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_ww3_circular(self):
        S, D = load_test_spec(math=True)

        list_of_conventions = [
            SpectralConvention.WW3,
            SpectralConvention.OCEAN,
            SpectralConvention.MATH,
            SpectralConvention.MATHVEC,
            SpectralConvention.MET,
            SpectralConvention.WW3,
        ]

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        S, D = load_test_spec(shifted=True, math=True)  # Spectra starting from 7.5 deg

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_math_circular(self):
        S, D = load_test_spec()

        list_of_conventions = [
            SpectralConvention.MATH,
            SpectralConvention.MATHVEC,
            SpectralConvention.WW3,
            SpectralConvention.OCEAN,
            SpectralConvention.MET,
            SpectralConvention.WW3,
            SpectralConvention.MATH,
        ]

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        S, D = load_test_spec(shifted=True)  # Spectra starting from 7.5 deg

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

    def test_mathvec_circular(self):
        S, D = load_test_spec(math=True)

        list_of_conventions = [
            SpectralConvention.MATHVEC,
            SpectralConvention.WW3,
            SpectralConvention.OCEAN,
            SpectralConvention.MET,
            SpectralConvention.MATH,
            SpectralConvention.WW3,
            SpectralConvention.MATHVEC,
        ]

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        S, D = load_test_spec(shifted=True, math=True)  # Spectra starting from 7.5 deg

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))

        list_of_conventions.reverse()

        Snew, Dnew = loop_conventions(list_of_conventions, S, D)
        self.assertIsNone(np.testing.assert_almost_equal(Dnew, D))
        self.assertIsNone(np.testing.assert_almost_equal(Snew, S))


if __name__ == "__main__":
    unittest.main()
