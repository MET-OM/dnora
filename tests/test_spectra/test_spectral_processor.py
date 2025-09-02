import dnora as dn
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


def test_ocean_to_ww3():
    S, D = load_test_spec()
    spectra = dn.spectra.Spectra(
        lon=(0, 1),
        lat=(10, 20),
        time=("2020-01-01 00:00", "2020-01-01 01:00"),
        freq=np.arange(30) / 10,
        dirs=D,
    )
    spectra.set_spec(S)
    spectra._mark_convention(SpectralConvention.OCEAN)
    spectra.set_convention("ww3")
    np.testing.assert_array_almost_equal(
        spectra.dirs(), np.mod(np.arange(90, -270, -15), 360)
    )


def test_ocean_to_ww3_shifted():
    S, D = load_test_spec(shifted=True)
    spectra = dn.spectra.Spectra(
        lon=(0, 1),
        lat=(10, 20),
        time=("2020-01-01 00:00", "2020-01-01 01:00"),
        freq=np.arange(30) / 10,
        dirs=D,
    )
    spectra.set_spec(S)

    spectra._mark_convention(SpectralConvention.OCEAN)

    spectra.set_convention("ww3")
    np.testing.assert_array_almost_equal(spectra.dirs()[0:6], np.flip(D[0:6]))
    np.testing.assert_array_almost_equal(spectra.dirs()[6:], np.flip(D[6:]))


def test_ocean_to_ww3_nora3():
    S, _ = load_test_spec(shifted=True)
    D = np.array(
        [
            7.5,
            22.5,
            37.499996,
            52.499996,
            67.5,
            82.49999,
            97.49999,
            112.49999,
            127.49999,
            142.49998,
            157.49998,
            172.49998,
            187.5,
            202.49998,
            217.49998,
            232.49998,
            247.49998,
            262.49997,
            277.49997,
            292.5,
            307.5,
            322.5,
            337.49997,
            352.49997,
        ]
    )
    spectra = dn.spectra.Spectra(
        lon=(0, 1),
        lat=(10, 20),
        time=("2020-01-01 00:00", "2020-01-01 01:00"),
        freq=np.arange(30) / 10,
        dirs=D,
    )
    spectra.set_spec(S)

    spectra._mark_convention(SpectralConvention.OCEAN)

    spectra.set_convention("ww3")
    np.testing.assert_array_almost_equal(spectra.dirs()[0:6], np.flip(D[0:6]))
    np.testing.assert_array_almost_equal(spectra.dirs()[6:], np.flip(D[6:]))
