import pytest
from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
import dnora as dn
import numpy as np


@pytest.fixture
def spec():
    spec = Spectra(
        lon=(1, 2, 3),
        lat=(10, 11, 12),
        time=("2020-01-01 00:00", "2020-01-01 02:00"),
        freq=[0.1, 0.2, 0.3],
        dirs=range(0, 360, 10),
    )
    spec.set_spec(1.0)
    return spec


@pytest.fixture
def spec1d():
    spec = Spectra1D(
        lon=(1, 2, 3),
        lat=(10, 11, 12),
        time=("2020-01-01 00:00", "2020-01-01 02:00"),
        freq=[0.1, 0.2, 0.3],
    )
    spec.set_spec(1.0)
    spec.set_dirm(90.0)
    spec.set_spr(10.0)
    return spec


def test_multiply(spec):
    np.testing.assert_almost_equal(np.max(spec.spec()), 1)
    spec.process(dn.process.spectra.Multiply(2))
    np.testing.assert_almost_equal(np.max(spec.spec()), 2)


def test_multiply1d(spec1d):
    np.testing.assert_almost_equal(np.max(spec1d.spec()), 1)
    spec1d.process(dn.process.spectra.Multiply(2))
    np.testing.assert_almost_equal(np.max(spec1d.spec()), 2)


def test_cutfreq(spec):
    np.testing.assert_almost_equal(spec.freq(), np.array([0.1, 0.2, 0.3]))
    spec.process(dn.process.spectra.CutFrequency(freq=(0.05, 0.25)))
    np.testing.assert_almost_equal(spec.freq(), np.array([0.1, 0.2]))


def test_cutfreq1d(spec1d):
    np.testing.assert_almost_equal(spec1d.freq(), np.array([0.1, 0.2, 0.3]))
    spec1d.process(dn.process.spectra.CutFrequency(freq=(0.05, 0.25)))
    np.testing.assert_almost_equal(spec1d.freq(), np.array([0.1, 0.2]))


def test_remove_empty(spec):
    np.testing.assert_almost_equal(spec.lon(), np.array([1, 2, 3]))
    np.testing.assert_almost_equal(spec.lat(), np.array([10, 11, 12]))
    data = spec.spec()
    data[:, 0, :, :] = np.nan
    data[:, 1, :, :] = 0.009
    spec.process(dn.process.spectra.RemoveEmpty(threshold=0.01))
    np.testing.assert_almost_equal(spec.lon(), np.array([3]))
    np.testing.assert_almost_equal(spec.lat(), np.array([12]))


def test_remove_empty_1d(spec1d):
    np.testing.assert_almost_equal(spec1d.lon(), np.array([1, 2, 3]))
    np.testing.assert_almost_equal(spec1d.lat(), np.array([10, 11, 12]))
    data = spec1d.spec()
    data[:, 0, :] = np.nan
    data[:, 1, :] = 0.009
    spec1d.process(dn.process.spectra.RemoveEmpty(threshold=0.01))
    np.testing.assert_almost_equal(spec1d.lon(), np.array([3]))
    np.testing.assert_almost_equal(spec1d.lat(), np.array([12]))


def test_remove_nantimes(spec):
    assert spec.time(datetime=False) == [
        "2020-01-01 00:00:00",
        "2020-01-01 01:00:00",
        "2020-01-01 02:00:00",
    ]
    data = spec.spec()
    data[1, 0, :, :] = np.nan
    spec.process(dn.process.spectra.RemoveNanTimes())
    assert spec.time(datetime=False) == [
        "2020-01-01 00:00:00",
        "2020-01-01 02:00:00",
    ]


def test_remove_nantimes_1d(spec1d):
    assert spec1d.time(datetime=False) == [
        "2020-01-01 00:00:00",
        "2020-01-01 01:00:00",
        "2020-01-01 02:00:00",
    ]
    data = spec1d.spec()
    data[1, 0, :] = np.nan
    spec1d.process(dn.process.spectra.RemoveNanTimes())
    assert spec1d.time(datetime=False) == [
        "2020-01-01 00:00:00",
        "2020-01-01 02:00:00",
    ]
