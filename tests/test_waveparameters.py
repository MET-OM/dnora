from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D

import numpy as np
import pandas as pd
import pytest

from dnora.wave_parameters.parameters import get_function

import geo_parameters as gp


@pytest.fixture
def spec1d():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    spec = Spectra1D(lon=0, lat=0, freq=freq, time=time)
    spec.set_spec(2)
    return spec


@pytest.fixture
def spec():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    dirs = np.linspace(0, 355, 36)
    spec = Spectra1D(lon=0, lat=0, freq=freq, time=time, dirs=dirs)
    spec.set_spec(2)
    return spec


def test_m0(spec1d, spec):
    func = get_function(gp.wave.M0)
    assert func is not None
    mom0 = 2 * 0.1 * 10

    m0 = func(spec1d)
    np.testing.assert_almost_equal(np.mean(m0), mom0)

    m0 = func(spec)
    np.testing.assert_almost_equal(np.mean(m0), mom0)


def test_moms(spec1d, spec):
    moms = [gp.wave.M0, gp.wave.M1, gp.wave.M2, gp.wave.M3, gp.wave.M4, gp.wave.M5]
    for n, m in enumerate(moms):
        val = np.mean(spec1d.spec())
        mm = np.sum(val * spec1d.freq() ** n) * spec1d.df()
        func = get_function(m)
        assert func is not None
        func_mom1d = func(spec1d)
        func_mom = func(spec)
        np.testing.assert_almost_equal(np.mean(func_mom1d), mm)
        np.testing.assert_almost_equal(np.mean(func_mom), mm)

    mm = np.sum(val * spec1d.freq() ** -1) * spec1d.df()
    func = get_function(gp.wave.M_1)
    assert func is not None
    func_mom1d = func(spec1d)
    func_mom = func(spec)
    np.testing.assert_almost_equal(np.mean(func_mom1d), mm)
    np.testing.assert_almost_equal(np.mean(func_mom), mm)


def test_hs(spec1d, spec):
    func = get_function(gp.wave.Hs)
    assert func is not None
    hsig = 4 * np.sqrt(2)

    hs = func(spec1d)
    np.testing.assert_almost_equal(np.mean(hs), hsig)

    hs = func(spec)
    np.testing.assert_almost_equal(np.mean(hs), hsig)


def test_tm01(spec1d, spec):
    func01 = get_function(gp.wave.Tm01)
    assert func01 is not None
    val = np.mean(spec1d.spec())
    m1 = np.sum(val * spec1d.freq() ** 1) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm01_1d = func01(spec1d)
    tm01 = func01(spec)
    np.testing.assert_almost_equal(np.mean(tm01_1d), m0 / m1)
    np.testing.assert_almost_equal(np.mean(tm01), m0 / m1)


def test_tm02(spec1d, spec):
    func02 = get_function(gp.wave.Tm02)
    assert func02 is not None
    val = np.mean(spec1d.spec())
    m2 = np.sum(val * spec1d.freq() ** 2) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm02_1d = func02(spec1d)
    tm02 = func02(spec)
    np.testing.assert_almost_equal(np.mean(tm02_1d), m0 / m2)
    np.testing.assert_almost_equal(np.mean(tm02), m0 / m2)


def test_tm_10(spec1d, spec):
    func_10 = get_function(gp.wave.Tm_10)
    assert func_10 is not None
    val = np.mean(spec1d.spec())
    m_1 = np.sum(val * spec1d.freq() ** -1) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm_10_1d = func_10(spec1d)
    tm_10 = func_10(spec)
    np.testing.assert_almost_equal(np.mean(tm_10_1d), m_1 / m0)
    np.testing.assert_almost_equal(np.mean(tm_10), m_1 / m0)
