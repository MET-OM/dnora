from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D

import numpy as np
import pandas as pd
import pytest

from dnora.wave_parameters.parameters import get_function
from dnora.type_manager.spectral_conventions import SpectralConvention
import geo_parameters as gp


@pytest.fixture
def spec1d():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    spec = Spectra1D(lon=0, lat=0, freq=freq, time=time)
    spec.set_spec(2)
    spec.set_dirm(0)
    spec.set_spr(30)
    return spec


@pytest.fixture
def spec():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    dirs = np.linspace(0, 350, 36)
    spec = Spectra(lon=0, lat=0, freq=freq, time=time, dirs=dirs)
    spec_val = np.zeros(spec.shape("spec"))

    spec_val[:, :, :, 9] = 2 / spec.dd(angular=True)
    spec.set_spec(spec_val)
    spec._mark_convention(SpectralConvention.OCEAN)
    return spec


@pytest.fixture
def spec1dmono():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    spec = Spectra1D(lon=0, lat=0, freq=freq, time=time)
    spec_val = np.zeros(spec.shape("spec"))
    spec_val[:, :, 2] = 2
    spec.set_spec(spec_val)

    dirm_val = np.zeros(spec.shape("spec"))
    dirm_val[:, :, 2] = 0
    spec.set_dirm(dirm_val)

    spr_val = np.zeros(spec.shape("spec"))
    spr_val[:, :, 2] = 30
    spec.set_spr(spr_val)

    return spec


@pytest.fixture
def specmono():
    time = pd.date_range(start="2020-01-01 00:00", end="2020-01-02 00:00", freq="1h")
    freq = np.linspace(0.1, 1, 10)
    dirs = np.linspace(0, 350, 36)
    spec = Spectra(lon=0, lat=0, freq=freq, time=time, dirs=dirs)
    spec_val = np.zeros(spec.shape("spec"))

    spec_val[:, :, 2, 9] = 2 / spec.dd(angular=True)
    spec.set_spec(spec_val)
    spec._mark_convention(SpectralConvention.OCEAN)
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


def test_tm01(spec1d, spec, spec1dmono, specmono):
    func01 = get_function(gp.wave.Tm01)
    assert func01 is not None
    val = np.mean(spec1d.spec())
    m1 = np.sum(val * spec1d.freq() ** 1) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm01_1d = func01(spec1d)
    tm01 = func01(spec)
    np.testing.assert_almost_equal(np.mean(tm01_1d), m0 / m1)
    np.testing.assert_almost_equal(np.mean(tm01), m0 / m1)

    tm_01_1d = func01(spec1dmono)
    tm_01 = func01(specmono)
    np.testing.assert_almost_equal(tm_01_1d, 1 / 0.3)
    np.testing.assert_almost_equal(tm_01, 1 / 0.3)


def test_tm02(spec1d, spec, spec1dmono, specmono):
    func02 = get_function(gp.wave.Tm02)
    assert func02 is not None
    val = np.mean(spec1d.spec())
    m2 = np.sum(val * spec1d.freq() ** 2) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm02_1d = func02(spec1d)
    tm02 = func02(spec)
    np.testing.assert_almost_equal(np.mean(tm02_1d), np.sqrt(m0 / m2))
    np.testing.assert_almost_equal(np.mean(tm02), np.sqrt(m0 / m2))

    tm_02_1d = func02(spec1dmono)
    tm_02 = func02(specmono)
    np.testing.assert_almost_equal(tm_02_1d, 1 / 0.3)
    np.testing.assert_almost_equal(tm_02, 1 / 0.3)


def test_tm_10(spec1d, spec, spec1dmono, specmono):
    func_10 = get_function(gp.wave.Tm_10)
    assert func_10 is not None
    val = np.mean(spec1d.spec())
    m_1 = np.sum(val * spec1d.freq() ** -1) * spec1d.df()
    m0 = np.sum(val * spec1d.freq() ** 0) * spec1d.df()
    tm_10_1d = func_10(spec1d)
    tm_10 = func_10(spec)
    np.testing.assert_almost_equal(np.mean(tm_10_1d), m_1 / m0)
    np.testing.assert_almost_equal(np.mean(tm_10), m_1 / m0)

    tm_10_1d = func_10(spec1dmono)
    tm_10 = func_10(specmono)
    np.testing.assert_almost_equal(tm_10_1d, 1 / 0.3)
    np.testing.assert_almost_equal(tm_10, 1 / 0.3)


def test_dirm(spec, spec1d, specmono, spec1dmono):
    func = get_function(gp.wave.Dirm)
    assert func is not None
    dirm = func(spec)
    np.testing.assert_almost_equal(np.mean(dirm), 270)
    dirm = func(spec1d)
    np.testing.assert_almost_equal(np.mean(dirm), 0)

    func = get_function(gp.wave.DirmTo)
    assert func is not None
    dirm = func(spec)
    np.testing.assert_almost_equal(np.mean(dirm), 90)
    dirm = func(spec1d)
    np.testing.assert_almost_equal(np.mean(dirm), 180)

    dirm = func(specmono)
    np.testing.assert_almost_equal(np.mean(dirm), 90)
    dirm = func(spec1dmono)
    np.testing.assert_almost_equal(np.mean(dirm), 180)


def test_sprm(spec, spec1d, specmono, spec1dmono):
    func = get_function(gp.wave.Spr)
    assert func is not None
    sprm = func(spec)
    np.testing.assert_almost_equal(np.mean(sprm), 0, decimal=5)
    sprm = func(spec1d)
    np.testing.assert_almost_equal(np.mean(sprm), 30, decimal=5)

    sprm = func(specmono)
    np.testing.assert_almost_equal(np.mean(sprm), 0)
    sprm = func(spec1dmono)
    np.testing.assert_almost_equal(np.mean(sprm), 30)


def test_fp(spec1dmono, specmono):
    func = get_function(gp.wave.Fp)
    assert func is not None
    fp = func(spec1dmono)
    np.testing.assert_almost_equal(fp, 0.3, decimal=5)
    fp = func(specmono)
    np.testing.assert_almost_equal(fp, 0.3, decimal=5)


def test_fm(spec1dmono, specmono):
    func = get_function(gp.wave.Fm)
    assert func is not None
    fm = func(spec1dmono)
    np.testing.assert_almost_equal(fm, 0.3, decimal=5)
    fm = func(specmono)
    np.testing.assert_almost_equal(fm, 0.3, decimal=5)


def test_wp(spec1dmono, specmono):
    func = get_function(gp.wave.Wp)
    assert func is not None
    wp = func(spec1dmono)
    np.testing.assert_almost_equal(wp, 2 * np.pi * 0.3, decimal=5)
    wp = func(specmono)
    np.testing.assert_almost_equal(wp, 2 * np.pi * 0.3, decimal=5)


def test_wm(spec1dmono, specmono):
    func = get_function(gp.wave.Wm)
    assert func is not None
    wm = func(spec1dmono)
    np.testing.assert_almost_equal(wm, 2 * np.pi * 0.3, decimal=5)
    wm = func(specmono)
    np.testing.assert_almost_equal(wm, 2 * np.pi * 0.3, decimal=5)


def test_tp(spec1dmono, specmono):
    func = get_function(gp.wave.Tp)
    assert func is not None
    tp = func(spec1dmono)
    np.testing.assert_almost_equal(tp, 1 / 0.3, decimal=5)

    tp = func(specmono)
    np.testing.assert_almost_equal(tp, 1 / 0.3, decimal=5)


def test_dirp(spec1dmono, specmono):
    func = get_function(gp.wave.Dirp)
    assert func is not None
    dirp = func(spec1dmono)
    np.testing.assert_almost_equal(dirp, 0, decimal=5)

    dirp = func(specmono)
    np.testing.assert_almost_equal(dirp, 270, decimal=5)

    func = get_function(gp.wave.DirpTo)
    assert func is not None
    dirp = func(spec1dmono)
    np.testing.assert_almost_equal(dirp, 180, decimal=5)

    dirp = func(specmono)
    np.testing.assert_almost_equal(dirp, 90, decimal=5)
