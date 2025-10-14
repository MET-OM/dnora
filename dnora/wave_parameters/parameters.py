import xarray as xr
import numpy as np

from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from dnora.type_manager.spectral_conventions import SpectralConvention
from typing import Union
import geo_parameters as gp
from typing import Union

from functools import partial

from copy import deepcopy


def trapz_or_sum(da):
    """Integrates using trapezodian integration if non-monotone frequency spacing and integrates by summing if monotone frequency spacing"""
    if np.std(np.diff(da.freq.data)) > 0.001:
        return da.integrate(coord="freq")
    else:
        return da.sum(dim="freq") * np.median(np.diff(da.freq.data))


def cos_sin(spec: Union[Spectra, Spectra1D]) -> tuple[xr.DataArray]:
    if isinstance(spec, Spectra):
        theta = np.deg2rad(spec.ds().dirs)
        dD = 360 / len(spec.dirs())
        # Normalizing here so that integration over direction becomes summing
        efth = dD * np.pi / 180 * spec.spec(data_array=True, squeeze=False)
        c1 = (np.cos(theta) * efth).sum(dim="dirs")  # Function of frequency
        s1 = (np.sin(theta) * efth).sum(dim="dirs")
    else:
        theta = np.deg2rad(spec.dirm(data_array=True, squeeze=False))
        c1 = np.cos(theta) * spec.spec(
            data_array=True, squeeze=False
        )  # Function of frequency
        s1 = np.sin(theta) * spec.spec(data_array=True, squeeze=False)

    return c1, s1


def first_fourier_coefficients(spec: Union[Spectra, Spectra1D]) -> tuple[xr.DataArray]:
    c1, s1 = cos_sin(spec)
    m0 = moment(spec, moment=0)
    a1m = trapz_or_sum(c1.fillna(0)) / m0
    b1m = trapz_or_sum(s1.fillna(0)) / m0

    return a1m, b1m


def moment(spec: Union[Spectra, Spectra1D], moment: float) -> np.ndarray:
    if isinstance(spec, Spectra):
        dD = 360 / len(spec.dirs())
        ds = dD * np.pi / 180 * spec.ds().sum(dim="dirs")
    else:
        ds = spec.ds()
    ds = trapz_or_sum(ds.spec * (ds.freq**moment))

    return ds.data


def hs(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func = get_function(gp.wave.M0)
    m0 = func(spec)
    return 4 * np.sqrt(m0)


def tm01(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func0 = get_function(gp.wave.M0)
    func1 = get_function(gp.wave.M1)
    m0 = func0(spec)
    m1 = func1(spec)
    return m0 / m1


def tm02(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func0 = get_function(gp.wave.M0)
    func2 = get_function(gp.wave.M2)
    m0 = func0(spec)
    m2 = func2(spec)
    return np.sqrt(m0 / m2)


def tm_10(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func0 = get_function(gp.wave.M0)
    func_1 = get_function(gp.wave.M_1)
    m0 = func0(spec)
    m_1 = func_1(spec)
    return m_1 / m0


def fm(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func0 = get_function(gp.wave.M0)
    func1 = get_function(gp.wave.M1)
    m0 = func0(spec)
    m1 = func1(spec)
    return m1 / m0


def wp(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func = get_function(gp.wave.Fp)
    return 2 * np.pi * func(spec)


def wm(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    func = get_function(gp.wave.Fm)
    return 2 * np.pi * func(spec)


def dirm_from(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    if isinstance(spec, Spectra):
        spec.set_convention(SpectralConvention.MET)
    a1m, b1m = first_fourier_coefficients(spec)
    thetam = np.arctan2(b1m, a1m)
    dirm = np.mod(thetam * 180 / np.pi, 360)

    return dirm.data


def dirm_to(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    return np.mod(dirm_from(spec) + 180, 360)


def sprm(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    m0 = moment(spec, moment=0)

    if isinstance(spec, Spectra):
        a1m, b1m = first_fourier_coefficients(spec)
        m1 = np.sqrt(b1m**2 + a1m**2)
        # Can be e-16 negative which leads to nan, so take absolute value
        sprm = np.sqrt(np.abs(2 - 2 * (m1))) * 180 / np.pi

    else:
        sprm = (
            trapz_or_sum(
                spec.spr(squeeze=False) ** 2 * spec.spec(data_array=True, squeeze=False)
            )
            / m0
        ) ** 0.5

    return sprm.data


def fp(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    if isinstance(spec, Spectra):
        dD = 360 / len(spec.dirs())
        ds = dD * np.pi / 180 * spec.ds().sum(dim="dirs")
    else:
        ds = spec.ds()
    inds = ds.spec.argmax(dim="freq")
    # To inherit DataArray to make slicing work
    freqs = deepcopy(ds.spec)
    freqs.data = np.tile(spec.freq(), [spec.shape("spec")[0], spec.shape("spec")[1], 1])
    fp = freqs[:, :, inds]
    return fp.data


def tp(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    fpeak = fp(spec)
    return 1 / fpeak


def dirp_from(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    if isinstance(spec, Spectra):
        spec.set_convention(SpectralConvention.MET)
        c1, s1 = cos_sin(spec)
        # efth normalization is wrong, but since it is just a constant it won't affect argmax
        efth = spec.spec(data_array=True, dask=False, squeeze=False).sum(dim="dirs")
        dirs = np.rad2deg(np.arctan2(s1, c1))
    else:
        efth = spec.spec(data_array=True, dask=False, squeeze=False)
        dirs = spec.dirm(data_array=True, squeeze=False)

    inds = efth.argmax(dim="freq")

    dirs = dirs[:, :, inds].data
    mask = dirs < 0
    dirs[mask] = dirs[mask] + 360
    return dirs


def dirp_to(spec: Union[Spectra, Spectra1D]) -> np.ndarray:
    return np.mod(dirp_from(spec) + 180, 360)


def get_function(param: gp.metaparameter.MetaParameter) -> callable:
    return dict_of_functions.get(param.standard_name)


dict_of_functions = {
    gp.wave.M0.standard_name: partial(moment, moment=0),
    gp.wave.M1.standard_name: partial(moment, moment=1),
    gp.wave.M_1.standard_name: partial(moment, moment=-1),
    gp.wave.M2.standard_name: partial(moment, moment=2),
    gp.wave.M3.standard_name: partial(moment, moment=3),
    gp.wave.M4.standard_name: partial(moment, moment=4),
    gp.wave.M5.standard_name: partial(moment, moment=5),
    gp.wave.Hs.standard_name: hs,
    gp.wave.Tm01.standard_name: tm01,
    gp.wave.Tm02.standard_name: tm02,
    gp.wave.Tm_10.standard_name: tm_10,
    gp.wave.Tp.standard_name: tp,
    gp.wave.Fp.standard_name: fp,
    gp.wave.Fm.standard_name: fm,
    gp.wave.Wp.standard_name: wp,
    gp.wave.Wm.standard_name: wm,
    gp.wave.Dirm.standard_name: dirm_from,
    gp.wave.DirmTo.standard_name: dirm_to,
    gp.wave.Spr.standard_name: sprm,
    gp.wave.Dirp.standard_name: dirp_from,
    gp.wave.DirpTo.standard_name: dirp_to,
}
