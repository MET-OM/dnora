import xarray as xr
import numpy as np

from dnora.spectra import Spectra
from dnora.spectra1d import Spectra1D
from typing import Union
import geo_parameters as gp
from typing import Union

from functools import partial


def moment(spec: Union[Spectra, Spectra1D], moment: float) -> xr.Dataset:
    if isinstance(spec, Spectra):
        dD = 360 / len(self._moment.dirs())
        ds = dD * np.pi / 180 * spec.ds().sum(dim="dirs")
    else:
        ds = spec.ds()

    if np.std(np.diff(ds.freq.data)) > 0.001:
        ds = (ds.spec * (ds.freq**moment)).integrate(coord="freq")
    else:
        ds = (ds.spec * (ds.freq**moment)).sum(dim="freq") * np.median(
            np.diff(ds.freq.data)
        )
    return ds.data


def hs(spec: Union[Spectra, Spectra1D]) -> xr.Dataset:
    func = get_function(gp.wave.M0)
    m0 = func(spec)
    return 4 * np.sqrt(m0)


def tm01(spec: Union[Spectra, Spectra1D]) -> xr.Dataset:
    func0 = get_function(gp.wave.M0)
    func1 = get_function(gp.wave.M1)
    m0 = func0(spec)
    m1 = func1(spec)
    return m0 / m1


def tm02(spec: Union[Spectra, Spectra1D]) -> xr.Dataset:
    func0 = get_function(gp.wave.M0)
    func2 = get_function(gp.wave.M2)
    m0 = func0(spec)
    m2 = func2(spec)
    return m0 / m2


def tm_10(spec: Union[Spectra, Spectra1D]) -> xr.Dataset:
    func0 = get_function(gp.wave.M0)
    func_1 = get_function(gp.wave.M_1)
    m0 = func0(spec)
    m_1 = func_1(spec)
    return m_1 / m0


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
}
