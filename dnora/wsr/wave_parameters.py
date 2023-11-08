from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from ..bnd.conventions import SpectralConvention
from ..bnd import Boundary
from ..spc import Spectra
from typing import Union
from pint import UnitRegistry

ureg = UnitRegistry()
ureg.default_format = "~S"


class WaveParameter(ABC):
    """Calculates a wave parameter from spectra xarray

    Spectra xarray Dataset contains the following data:

    freq:   Frequency
    x:      Point index [int starting from 0]
    lon:    Longitude with length of x
    lat:    Latitude with length of x
    time:   Times
    spec:   Variance density [m^2/Hz]
    dirm:   Mean direction as a function of frequency [deg]
    spr:    Directional spreading as a function of frequency [deg]

    Boundary xarray Dataset contains the following data:

    freq:   Frequency
    dirs:   Directions [deg]
    x:      Point index [int starting from 0]
    lon:    Longitude with length of x
    lat:    Latitude with length of x
    time:   Times
    spec:   Variance density [m^2/Hz/rad]
    """

    @abstractmethod
    def __call__(self, spec: xr.Dataset):
        """Calculates parameter from a 1D-spectra"""
        pass

    @abstractmethod
    def name(self) -> str:
        pass

    @abstractmethod
    def unit(self) -> str:
        pass

    def standard_name(self) -> str:
        """"""


class Moment(WaveParameter):
    """Spectral moment from spectra"""

    def __init__(self, moment: int) -> None:
        self._moment = moment
        pass

    def __call__(self, spec: Union[Spectra, Boundary]) -> np.ndarray:
        if isinstance(spec, Boundary):
            dD = 360 / len(self._moment.dirs())
            ds = dD * np.pi / 180 * spec.ds().sum(dim="dirs")
        else:
            ds = spec.ds()

        ds = (ds.spec * (ds.freq**self._moment)).integrate(coord="freq")

        return ds.values

    def unit(self):
        return ureg.meter**2 * pow(ureg.second, -self._moment)

    def name(self):
        """1 returns 'm1', 0.5 returns 'm05' etc."""
        strs = f"{self._moment:.1f}".split(".")
        if strs[1] == "0":
            return f"m{strs[0]}"
        else:
            return f"m{strs[0]}{strs[1]}"


class PowerMoment(WaveParameter):
    """Spectral moment from spectra meighted by spectra any power"""

    def __init__(self, moment: float, power: float) -> None:
        self._moment = moment
        self._power = power

    def __call__(self, spec: Union[Spectra, Boundary]) -> np.ndarray:
        # moment = spec
        if isinstance(spec, Boundary):
            dD = 360 / len(spec.dirs())
            ds = dD * np.pi / 180 * spec.ds().sum(dim="dirs")
        else:
            ds = spec.ds()

        ds = (ds.spec**self._power * (ds.freq**self._moment)).integrate(
            coord="freq"
        )
        # breakpoint()
        # moment = self._format_dataset(moment, spec)

        return ds.values

    def unit(self):
        return (
            pow(ureg.meter, self._power * 2)
            * pow(ureg.second, self._power)
            * ureg.second**-1
            * pow(ureg.second, -self._moment)
        )

    def name(self):
        """1 returns 'm1', 0.5 returns 'm05' etc."""
        strs = f"{self._moment:.1f}".split(".")
        if strs[1] == "0":
            strs = f"m{strs[0]}"
        else:
            strs = f"m{strs[0]}{strs[1]}"

        strs2 = f"{self._power:.1f}".split(".")
        if strs2[1] == "0":
            strs += f"p{strs2[0]}"
        else:
            strs += f"p{strs2[0]}{strs2[1]}"

        return strs


class Hs(WaveParameter):
    """Singificant wave height from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 4 * Moment(0)(spec) ** 0.5

    def name(self):
        return "hs"

    def unit(self):
        return Moment(0).unit() ** 0.5

    def standard_name(self):
        return "sea_surface_wave_significant_height"


class Tm01(WaveParameter):
    """Mean wave period from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return Moment(0)(spec) / Moment(1)(spec)

    def name(self):
        return "tm01"

    def unit(self):
        return Moment(0).unit() / Moment(1).unit()

    def standard_name(self):
        return "sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"


class Tm_10(WaveParameter):
    """Mean wave period based on inverse moment from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return Moment(-1)(spec) / Moment(0)(spec)

    def name(self):
        return "tm_10"

    def unit(self):
        return "s"

    def standard_name(self):
        return "sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"


class Tm02(WaveParameter):
    """Mean wave period based on second moment from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return (Moment(0)(spec) / Moment(2)(spec)) ** (0.5)

    def name(self):
        return "tm02"

    def unit(self):
        return (Moment(0).unit() / Moment(2).unit()) ** 0.5

    def standard_name(self):
        return "sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment"


class Tp(WaveParameter):
    """Peak period from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        spec.set_convention(SpectralConvention.MET)

        if isinstance(spec, Boundary):
            theta = np.deg2rad(spec.mdir(data_array=True))
            dD = 360 / len(spec.dirs())
            # Normalizing here so that integration over direction becomes summing
            efth = dD * np.pi / 180 * spec.spec(data_array=True)

            c1 = (np.cos(theta) * efth).sum(dim="dirs")  # Function of frequency
            s1 = (np.sin(theta) * efth).sum(dim="dirs")
            efth = efth.sum(dim="dirs")
            dirs = np.rad2deg(np.arctan2(s1, c1))
        else:
            efth = spec.spec(data_array=True)
            dirs = spec.mdir(data_array=True)

        inds = efth.argmax(dim="freq")

        freqs = np.tile(efth.freq.values, [efth.shape[0], efth.shape[1], 1])
        dirs.values = freqs  # Inherit DataArray structure

        return 1 / dirs[:, :, inds].values

    def name(self):
        return "tp"

    def unit(self):
        return ureg.second

    def standard_name(self):
        return "sea_surface_wave_period_at_variance_spectral_density_maximum"


class TpI(WaveParameter):
    """Peak period from spectra from integration"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        m_1 = PowerMoment(-1, 10)(spec)
        mask = np.where(m_1 < 0.0000000001)
        m_1[mask] = -999
        tp = m_1 / PowerMoment(0, 10)(spec)
        tp[mask] = 0
        return tp

    def name(self):
        return "tpi"

    def unit(self):
        return PowerMoment(-1, 10).unit() / PowerMoment(0, 10).unit()

    def standard_name(self):
        return "sea_surface_wave_period_at_variance_spectral_density_maximum"


class Fm(WaveParameter):
    """Mean freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        m0 = Moment(0)(spec)
        mask = np.where(m0 < 0.0000000001)
        m0[mask] = -999
        fm = Moment(1)(spec) / m0
        fm[mask] = 0
        return fm

    def name(self):
        return "fm"

    def unit(self):
        return Moment(1).unit() / Moment(0).unit()


class Wm(WaveParameter):
    """Mean angular freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return (Moment(1)(spec) / Moment(0)(spec)) * 2 * np.pi

    def name(self):
        return "wm"

    def unit(self):
        return Moment(1).unit() / Moment(0).unit() * ureg.rad


class Fc(WaveParameter):
    """Characteristic freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return PowerMoment(1, 4)(spec) / PowerMoment(0, 4)(spec)

    def name(self):
        return "fc"

    def unit(self):
        return PowerMoment(1, 4).unit() / PowerMoment(0, 4).unit()


class Wc(WaveParameter):
    """Characteristic angular freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 2 * np.pi * PowerMoment(1, 4)(spec) / PowerMoment(0, 4)(spec)

    def name(self):
        return "wc"

    def unit(self):
        return PowerMoment(1, 4).unit() / PowerMoment(0, 4).unit() * ureg.rad


class Fp(WaveParameter):
    """Peak frequency from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return PowerMoment(1, 10)(spec) / PowerMoment(0, 10)(spec)

    def name(self):
        return "fp"

    def unit(self):
        return PowerMoment(1, 10).unit() / PowerMoment(0, 10).unit()

    def standard_name(self):
        return "sea_surface_wave_frequency_at_variance_spectral_density_maximum"


class Wp(WaveParameter):
    """Peak angular frequency from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 2 * np.pi * PowerMoment(1, 10)(spec) / PowerMoment(0, 10)(spec)

    def name(self):
        return "wp"

    def unit(self):
        return PowerMoment(1, 10).unit() / PowerMoment(0, 10).unit() * ureg.rad

    def standard_name(self):
        return "sea_surface_wave_angular_frequency_at_variance_spectral_density_maximum"


class Dirm(WaveParameter):
    """Mean wave direction"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        spec.set_convention(SpectralConvention.MET)

        if isinstance(spec, Boundary):
            theta = np.deg2rad(spec.mdir(data_array=True))
            dD = 360 / len(spec.dirs())
            # Normalizing here so that integration over direction becomes summing
            efth = dD * np.pi / 180 * spec.spec(data_array=True)

            c1 = (np.cos(theta) * efth).sum(dim="dirs")  # Function of frequency
            s1 = (np.sin(theta) * efth).sum(dim="dirs")
        else:
            theta = np.deg2rad(spec.mdir(data_array=True))

            c1 = np.cos(theta) * spec.spec(data_array=True)  # Function of frequency
            s1 = np.sin(theta) * spec.spec(data_array=True)

        m0 = Moment(0)(spec)

        a1m = c1.integrate(coord="freq") / m0  # Mean parameters
        b1m = s1.integrate(coord="freq") / m0

        thetam = np.arctan2(b1m, a1m)
        dirm = np.mod(thetam * 180 / np.pi, 360)

        return dirm.values

    def name(self):
        return "dirm"

    def unit(self):
        return ureg.deg

    def standard_name(self):
        return "sea_surface_wave_from_direction"


class Sprm(WaveParameter):
    """Mean wave spreading"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        m0 = Moment(0)(spec)

        if isinstance(spec, Boundary):
            theta = np.deg2rad(spec.dirs())
            dD = 360 / len(spec.dirs())
            # Normalizing here so that integration over direction becomes summing
            efth = dD * np.pi / 180 * spec.ds()

            c1 = (np.cos(theta) * efth).sum(dim="dirs")  # Function of frequency
            s1 = (np.sin(theta) * efth).sum(dim="dirs")

            a1m = c1.integrate(coord="freq") / m0  # Mean parameters
            b1m = s1.integrate(coord="freq") / m0

            m1 = np.sqrt(b1m**2 + a1m**2)
            sprm = np.sqrt(2 - 2 * (m1)) * 180 / np.pi

        else:
            sprm = (
                (spec.spr() ** 2 * spec.spec(data_array=True)).integrate(coord="freq")
                / m0
            ) ** 0.5

        # sprm = self._format_dataset(sprm, spec)
        return sprm.values

    def name(self):
        return "sprm"

    def unit(self):
        return ureg.deg

    def standard_name(self):
        return "sea_surface_wave_directional_spread"


class Dirp(WaveParameter):
    """Peak wave direction"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        spec.set_convention(SpectralConvention.MET)

        if isinstance(spec, Boundary):
            theta = np.deg2rad(spec.mdir(data_array=True))
            dD = 360 / len(spec.dirs())
            # Normalizing here so that integration over direction becomes summing
            efth = dD * np.pi / 180 * spec.spec(data_array=True)

            c1 = (np.cos(theta) * efth).sum(dim="dirs")  # Function of frequency
            s1 = (np.sin(theta) * efth).sum(dim="dirs")
            efth = efth.sum(dim="dirs")
            dirs = np.rad2deg(np.arctan2(s1, c1))
        else:
            efth = spec.spec(data_array=True)
            dirs = spec.mdir(data_array=True)

        inds = efth.argmax(dim="freq")

        return dirs[:, :, inds].values

    def name(self):
        return "dirp"

    def unit(self):
        return ureg.deg

    def standard_name(self):
        return "sea_surface_wave_from_direction"


class FF(WaveParameter):
    """Wind speed"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        raise NotImplementedError

    def name(self):
        return "ff"

    def unit(self):
        return ureg.m * ureg.s**-1

    def standard_name(self):
        return "wind_speed"


class DD(WaveParameter):
    """Wind direction"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        raise NotImplementedError

    def name(self):
        return "dd"

    def unit(self):
        return ureg.deg

    def standard_name(self):
        return "wind_from_direction"
