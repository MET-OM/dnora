from abc import ABC, abstractmethod
import xarray as xr
from typing import List
import numpy as np
from ..bnd.conventions import SpectralConvention
from ..bnd import Boundary
from ..spc import Spectra
from typing import Union
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
        ''

    def _is_boundary(self, spec: xr.Dataset) -> bool:
        return 'dirs' in list(spec.coords)

    # def _format_dataset(self, wave_parameter: xr.Dataset, spec: xr.Dataset, data_var=None) -> xr.Dataset:
    #     if data_var is None:
    #         data_var = list(wave_parameter.data_vars)[0]
    #
    #     # Data variable will be name e.g. 'm0' is it has been calculated using Moment(0)
    #     # Change name of data variable to name of the parameter
    #     wave_parameter = wave_parameter.rename({data_var: self.name()})
    #
    #     # Drop extra variables from 1D spectra Dataset
    #     if 'mdir' in list(wave_parameter.data_vars):
    #         wave_parameter = wave_parameter.drop('mdir')
    #     if 'spr' in list(wave_parameter.data_vars):
    #         wave_parameter = wave_parameter.drop('spr')
    #
    #     # Copy over global attributes from the spectral Dataset
    #     #wave_parameter = wave_parameter.assign_attrs({'name': spec.name, 'source': spec.source})
    #
    #     # Set parameter specific attributes
    #     wave_parameter[self.name()].attrs['unit'] = self.unit()
    #     wave_parameter[self.name()].attrs['standard_name'] = self.standard_name()  # None if not defined
    #
    #     return wave_parameter

class Moment(WaveParameter):
    """Spectral moment from spectra"""
    def __init__(self, moment: float) -> None:
        self._moment = moment
        pass

    def __call__(self, spec: Union[Spectra, Boundary]) -> np.ndarray:
        if isinstance(spec, Boundary):
            dD = 360/len(moment.dirs())
            ds = dD*np.pi/180*spec.ds().sum(dim='dirs')
        else:
            ds = spec.ds()

        ds = (ds.spec*(ds.freq**self._moment)).integrate(coord='freq')

        return ds.values

    def unit(self):
        ''

    def name(self):
        """ 1 returns 'm1', 0.5 returns 'm05' etc."""
        strs = f'{self._moment:.1f}'.split('.')
        if strs[1] == '0':
            return f'm{strs[0]}'
        else:
            return f'm{strs[0]}{strs[1]}'

class PowerMoment(WaveParameter):
    """Spectral moment from spectra meighted by spectra any power"""
    def __init__(self, moment: float, power: float) -> None:
        self._moment = moment
        self._power = power

    def __call__(self, spec: Union[Spectra, Boundary]) -> np.ndarray:
        moment = spec
        if self._is_boundary(moment):
            dD = 360/len(moment.dirs)
            moment = dD*np.pi/180*moment.sum(dim='dirs')

        moment = (moment.spec**self._power*(moment.freq**self._moment)).integrate(coord='freq')
        #breakpoint()
        #moment = self._format_dataset(moment, spec)

        return moment.values

    def unit(self):
        ''

    def name(self):
        """ 1 returns 'm1', 0.5 returns 'm05' etc."""
        strs = f'{self._moment:.1f}'.split('.')
        if strs[1] == '0':
            strs = f'm{strs[0]}'
        else:
            strs = f'm{strs[0]}{strs[1]}'

        strs2 = f'{self._power:.1f}'.split('.')
        if strs2[1] == '0':
            strs += f'p{strs2[0]}'
        else:
            strs += f'p{strs2[0]}{strs2[1]}'

        return strs

class Hs(WaveParameter):
    """Singificant wave height from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 4*Moment(0)(spec)**0.5

    def name(self):
        return 'hs'

    def unit(self):
        return 'm'

    def standard_name(self):
        return 'sea_surface_wave_significant_height'

class Tm01(WaveParameter):
    """Mean wave period from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return Moment(0)(spec)/Moment(1)(spec)

    def name(self):
        return 'tm01'

    def unit(self):
        return 's'

    def standard_name(self):
        return 'sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment'

class Tm_10(WaveParameter):
    """Mean wave period based on inverse moment from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return Moment(-1)(spec)/Moment(0)(spec)

    def name(self):
        return 'tm_10'

    def unit(self):
        return 's'

    def standard_name(self):
        return 'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'

class Tm02(WaveParameter):
    """Mean wave period based on second moment from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return (Moment(0)(spec)/Moment(2)(spec))**(0.5)

    def name(self):
        return 'tm02'

    def unit(self):
        return 's'

    def standard_name(self):
        return 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'

class Tp(WaveParameter):
    """Peak period from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return PowerMoment(0,10)(spec)/PowerMoment(-1,10)(spec)

    def name(self):
        return 'tp'

    def unit(self):
        return 's'

    def standard_name(self):
        return 'sea_surface_wave_period_at_variance_spectral_density_maximum'

class Fm(WaveParameter):
    """Mean freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return Moment(1)(spec)/Moment(0)(spec)

    def name(self):
        return 'fm'

    def unit(self):
        return 'Hz'

class Wm(WaveParameter):
    """Mean angular freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return (Moment(1)(spec)/Moment(0)(spec).m0)*2*np.pi

    def name(self):
        return 'wm'

    def unit(self):
        return 'rad/s'


class Fc(WaveParameter):
    """Characteristic freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return PowerMoment(1,4)(spec)/PowerMoment(0,4)(spec)

    def name(self):
        return 'fc'

    def unit(self):
        return 'Hz'

class Wc(WaveParameter):
    """Characteristic angular freqeuncy from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 2*np.pi*PowerMoment(1,4)(spec)/PowerMoment(0,4)(spec)

    def name(self):
        return 'wc'

    def unit(self):
        return 'rad/s'

class Fp(WaveParameter):
    """Peak frequency from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return PowerMoment(1,10)(spec)/PowerMoment(0,10)(spec)

    def name(self):
        return 'fp'

    def unit(self):
        return 'Hz'

    def standard_name(self):
        return 'sea_surface_wave_frequency_at_variance_spectral_density_maximum'


class Fp(WaveParameter):
    """Peak angular frequency from spectra"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        return 2*np.pi*PowerMoment(1,10)(spec)/PowerMoment(0,10)(spec)

    def name(self):
        return 'wp'

    def unit(self):
        return 'Hz'

    def standard_name(self):
        return 'sea_surface_wave_angular_frequency_at_variance_spectral_density_maximum'

class Dirm(WaveParameter):
    """Mean wave direction"""
    def __call__(self, spec: Union[Spectra, Boundary]):
        spec._set_convention(SpectralConvention.MET)

        if isinstance(spec, Boundary):
            theta = np.deg2rad(spec.mdir(data_array=True))
            dD = 360/len(spec.dirs())
            # Normalizing here so that integration over direction becomes summing
            efth = dD*np.pi/180*spec.spec(data_array=True)

            c1 = ((np.cos(theta)*efth).sum(dim='dirs'))  # Function of frequency
            s1 = ((np.sin(theta)*efth).sum(dim='dirs'))
        else:
            theta = np.deg2rad(spec.mdir(data_array=True))

            c1 = np.cos(theta)*spec.spec(data_array=True)  # Function of frequency
            s1 = np.sin(theta)*spec.spec(data_array=True)

        m0 =Moment(0)(spec)

        a1m = c1.integrate(coord='freq')/m0  # Mean parameters
        b1m = s1.integrate(coord='freq')/m0

        thetam = np.arctan2(b1m,a1m)
        dirm = np.mod(thetam*180/np.pi, 360)

        return dirm.values

    def name(self):
        return 'dirm'

    def unit(self):
        return 'deg'

    def standard_name(self):
        return 'sea_surface_wave_from_direction'

class Sprm(WaveParameter):
    """Mean wave spreading"""

    def __call__(self, spec: Union[Spectra, Boundary]):
        m0 =Moment(0)(spec)

        if self._is_boundary(spec):
            theta = np.deg2rad(spec.dirs.values)
            dD = 360/len(spec.dirs)
            # Normalizing here so that integration over direction becomes summing
            efth = dD*np.pi/180*spec

            c1 = ((np.cos(theta)*efth).sum(dim='dirs'))  # Function of frequency
            s1 = ((np.sin(theta)*efth).sum(dim='dirs'))

            a1m = c1.integrate(coord='freq')/m0  # Mean parameters
            b1m = s1.integrate(coord='freq')/m0

            m1 = np.sqrt(b1m**2 + a1m**2)
            sprm = np.sqrt(2-2*(m1))*180/np.pi

        else:
            sprm = ((spec.spr.values**2*spec).integrate(coord='freq')/m0)**0.5


        #sprm = self._format_dataset(sprm, spec)
        return sprm.values

    def name(self):
        return 'sprm'

    def unit(self):
        return 'deg'

    def standard_name(self):
        return 'sea_surface_wave_directional_spread'
