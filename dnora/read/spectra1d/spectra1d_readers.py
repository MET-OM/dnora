from copy import deepcopy
from dnora.spectra import Spectra
import numpy as np

from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora import msg
from dnora.type_manager.data_sources import DataSource
from dnora.read.abstract_readers import PointDataReader
from dnora.utils.spec import jonswap1d


class SpectraTo1D(PointDataReader):
    """Integrates boundary spectra to omnidairectional spectra"""

    def __init__(self, spectra: Spectra) -> None:
        self._boundary = deepcopy(spectra)
        if self._boundary is not None:
            self._boundary.set_convention(SpectralConvention.MET)

    def default_data_source(self) -> DataSource:
        return DataSource.CREATION

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        all_points = {
            "lon": self._boundary.lon(strict=True),
            "lat": self._boundary.lat(strict=True),
            "x": self._boundary.x(strict=True),
            "y": self._boundary.y(strict=True),
        }
        return all_points

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        **kwargs,
    ) -> tuple:
        time = (
            self._boundary.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._boundary.lon(strict=True)
        lat = self._boundary.lat(strict=True)
        x = self._boundary.x(strict=True)
        y = self._boundary.y(strict=True)

        freq = self._boundary.freq()
        theta = np.deg2rad(self._boundary.dirs())

        dD = 360 / len(self._boundary.dirs())

        efth = (
            self._boundary.spec(data_array=True, squeeze=False).sel(
                time=slice(start_time, end_time), inds=inds
            )
            * dD
            * np.pi
            / 180
        )
        ef = efth.sum(dim="dirs")
        eth = efth.integrate(coord="freq")

        b1 = ((np.sin(theta) * efth).sum(dim="dirs")) / ef  # Function of frequency
        a1 = ((np.cos(theta) * efth).sum(dim="dirs")) / ef
        thetam = np.arctan2(b1, a1)
        m1 = np.sqrt(b1**2 + a1**2)
        spr = np.sqrt(2 - 2 * (m1)).values * 180 / np.pi

        mdir = np.mod(thetam.values * 180 / np.pi, 360)
        spec = ef.values

        coord_dict = {
            "lon": lon,
            "lat": lat,
            "x": x,
            "y": y,
            "time": time,
            "freq": freq,
        }
        data_dict = {"spec": spec, "dirm": mdir, "spr": spr}
        meta_dict = deepcopy(self._boundary.ds().attrs)

        return coord_dict, data_dict, meta_dict

    def name(self):
        if self._boundary is None:
            return "EmptyData"
        return self._boundary.name


class WaveSeriesToJONSWAP1D(PointDataReader):
    """Integrates boundary spectra to omnidairectional spectra"""

    def __init__(self, waveseries, freq) -> None:
        self._waveseries = waveseries
        self._freq = freq

    def default_data_source(self) -> DataSource:
        return DataSource.CREATION

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        all_points = {
            "lon": self._waveseries.lon(strict=True),
            "lat": self._waveseries.lat(strict=True),
            "x": self._waveseries.x(strict=True),
            "y": self._waveseries.y(strict=True),
        }
        return all_points

    def __call__(
        self,
        obj_type,
        grid,
        start_time,
        end_time,
        source,
        folder,
        filename,
        inds,
        dnora_class,
        **kwargs,
    ) -> tuple:
        time = (
            self._waveseries.time(data_array=True)
            .sel(time=slice(start_time, end_time))
            .values
        )
        lon = self._waveseries.lon(strict=True)
        lat = self._waveseries.lat(strict=True)
        x = self._waveseries.x(strict=True)
        y = self._waveseries.y(strict=True)
        freq = self._freq
        obj = dnora_class(time=time, lon=lon, lat=lat, x=x, y=y, freq=freq)

        if self._waveseries.tp(strict=True) is None:
            raise ValueError("No peak period defined!")

        if self._waveseries.hs(strict=True) is None:
            raise ValueError("No significant wave height defined!")

        fp = 1 / self._waveseries.tp(squeeze=False)
        m0 = (self._waveseries.hs(squeeze=False) / 4) ** 2

        msg.plain("Calculating JONSWAP spectra with given Tp and Hs...")
        E = jonswap1d(fp=fp, m0=m0, freq=freq)

        obj.set_spec(E)
        obj = obj.sel(inds=inds)
        return obj.ds()

    def name(self):
        readername = "JONSWAP"
        if self._waveseries is not None:
            readername += f"_{self._waveseries.name}"
        return readername
