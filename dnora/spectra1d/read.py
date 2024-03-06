from copy import copy
from dnora.spectra import Spectra
import numpy as np

from dnora.spectral_conventions import convert_2d_to_1d

from dnora.dnora_type_manager.data_sources import DataSource
from dnora.readers.abstract_readers import SpectralDataReader


class SpectraTo1D(SpectralDataReader):
    """Integrates boundary spectra to omnidairectional spectra"""

    def __init__(self, spectra: Spectra) -> None:
        self._boundary = copy(spectra)
        # self._boundary._set_convention(SpectralConvention.OCEAN)

    def convention(self):
        return convert_2d_to_1d(self._boundary._convention)

    def get_coordinates(
        self, grid, start_time: str, source: DataSource, folder: str, **kwargs
    ) -> dict[str : np.ndarray]:
        all_points = {
            "lon": self._boundary.lon(strict=True),
            "lat": self._boundary.lat(strict=True),
            "x": self._boundary.x(strict=True),
            "y": self._boundary.y(strict=True),
        }
        return all_points

    def __call__(
        self, grid, start_time, end_time, inds, source: DataSource, **kwargs
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
            self._boundary.spec(data_array=True).sel(
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
        meta_dict = self._boundary.ds().attrs
        metaparameter_dict = {}
        return coord_dict, data_dict, meta_dict, metaparameter_dict

    def name(self):
        if self._boundary is None:
            return "EmptyData"
        return self._boundary.name
