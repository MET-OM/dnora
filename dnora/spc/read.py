from abc import ABC, abstractmethod
from typing import Tuple
from copy import copy
from ..bnd.bnd_mod import Boundary
import numpy as np

from ..bnd.conventions import SpectralConvention
from ..bnd.conventions import convert_2d_to_1d


class SpectralReader(ABC):
    """Reads boundary spectra from some source and provide it to the object."""
    def __init__(self):
        pass

    @abstractmethod
    def get_coordinates(self, grid, start_time):
        """Return a list of all the available coordinated in the source.

        These are needed fo the PointPicker object to choose the relevant
        point to actually read in.

        Provide the result as two equally long nump arrays.
        """
        return lon, lat

    @abstractmethod
    def convention(self) -> SpectralConvention:
        """Return the convention of the spectra returned to the object.

        The conventions to choose from are predetermined:

        OCEAN:    Oceanic convention
                    Direction to. North = 0, East = 90.

        MET:      Meteorological convention
                    Direction from. North = 0, East = 90.

        MATH:     Mathematical convention
                    Direction to. North = 90, East = 0.
        """
        return convention

    @abstractmethod
    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        """Reads in the spectra from inds between start_time and end_time.

        The variables needed to be returned are:

        time:   Time stamps as numpy.datetime64 array
        freq:   Frequency vector as numpy array
        spec:   Sectra [time, station, freq] as numpy array
        dirm:   Mean direction as numpy array [None if not calculated]
        spr:    Spreading as a numpy array [None if not calcultated]
        lon:    Longitude vector as numpy array
        lat:    Latitude vector as numpy array
        source: Source of the data as String
        """

        return time, freq, spec, mdir, spr, lon, lat, x, y, source

class BoundaryToSpectra(SpectralReader):
    """Integrates boundary spectra to omnidairectional spectra"""
    def __init__(self, boundary: Boundary) -> None:
        self._boundary = copy(boundary)

    def convention(self):
        return convert_2d_to_1d(self._boundary._convention)

    def get_coordinates(self, grid, start_time: str) -> Tuple[np.ndarray, np.ndarray]:
        return self._boundary.lon(), self._boundary.lat()
        #return self._boundary.data.lon.values, self._boundary.data.lat.values

    def __call__(self, grid, start_time, end_time, inds) -> Tuple:
        self.name = self._boundary.name
        #source = self._boundary.data.source
        time = self._boundary.time(data_array=True).sel(time=slice(start_time, end_time)).values
        lon = self._boundary.lon(strict=True)
        lat = self._boundary.lat(strict=True)
        x = self._boundary.x(strict=True)
        y = self._boundary.y(strict=True)

        freq = self._boundary.freq()
        theta = np.deg2rad(self._boundary.dirs())
        dD = 360/len(self._boundary.dirs())
        # Normalizing here so that integration over direction becomes summing
        #efth = self._boundary.data.sel(time=slice(start_time, end_time), x=inds).spec*dD*np.pi/180
        efth = self._boundary.spec(data_array=True).sel(time=slice(start_time, end_time), inds=inds)*dD*np.pi/180
        ef = efth.sum(dim='dirs')
        eth = efth.integrate(coord='freq')
        # m0 = ef.integrate(coord='freq')

        # b1 = np.sin(theta)*efth  # Function of theta and frequency
        # a1 = np.cos(theta)*efth
        # thetam = np.arctan2(b1.sum(dim='dirs'),a1.sum(dim='dirs')) # Function of frequency
        # spr = np.sqrt(2*np.sin(0.5*(thetam-theta))**2*eth).sum(dim='dirs').values*180/np.pi

        b1 = ((np.sin(theta)*efth).sum(dim='dirs'))/ef  # Function of frequency
        a1 = ((np.cos(theta)*efth).sum(dim='dirs'))/ef
        thetam = np.arctan2(b1,a1)
        m1 = np.sqrt(b1**2+a1**2)
        spr = np.sqrt(2-2*(m1)).values*180/np.pi

        mdir = np.mod(thetam.values*180/np.pi, 360)
        spec = ef.values

        return time, freq, spec, mdir, spr, lon, lat, x, y, self._boundary.ds().attrs
