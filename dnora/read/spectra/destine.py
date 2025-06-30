import pandas as pd
from dnora.type_manager.data_sources import DataSource
from dnora.type_manager.dnora_types import DnoraDataType
import xarray as xr
import numpy as np
from dnora import msg
from dnora import utils
from dnora.read.ds_read_functions import setup_temp_dir
from dnora.cacher.caching_strategies import CachingStrategy
from dnora.type_manager.spectral_conventions import SpectralConvention
from dnora.read.abstract_readers import SpectralDataReader
from pathlib import Path
from dnora.process.spectra import RemoveEmpty
import glob


def download_ecmwf_from_destine(start_time, end_time, filename: str) -> str:
    """Downloads 10 m wind data DestinE ClimateDT data portfolio data from the Destine Earth Store System  for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        from polytope.api import Client
    except ImportError as e:
        msg.advice(
            "The polytope package is required to acces these data! Install by e.g. 'python -m pip install polytope-client' and 'conda install cfgrib eccodes=2.41.0'"
        )
        raise e
    c = Client(address="polytope.lumi.apps.dte.destination-earth.eu")

    request_waves = {
        "class": "d1",
        "expver": "0001",
        "dataset": "extremes-dt",
        "stream": "wave",
        "type": "fc",
        "levtype": "sfc",
        # Tm02/Hs/Dirm/Tp/Tm
        # "param": "140221/140229/140230/140231/140232",
        "param": "140229/140230/140231",
        "time": "00",
        "step": "0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23",
    }
    date_str = start_time.strftime("%Y%m%d")
    request_waves["date"] = date_str

    c.retrieve("destination-earth", request_waves, filename)


def download_ecmwf_coordinates_from_destine(start_time, end_time, filename: str) -> str:
    """Downloads 10 m wind data DestinE ClimateDT data portfolio data from the Destine Earth Store System  for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        from polytope.api import Client
    except ImportError as e:
        msg.advice(
            "The polytope package is required to acces these data! Install by e.g. 'python -m pip install polytope-client' and 'conda install cfgrib eccodes=2.41.0'"
        )
        raise e
    c = Client(address="polytope.lumi.apps.dte.destination-earth.eu")

    request_waves = {
        "class": "d1",
        "expver": "0001",
        "dataset": "extremes-dt",
        "stream": "wave",
        "type": "fc",
        "levtype": "sfc",
        # Tm02/Hs/Dirm/Tp/Tm
        "param": "140221",
        # "param": "140229/140229",
        "time": "00",
        "step": "0",
    }
    date_str = start_time.strftime("%Y%m%d")
    request_waves["date"] = date_str

    c.retrieve("destination-earth", request_waves, filename)


class ECMWF(SpectralDataReader):
    def convention(self) -> str:
        return SpectralConvention.MET

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def post_processing(self):
        return RemoveEmpty()

    def caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.SinglePatch

    def get_coordinates(
        self, grid, start_time, source: DataSource, folder: str, **kwargs
    ) -> dict:
        """The download of the file take place already here so we don't have to do it twice"""
        if not folder:
            folder = setup_temp_dir(
                DnoraDataType.SPECTRA, self.name(), clean_old_files=False
            )
        grib_file = f"{folder}/coordinates_ECMWF_destine.grib"

        if glob.glob(grib_file):
            msg.from_file(grib_file)
        else:
            download_ecmwf_coordinates_from_destine(start_time, start_time, grib_file)
        ds = xr.open_dataset(grib_file, engine="cfgrib", decode_timedelta=True)
        return {"lat": ds.latitude.values, "lon": ds.longitude.values}

    def __call__(
        self,
        obj_type,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        dnora_class=None,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        msg.info(
            f"Getting Destine boundary spectra using JONSWAP fits from {start_time} to {end_time}"
        )

        if not folder:
            folder = setup_temp_dir(obj_type, self.name(), clean_old_files=not filename)
        temp_file = filename or f"{self.name()}_temp.grib"
        grib_file = f"{folder}/{temp_file}"

        # If a filename is not given, then call the API to download data
        if filename is None:
            download_ecmwf_from_destine(start_time, end_time, grib_file)
        else:
            msg.from_file(grib_file)

        ds = xr.open_dataset(grib_file, engine="cfgrib", decode_timedelta=True)
        ds = ds.isel(values=inds)

        msg.plain("Calculating JONSWAP spectra with given Hs and Tp...")
        fp = 1 / ds.pp1d.values
        m0 = ds.swh.values**2 / 16
        freq0: float = 0.04118
        nfreq: int = 32
        finc: float = 1.1
        freq = np.array([freq0 * finc**n for n in np.linspace(0, nfreq - 1, nfreq)])
        dirs = np.linspace(0, 350, 36)
        E = utils.spec.jonswap1d(fp=fp, m0=m0, freq=freq)

        msg.plain(
            "Expanding to cos**2s directinal distribution around mean direction..."
        )
        Ed = utils.spec.expand_to_directional_spectrum(
            E, freq, dirs, dirp=ds.mwd.values
        )
        obj = dnora_class.from_ds(ds, freq=freq, dirs=dirs, time=ds.valid_time.values)
        obj.set_spec(Ed)

        return obj.ds()
