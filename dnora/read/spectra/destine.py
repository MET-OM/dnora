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
from dnora.read.product_readers import SpectralProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.file_structure import FileStructure
from pathlib import Path
from dnora.process.spectra import RemoveEmpty
import glob
from functools import partial
from dnora.spectra import Spectra


def download_ecmwf_from_destine(start_time, filename: str, end_time=None) -> None:
    """Downloads wave data from DestinE. If no end_time is given, a minimal query is done to get the coordinates of the points"""

    start_time = pd.Timestamp(start_time)

    if end_time is None:
        params = "140221"
        steps = "0"
        end_time = start_time
    else:
        params = "140229/140230/140231"
        # Because the reader is set up to do daily chunks, the start and end time will always be in the same day
        days = [str(l) for l in range(start_time.hour, end_time.hour + 1)]
        steps = "/".join(days)
        # steps = "0/1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23"
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
        "param": params,
        "time": "00",
        "step": steps,
    }
    date_str = start_time.strftime("%Y%m%d")
    request_waves["date"] = date_str

    c.retrieve("destination-earth", request_waves, filename)


def destine_wave_ds_read(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    url: str,
    ## Partial variables from ProductReader
    inds: list[int],
    ## Partial variables in ProductConfiguration
    freq0: float,
    nfreq: int,
    finc: float,
    ndirs: int,
    **kwargs,
):
    name = "ECMWF"
    folder = setup_temp_dir(DnoraDataType.SPECTRA, name)

    temp_file = f"{name}_temp.grib"
    grib_file = f"{folder}/{temp_file}"
    download_ecmwf_from_destine(start_time, grib_file, end_time)

    ds = xr.open_dataset(grib_file, engine="cfgrib", decode_timedelta=True)
    ds = ds.isel(values=inds)
    ii = np.where(
        np.logical_and(
            ds.valid_time.values >= start_time, ds.valid_time.values <= end_time
        )
    )[0]
    ds = ds.isel(step=ii)
    msg.plain("Calculating JONSWAP spectra with given Hs and Tp...")
    fp = 1 / ds.pp1d.values
    m0 = ds.swh.values**2 / 16

    freq = np.array([freq0 * finc**n for n in np.linspace(0, nfreq - 1, nfreq)])
    dD = 360 / ndirs
    dirs = np.linspace(0, 360 - dD, ndirs)
    E = utils.spec.jonswap1d(fp=fp, m0=m0, freq=freq)

    msg.plain("Expanding to cos**2s directional distribution around mean direction...")
    Ed = utils.spec.expand_to_directional_spectrum(E, freq, dirs, dirp=ds.mwd.values)
    obj = Spectra.from_ds(ds, freq=freq, dirs=dirs, time=ds.valid_time.values)
    obj.set_spec(Ed)

    return obj.ds()


class ECMWF(SpectralProductReader):
    product_configuration = ProductConfiguration(
        ds_creator_function=partial(
            destine_wave_ds_read, freq0=0.04118, nfreq=32, finc=1.1, ndirs=36
        ),
        convention=SpectralConvention.MET,
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )

    def post_processing(self):
        return RemoveEmpty()

    def caching_strategy(self) -> CachingStrategy:
        return CachingStrategy.SinglePatch

    def get_coordinates(
        self, grid, start_time, source: DataSource, folder: str, **kwargs
    ) -> dict:
        """We only want to do a minimal download to get the coordinates"""
        folder = setup_temp_dir(
            DnoraDataType.SPECTRA, self.name(), clean_old_files=False
        )
        grib_file = f"{folder}/coordinates_ECMWF_destine.grib"

        if glob.glob(grib_file):
            msg.from_file(grib_file)
        else:
            download_ecmwf_from_destine(start_time, grib_file)
        ds = xr.open_dataset(grib_file, engine="cfgrib", decode_timedelta=True)
        return {"lat": ds.latitude.values, "lon": ds.longitude.values}
