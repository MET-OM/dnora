from dnora import modelrun, pick
from dnora.read.generic import ConstantData
import numpy as np
import dnora as dn
import pandas as pd


def test_import_constant_waveseries_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = modelrun.ModelRun(
        grid, start_time="2020-01-01 00:00", end_time="2020-01-02 00:00"
    )

    model.import_spectra(ConstantData(), point_picker=pick.Trivial())
    model.spectra_to_waveseries()

    assert model.waveseries().size() == (25, 1)
    manual_hs = 4 * np.sqrt(
        np.sum(model.spectra().isel(time=1).spec())
        * model.spectra().dd(angular=True)
        * model.spectra().df()
    )
    np.testing.assert_almost_equal(np.mean(model.waveseries().hs()), manual_hs)
