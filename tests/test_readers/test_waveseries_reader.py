import dnora as dn
import numpy as np


def test_spectra_to_waveseries():
    model = dn.modelrun.Constant(
        start_time="2018-01-01 00:00", end_time="2018-05-01 00:00"
    )
    model.import_spectra()
    model.spectra_to_waveseries()
    np.testing.assert_almost_equal(model.waveseries().tp()[0], 1 / 0.3)
    np.testing.assert_almost_equal(
        model.waveseries().hs()[0], 4 * (1 * 0.1 * np.deg2rad(10)) ** 0.5
    )
    np.testing.assert_almost_equal(model.waveseries().dirp(dir_type="to")[0], 0)
