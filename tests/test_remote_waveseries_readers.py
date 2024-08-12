import dnora as dn
import pytest
import pandas as pd
import numpy as np
from datetime import datetime, timedelta


@pytest.mark.remote
def test_e39():
    model = dn.modelrun.ModelRun(
        start_time="2018-01-01 00:00", end_time="2018-05-01 00:00"
    )
    model.import_waveseries(dn.waveseries.read.E39(loc="D"))
    assert model.waveseries().time()[0] == pd.Timestamp("2018-01-01 00:00")
    assert model.waveseries().time()[-1] == pd.Timestamp("2018-05-01 00:00")
