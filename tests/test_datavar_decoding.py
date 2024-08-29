from dnora.read.waveseries.waveseries_readers import read_data_vars
from geo_skeletons import PointSkeleton
import geo_parameters as gp
import numpy as np


def test_strings():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar("hs")
    data.set_hs(1)
    data.add_datavar("tp")
    data.set_tp(2)

    data_dict = read_data_vars(
        ["hs", "tp"],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=False,
        decode_cf=False,
    )

    assert set(data_dict.keys()) == {"hs", "tp"}
    np.testing.assert_array_almost_equal(data_dict["hs"], data.hs())
    np.testing.assert_array_almost_equal(data_dict["tp"], data.tp())


def test_gp_uninit():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar("hs")
    data.set_hs(1)
    data.add_datavar("tp")
    data.set_tp(2)

    data_dict = read_data_vars(
        [gp.wave.Hs, gp.wave.Tp],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=False,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hs", "tp"}


def test_gp_init():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar("hsig")
    data.set_hsig(1)
    data.add_datavar("tpeak")
    data.set_tpeak(2)

    data_dict = read_data_vars(
        [gp.wave.Hs("hsig"), gp.wave.Tp("tpeak")],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=False,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hsig", "tpeak"}


def test_gp_init_keep_gp_names():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar("hsig")
    data.set_hsig(1)
    data.add_datavar("tpeak")
    data.set_tpeak(2)

    data_dict = read_data_vars(
        [gp.wave.Hs("hsig"), gp.wave.Tp("tpeak")],
        data.ds(),
        keep_gp_names=True,
        keep_source_names=False,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hs", "tp"}


def test_gp_init_keep_source_names():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar("hsig")
    data.set_hsig(1)
    data.add_datavar("tpeak")
    data.set_tpeak(2)

    data_dict = read_data_vars(
        [gp.wave.Hs("hsig"), gp.wave.Tp("tpeak")],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=True,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hsig", "tpeak"}


def test_gp_init_keep_source_names2():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar(gp.wave.Hs("hsig"))
    data.set_hsig(1)
    data.add_datavar(gp.wave.Tp("tpeak"))
    data.set_tpeak(2)

    data_dict = read_data_vars(
        [gp.wave.Hs("hs"), gp.wave.Tp("tp")],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=True,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hsig", "tpeak"}


def test_gp_init_keep_source_names3():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar(gp.wave.Hs("hsig"))
    data.set_hsig(1)
    data.add_datavar(gp.wave.Tp("tpeak"))
    data.set_tpeak(2)

    data_dict = read_data_vars(
        [gp.wave.Hs("hsdd"), gp.wave.Tp("tpss")],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=True,
        decode_cf=False,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hsig", "tpeak"}


def test_gp_decode_cf():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar(gp.wave.Hs("hs"))
    data.set_hs(1)
    data.add_datavar(gp.wave.Tp("tp"))
    data.set_tp(2)

    data_dict = read_data_vars(
        ["hs", "tp"],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=False,
        decode_cf=True,
    )
    key_names = [p.name for p in data_dict.keys()]
    assert set(key_names) == {"hs", "tp"}


def test_gp_decode_cf_false():
    data = PointSkeleton(lon=0, lat=0)
    data.add_datavar(gp.wave.Hs("hs"))
    data.set_hs(1)
    data.add_datavar(gp.wave.Tp("tp"))
    data.set_tp(2)

    data_dict = read_data_vars(
        ["hs", "tp"],
        data.ds(),
        keep_gp_names=False,
        keep_source_names=False,
        decode_cf=False,
    )
    assert set(data_dict.keys()) == {"hs", "tp"}
