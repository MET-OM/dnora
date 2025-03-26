import dnora as dn
import glob
import os
import shutil


def test_write_spectral_cache():
    grid = dn.grid.Grid(lon=(4, 11), lat=(60, 65))
    grid.set_spacing(dlon=1, dlat=1)

    model = dn.modelrun.ModelRun(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )
    if os.path.isdir("spectra_cache"):
        shutil.rmtree("spectra_cache")

    model.import_spectra(
        dn.read.generic.ConstantData(debug_cache=True), write_cache=True
    )
    assert glob.glob("spectra_cache/constantdata/*") == [
        "spectra_cache/constantdata/2020"
    ]
    assert set(glob.glob("spectra_cache/constantdata/2020/*")) == {
        "spectra_cache/constantdata/2020/01",
        "spectra_cache/constantdata/2020/02",
    }
    jan_files = glob.glob("spectra_cache/constantdata/2020/01/*")
    days = 1
    # Area: lon: (4.00, 11.00), (60.00, 65.00)
    # >>> Using expansion_factor = 1.50 <<<
    # Area: lon: (2.25, 12.75), (58.75, 66.25)
    lat_tiles = 3  # 60-65, 65-70 + 55-60 from the expansion!
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(jan_files) == days * lat_tiles * lon_tiles

    feb_files = glob.glob("spectra_cache/constantdata/2020/02/*")
    days = 2
    lat_tiles = 3  # 55-60, 60-65, 65-70
    lon_tiles = 3  # 0-5, 5-10, 10-15
    assert len(feb_files) == days * lat_tiles * lon_tiles

    if os.path.isdir("spectra_cache"):
        shutil.rmtree("spectra_cache")