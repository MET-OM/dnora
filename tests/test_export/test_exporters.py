import dnora as dn
import xarray as xr
import numpy as np
import os
import shutil
import stat


def handle_remove_readonly(func, path, exc_info):
    """Clear the read-only bit and reattempt the removal."""
    os.chmod(path, stat.S_IWRITE)  # Change to writable
    func(path)


def test_ww3_spectral_export_one_file():
    grid = dn.grid.TriGrid(lon=(0, 1, 2), lat=(10, 11, 12), name="TestGrid")
    grid.set_boundary_points(dn.grid.mask.All())
    model = dn.modelrun.Constant(grid, year=2019, month=1, day=1)
    model.import_spectra()
    exporter = dn.export.WW3(model)
    exporter.export_spectra(filename="pytest_ww3_spectra")

    ds = xr.open_dataset("TestGrid_WW3/pytest_ww3_spectra.nc")
    np.testing.assert_array_almost_equal(
        np.full((24, 3), [0, 1, 2]), ds.longitude.values
    )
    np.testing.assert_array_almost_equal(
        np.full((24, 3), [10, 11, 12]), ds.latitude.values
    )
    if os.path.isdir("TestGrid_WW3"):
        try:
            shutil.rmtree("TestGrid_WW3", onerror=handle_remove_readonly)
        except PermissionError:
            pass


def test_ww3_spectral_export_three_files():
    grid = dn.grid.TriGrid(lon=(0, 1, 2), lat=(10, 11, 12), name="TestGrid")
    grid.set_boundary_points(dn.grid.mask.All())
    model = dn.modelrun.Constant(grid, year=2019, month=1, day=1)
    model.import_spectra()
    exporter = dn.export.WW3(model)
    exporter.export_spectra(
        one_file=False,
        filename="pytest_ww3_spectra_#LON0_#LAT0",
    )
    for lon, lat in zip([0, 1, 2], [10, 11, 12]):
        ds = xr.open_dataset(
            f"TestGrid_WW3/pytest_ww3_spectra_{lon:010.7f}_{lat:010.7f}.nc"
        )
        np.testing.assert_array_almost_equal(
            np.full((24, 1), [lon]), ds.longitude.values
        )
        np.testing.assert_array_almost_equal(
            np.full((24, 1), [lat]), ds.latitude.values
        )
    if os.path.isdir("TestGrid_WW3"):
        try:
            shutil.rmtree("TestGrid_WW3", onerror=handle_remove_readonly)
        except PermissionError:
            pass


def test_ww3_spectral_export_squeeze_lonlat():
    grid = dn.grid.TriGrid(lon=(0, 1, 2), lat=(10, 11, 12), name="TestGrid")
    grid.set_boundary_points(dn.grid.mask.All())
    model = dn.modelrun.Constant(grid, year=2019, month=1, day=1)
    model.import_spectra()
    exporter = dn.export.WW3(model)
    exporter.export_spectra(squeeze_lonlat=True, filename="pytest_ww3_spectra")

    ds = xr.open_dataset("TestGrid_WW3/pytest_ww3_spectra.nc")
    np.testing.assert_array_almost_equal([0, 1, 2], ds.longitude.values)
    np.testing.assert_array_almost_equal([10, 11, 12], ds.latitude.values)
    if os.path.isdir("TestGrid_WW3"):
        try:
            shutil.rmtree("TestGrid_WW3", onerror=handle_remove_readonly)
        except PermissionError:
            pass
