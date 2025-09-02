import dnora as dn
import glob
import os
import shutil
import pytest
import numpy as np
from dnora.type_manager.spectral_conventions import SpectralConvention
import stat


def handle_remove_readonly(func, path, exc_info):
    """Clear the read-only bit and reattempt the removal."""
    os.chmod(path, stat.S_IWRITE)  # Change to writable
    func(path)


def cleanup():
    if os.path.isdir("ww3_pytest"):
        shutil.rmtree("ww3_pytest", onerror=handle_remove_readonly)


@pytest.fixture(scope="session")
def grid():
    grid = dn.grid.Grid(lon=(10, 11), lat=(70, 71))
    grid.set_spacing(dlon=1, dlat=1)
    grid.set_boundary_points(dn.grid.mask.All())
    return grid


@pytest.fixture(scope="session")
def model(grid):
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )
    return model


def check_that_spectra_is_consistent(spectra, convention):
    assert spectra.isel(inds=0, time=0).spec().shape == (10, 36)
    assert spectra.freq().shape == (10,)
    assert spectra.dirs().shape == (36,)
    assert spectra.convention() == convention


def check_that_peak_dir_is_right(spectra, expected_peak_dir: float):
    one_1dspec = spectra.isel(inds=0, time=0, freq=2)
    dirs = one_1dspec.dirs()
    ipeak = np.argmax(one_1dspec.spec())
    np.testing.assert_almost_equal(dirs[ipeak], expected_peak_dir)


def dirs_are_normal(spectra):
    np.testing.assert_almost_equal(spectra.dirs(), np.arange(0, 360, 10))


def dirs_are_math(spectra):
    np.testing.assert_almost_equal(
        spectra.dirs(), np.mod(np.arange(360, 0, -10) - 270, 360)
    )


def test_ww3_writes_ww3(model):
    cleanup()
    model.import_spectra()
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    exp = dn.export.WW3(model)
    dirs_are_normal(model.spectra())

    # Export to file and check that the spectra is converted right
    exp.export_spectra(
        folder="ww3_pytest", filename="spec_in_ww3"
    )  # Default exports (i.e. converts) to WW3 format
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.WW3)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    dirs_are_math(model.spectra())

    # Reread spectra to see that it is correctly written
    model["spectra"] = None
    model.import_spectra(
        dn.read.generic.PointNetcdf(),
        folder="ww3_pytest",
        filename="spec_in_ww3.nc",
        convention="ww3",
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.WW3)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    dirs_are_math(model.spectra())
    cleanup()


def test_ww3_writes_ocean(model):
    cleanup()
    model.import_spectra()
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    exp = dn.export.WW3(model)
    dirs_are_normal(model.spectra())

    # Export to file and check that the spectra is converted right
    # exp.export_spectra(dn.export.spectra_writers.WW3(convention='ocean'),folder='ww3_pytest', filename='spec_in_ocean')  # Equivalent to the line below
    exp.export_spectra(
        folder="ww3_pytest", filename="spec_in_ocean", convention="ocean"
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    dirs_are_normal(model.spectra())

    # Reread spectra to see that it is correctly written
    model["spectra"] = None
    model.import_spectra(
        dn.read.generic.PointNetcdf(),
        folder="ww3_pytest",
        filename="spec_in_ocean.nc",
        convention="ocean",
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    dirs_are_normal(model.spectra())
    cleanup()


def test_ww3_writes_met(model):
    cleanup()
    model.import_spectra()
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    exp = dn.export.WW3(model)
    dirs_are_normal(model.spectra())

    # Export to file and check that the spectra is converted right
    exp.export_spectra(folder="ww3_pytest", filename="spec_in_met", convention="met")
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MET)
    check_that_peak_dir_is_right(model.spectra(), 180.0)
    dirs_are_normal(model.spectra())

    # Reread spectra to see that it is correctly written
    model["spectra"] = None
    model.import_spectra(
        dn.read.generic.PointNetcdf(),
        folder="ww3_pytest",
        filename="spec_in_met.nc",
        convention="met",
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MET)
    check_that_peak_dir_is_right(model.spectra(), 180.0)
    dirs_are_normal(model.spectra())
    cleanup()


def test_ww3_writes_math(model):
    cleanup()
    model.import_spectra()
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    exp = dn.export.WW3(model)
    dirs_are_normal(model.spectra())

    # Export to file and check that the spectra is converted right
    exp.export_spectra(folder="ww3_pytest", filename="spec_in_math", convention="math")
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MATH)
    check_that_peak_dir_is_right(model.spectra(), 90.0)
    dirs_are_normal(model.spectra())

    # Reread spectra to see that it is correctly written
    model["spectra"] = None
    model.import_spectra(
        dn.read.generic.PointNetcdf(),
        folder="ww3_pytest",
        filename="spec_in_math.nc",
        convention="math",
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MATH)
    check_that_peak_dir_is_right(model.spectra(), 90.0)
    dirs_are_normal(model.spectra())
    cleanup()


def test_ww3_writes_mathvec(model):
    cleanup()
    model.import_spectra()
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.OCEAN)
    check_that_peak_dir_is_right(model.spectra(), 0.0)
    exp = dn.export.WW3(model)
    dirs_are_normal(model.spectra())

    # Export to file and check that the spectra is converted right
    exp.export_spectra(
        folder="ww3_pytest", filename="spec_in_mathvec", convention="mathvec"
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MATHVEC)
    check_that_peak_dir_is_right(model.spectra(), 90.0)
    dirs_are_math(model.spectra())

    # Reread spectra to see that it is correctly written
    model["spectra"] = None
    model.import_spectra(
        dn.read.generic.PointNetcdf(),
        folder="ww3_pytest",
        filename="spec_in_mathvec.nc",
        convention="mathvec",
    )
    check_that_spectra_is_consistent(model.spectra(), SpectralConvention.MATHVEC)
    check_that_peak_dir_is_right(model.spectra(), 90.0)
    dirs_are_math(model.spectra())
    cleanup()
