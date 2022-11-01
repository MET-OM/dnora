from dnora.skeletons.skeleton import will_grid_be_spherical_or_cartesian as func
import numpy as np

def test_lon_lat_tuple():
    native_x, native_y, xvec, yvec = func(lon=(0.,1.), lat=(2.,3.), x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_lon_lat_tuple_none_tuple():
    native_x, native_y, xvec, yvec = func(lon=(0.,1.), lat=(2.,3.), x=(None,None), y=(None,None))
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_lon_lat_int_tuple():
    native_x, native_y, xvec, yvec = func(lon=(0,1), lat=(2,3), x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_lon_lat_single_tuple():
    native_x, native_y, xvec, yvec = func(lon=(0.), lat=(2.), x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == np.array([0.]))
    assert np.all(yvec == np.array([2.]))

def test_lon_lat_single_value():
    native_x, native_y, xvec, yvec = func(lon=0., lat=2., x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == np.array([0.]))
    assert np.all(yvec == np.array([2.]))

def test_lon_lat_array():
    lonvec = np.array([0.,1.,2.,3.])
    latvec = np.array([2.,3.,4.,5.])
    native_x, native_y, xvec, yvec = func(lon=lonvec, lat=latvec, x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == lonvec)
    assert np.all(yvec == latvec)

def test_lon_lat_int_array():
    lonvec = np.array([0.,1.,2.,3.])
    latvec = np.array([2.,3.,4.,5.])
    native_x, native_y, xvec, yvec = func(lon=lonvec.astype(int), lat=latvec.astype(int), x=None, y=None)
    assert native_x == 'lon'
    assert native_y == 'lat'
    assert np.all(xvec == lonvec)
    assert np.all(yvec == latvec)

def test_x_y_tuple():
    native_x, native_y, xvec, yvec = func(x=(0.,1.), y=(2.,3.), lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_x_y_tuple_none_tuple():
    native_x, native_y, xvec, yvec = func(x=(0.,1.), y=(2.,3.), lon=(None,None), lat=(None,None))
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_x_y_int_tuple():
    native_x, native_y, xvec, yvec = func(x=(0,1), y=(2,3), lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == np.array([0.,1.]))
    assert np.all(yvec == np.array([2.,3.]))

def test_x_y_single_tuple():
    native_x, native_y, xvec, yvec = func(x=(0.), y=(2.), lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == np.array([0.]))
    assert np.all(yvec == np.array([2.]))

def test_x_y_single_value():
    native_x, native_y, xvec, yvec = func(x=0., y=2., lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == np.array([0.]))
    assert np.all(yvec == np.array([2.]))

def test_x_y_array():
    lonvec = np.array([0.,1.,2.,3.])
    latvec = np.array([2.,3.,4.,5.])
    native_x, native_y, xvec, yvec = func(x=lonvec, y=latvec, lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == lonvec)
    assert np.all(yvec == latvec)

def test_x_y_int_array():
    lonvec = np.array([0.,1.,2.,3.])
    latvec = np.array([2.,3.,4.,5.])
    native_x, native_y, xvec, yvec = func(x=lonvec.astype(int), y=latvec.astype(int), lon=None, lat=None)
    assert native_x == 'x'
    assert native_y == 'y'
    assert np.all(xvec == lonvec)
    assert np.all(yvec == latvec)
