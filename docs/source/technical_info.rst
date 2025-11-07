Technical information
====================================

Here you can find technical information that is probably most useful if you are doing some kind of development of dnora, but it also serves as a reference manual for available parameters etc.


List of data types
+++++++++++++++++++

dnora has defined classes for each type of data. E.g. the class for representing wind data is ``dnora.wind.Wind``, and the corresponding methods are called ``import_wind``, ``export_wind`` etc. 

All classes are built using the geo-skeletons package (https://pypi.org/project/geo-skeletons/) to ensure consistency. This means that gridded data inherits from ``geo_skeletons.GriddedSkeleton`` and non-gridded data inherits from ``geo_skeletons.PointSkeleton``. The variables are represented gy geo-parameter (https://pypi.org/project/geo-parameters/) classes containin meta data that can be used and defines e.g. directional conventions.

The main advantages of have these pre-defined classes is that you have all the metadata etc. in the class. That means that you can inspect the class using the class ``core`` property. If you use a data object (such as ``xarray.Dataset``), then you need to first set the data, and the data itself defines the structure, metadata etc.

The actual data of an instance that are initialized from the objects are stored in an ``xarray.Dataset`` that is accessable by the  ``.ds()`` method. This dataset is created based on the specifications that are imposed on the class.

The available data types are:


Wind 
-----

Wind on one height level (typically 10m, but strictly speaking not specified) and is the main forcing object for wave models. It stores the u and v components of the wind, but contains wrappers for accessing the magnitude and direction.

.. code-block:: python

    >>> from dnora.wind import Wind
    >>> Wind
    <class 'dnora.wind.wind.Wind'>
    >>> Wind.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        u  (time, y, x):  0.1 [m/s] x_wind
        v  (time, y, x):  0.1 [m/s] y_wind
    Masks:
        *empty*
    Magnitudes:
    mag: magnitude of (u,v) [m/s] wind_speed
    Directions:
    dir: direction of (u,v) [deg] wind_from_direction
    --------------------------------------------------------------------------------


Spectra 
-----------------------------------

Directional (2D) spectra that are used for boundary conditions in spectral wave models. The spectral conventions differs between models, but the handeling of this is mostly automatic in dnora. For more information, see the section on spectral conventions.

.. code-block:: python

    >>> from dnora.spectra import Spectra
    >>> Spectra
    <class 'dnora.spectra.spectra.Spectra'>
    >>> Spectra.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (inds)
    Grid:       (time, inds)
    Gridpoint:  (freq, dirs)
    All:        (time, inds, freq, dirs)
    ------------------------------------- Data -------------------------------------
    Variables:
        spec  (time, inds, freq, dirs):  0.0 [m**2/Hz/rad] sea_surface_wave_directional_variance_spectral_density
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------

Spectra1D
------------------------------------------

Omni-directional (1D) spectral data, and is therefore an integrated version of the Spectra class. Note, that you can use the ``'dir_type'`` keyword to get the mean direction in some other convention (see the ``Wind`` class for an explanation).

.. code-block:: python

    >>> from dnora.spectra1d import Spectra1D
    >>> Spectra1D
    <class 'dnora.spectra1d.spectra1d.Spectra1D'>
    >>> Spectra1D.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (inds)
    Grid:       (time, inds)
    Gridpoint:  (freq)
    All:        (time, inds, freq)
    ------------------------------------- Data -------------------------------------
    Variables:
        spr   (time, inds, freq):  0.0 [deg] sea_surface_wave_directional_spread
        dirm  (time, inds, freq):  0.0 [deg] sea_surface_wave_from_direction
        spec  (time, inds, freq):  0.0 [m**2/Hz] sea_surface_wave_variance_spectral_density
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------

WaveSeries
---------------------------------------------

Integrated wave paramers in a non-gridded format. Some on the parameters (wind, current) are strictly not wave parameters, but are included since they are typically present in wave model output files, and it is conventient to have them in the same class. For a gridded representation of wave parameters, see the ``WaveGrid`` class.

For an explanation of how to get different from/to directions for the directional parameters, see the ``Wind`` class.

.. code-block:: python

    >>> from dnora.waveseries import WaveSeries
    >>> WaveSeries
    <class 'dnora.waveseries.waveseries.WaveSeries'>
    >>> WaveSeries.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (inds)
    Grid:       (inds)
    Gridpoint:  (time)
    All:        (time, inds)
    ------------------------------------- Data -------------------------------------
    Variables:
        hs           (time, inds):  0.0 [m] sea_surface_wave_significant_height
        tp           (time, inds):  0.0 [s] sea_surface_wave_period_at_variance_spectral_density_maximum
        dirp         (time, inds):  0.0 [deg] sea_surface_wave_from_direction_at_variance_spectral_density_maximum
        dirm         (time, inds):  0.0 [deg] sea_surface_wave_from_direction
        tm01         (time, inds):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10        (time, inds):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02         (time, inds):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        hs_swell     (time, inds):  0.0 [m] sea_surface_swell_wave_significant_height
        tp_swell     (time, inds):  0.0 [s] sea_surface_swell_wave_period_at_variance_spectral_density_maximum
        dirp_swell   (time, inds):  0.0 [deg] sea_surface_swell_wave_from_direction_at_variance_spectral_density_maximum
        dirm_swell   (time, inds):  0.0 [deg] sea_surface_swell_wave_from_direction
        tm01_swell   (time, inds):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10_swell  (time, inds):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02_swell   (time, inds):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        hs_sea       (time, inds):  0.0 [m] sea_surface_wind_wave_significant_height
        tp_sea       (time, inds):  0.0 [s] sea_surface_wind_wave_period_at_variance_spectral_density_maximum
        dirp_sea     (time, inds):  0.0 [deg] sea_surface_wind_wave_from_direction_at_variance_spectral_density_maximum
        dirm_sea     (time, inds):  0.0 [deg] sea_surface_wind_wave_from_direction
        tm01_sea     (time, inds):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10_sea    (time, inds):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02_sea     (time, inds):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        depth        (time, inds):  0.0 [m] sea_floor_depth_below_sea_surface
        x_wind       (time, inds):  0.1 [m/s] x_wind
        y_wind       (time, inds):  0.1 [m/s] y_wind
        x_current    (time, inds):  0.1 [m/s] sea_water_x_velocity
        y_current    (time, inds):  0.1 [m/s] sea_water_y_velocity
    Masks:
        buoy_mask  (inds):  True
    Magnitudes:
    ff: magnitude of (x_wind,y_wind) [m/s] wind_speed
    current: magnitude of (x_current,y_current) [m/s] sea_water_speed
    Directions:
    dd: direction of (x_wind,y_wind) [deg] wind_from_direction
    current_dir: direction of (x_current,y_current) [deg] sea_water_velocity_to_direction
    --------------------------------------------------------------------------------

Current
---------------------------------------------

Surface current and is used to force wave model runs. It stores the u and v components of the wind, but contains wrappers for accessing the magnitude and direction.  

For an explanation of how to get different from/to directions for the directional parameters, see the ``Wind`` class.

.. code-block:: python

    >>> from dnora.current import Current
    >>> Current
    <class 'dnora.current.current.Current'>
    >>> Current.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        u  (time, y, x):  0.1 [m/s] sea_water_x_velocity
        v  (time, y, x):  0.1 [m/s] sea_water_y_velocity
    Masks:
        *empty*
    Magnitudes:
    mag: magnitude of (u,v) [m/s] sea_water_speed
    Directions:
    dir: direction of (u,v) [deg] sea_water_velocity_to_direction
    --------------------------------------------------------------------------------


WaterLevel
---------------------------------------------

Sea surface elevation that can be used as a forcing for some wave models. Nota, that is can also represent the output from a phase resolving wave model.

.. code-block:: python

    >>> from dnora.waterlevel import WaterLevel
    >>> WaterLevel
    <class 'dnora.waterlevel.waterlevel.WaterLevel'>
    >>> WaterLevel.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        eta  (time, y, x):  0.0 [m] sea_surface_elevation
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------

Ice
---------------------------------------------

Sea ice thickness and concentration, which can be used to force some wave models.

.. code-block:: python

    >>> from dnora.ice import Ice
    >>> Ice
    <class 'dnora.ice.ice.Ice'>
    >>> Ice.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        sit  (time, y, x):  0.0 [m] sea_ice_thickness
        sic  (time, y, x):  0.0 [%] sea_ice_area_fraction
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------

Ocean
---------------------------------------------

Sea surface temperature and salinity.

.. code-block:: python

    >>> from dnora.ocean import Ocean
    >>> Ocean
    <class 'dnora.ocean.ocean.Ocean'>
    >>> Ocean.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        sst  (time, y, x):  nan [K] sea_surface_temperature
        sss  (time, y, x):  nan [g/kg] sea_surface_salinity
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------


WaveGrid
---------------------------------------------

Integrated parameters in a gridded format. Essentially the gridded version of the ``WaveSeries`` class. Can also be used to represent spectral model output data.

.. code-block:: python

    >>> from dnora.wavegrid import WaveGrid
    >>> WaveGrid
    <class 'dnora.wavegrid.wavegrid.WaveGrid'>
    >>> WaveGrid.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (y, x)
    Gridpoint:  (time)
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        hs           (time, y, x):  0.0 [m] sea_surface_wave_significant_height
        tp           (time, y, x):  0.0 [s] sea_surface_wave_period_at_variance_spectral_density_maximum
        dirp         (time, y, x):  0.0 [deg] sea_surface_wave_from_direction_at_variance_spectral_density_maximum
        dirm         (time, y, x):  0.0 [deg] sea_surface_wave_from_direction
        tm01         (time, y, x):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10        (time, y, x):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02         (time, y, x):  0.0 [s] sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        hs_swell     (time, y, x):  0.0 [m] sea_surface_swell_wave_significant_height
        tp_swell     (time, y, x):  0.0 [s] sea_surface_swell_wave_period_at_variance_spectral_density_maximum
        dirp_swell   (time, y, x):  0.0 [deg] sea_surface_swell_wave_from_direction_at_variance_spectral_density_maximum
        dirm_swell   (time, y, x):  0.0 [deg] sea_surface_swell_wave_from_direction
        tm01_swell   (time, y, x):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10_swell  (time, y, x):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02_swell   (time, y, x):  0.0 [s] sea_surface_swell_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        hs_sea       (time, y, x):  0.0 [m] sea_surface_wind_wave_significant_height
        tp_sea       (time, y, x):  0.0 [s] sea_surface_wind_wave_period_at_variance_spectral_density_maximum
        dirp_sea     (time, y, x):  0.0 [deg] sea_surface_wind_wave_from_direction_at_variance_spectral_density_maximum
        dirm_sea     (time, y, x):  0.0 [deg] sea_surface_wind_wave_from_direction
        tm01_sea     (time, y, x):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_first_frequency_moment
        tm_10_sea    (time, y, x):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment
        tm02_sea     (time, y, x):  0.0 [s] sea_surface_wind_wave_mean_period_from_variance_spectral_density_second_frequency_moment
        depth        (time, y, x):  0.0 [m] sea_floor_depth_below_sea_surface
        x_wind       (time, y, x):  0.1 [m/s] x_wind
        y_wind       (time, y, x):  0.1 [m/s] y_wind
        x_current    (time, y, x):  0.1 [m/s] sea_water_x_velocity
        y_current    (time, y, x):  0.1 [m/s] sea_water_y_velocity
    Masks:
        buoy_mask  (y, x):  True
    Magnitudes:
    ff: magnitude of (x_wind,y_wind) [m/s] wind_speed
    current: magnitude of (x_current,y_current) [m/s] sea_water_speed
    Directions:
    dd: direction of (x_wind,y_wind) [deg] wind_from_direction
    current_dir: direction of (x_current,y_current) [deg] sea_water_velocity_to_direction
    --------------------------------------------------------------------------------

Atmoshphere
---------------------

Basic atmospheric variables: Temperature at 2m height, relative humidity and air pressure at mean sea level. Used for the VesselIcing model.

.. code-block:: python

    >>> from dnora.atmosphere import Atmosphere
    >>> Atmosphere
    <class 'dnora.atmosphere.atmosphere.Atmosphere'>
    >>> Atmosphere.core
    ------------------------------ Coordinate groups -------------------------------
    Spatial:    (y, x)
    Grid:       (time, y, x)
    Gridpoint:  *empty*
    All:        (time, y, x)
    ------------------------------------- Data -------------------------------------
    Variables:
        t2m   (time, y, x):  nan [K] air_temperature
        r     (time, y, x):  nan [%] relative_humidity
        mslp  (time, y, x):  nan [N/m**2] air_pressure_at_mean_sea_level
    Masks:
        *empty*
    Magnitudes:
        *empty*
    Directions:
        *empty*
    --------------------------------------------------------------------------------


Directional variables and conventions
++++++++++++++++++++++++++++++++++++++

Directional variables can be added either directly (liek wave direction) or through a wrapper that calculates then on the fly from the stored components. Both have a default direction specified in the standard_name, but can compute another directional type when requested. This can be useful when needing data in a specific format when e.g. exporting it.

We take the WaveSeries object as an example, since it contains both direct variables and wrappers. 

.. code-block:: python

    >>> data = WaveSeries(lon=0, lat=0, time=('2020-01-01 00:00', '2020-01-01 03:00'))
    >>> data.set_dirm(180) # Stored directly, direction from
    >>> data.set_x_wind(0) # components fro ff and dd (direction from)
    >>> data.set_y_wind(10)
    >>> data.set_current_dir(190) # Stored as components, but can set direction directly also (direction to)
    >>> data.set_current(0.2) # Set magnitude that is stored as components

    >>> data.dirm()
    array([180., 180., 180., 180.])
    >>> data.dirm(dir_type='from')
    array([180., 180., 180., 180.])
    >>> data.dirm(dir_type='to')
    array([0., 0., 0., 0.])

    >>> data.dd()
    array([180., 180., 180., 180.])
    >>> data.dd(dir_type='from')
    array([180., 180., 180., 180.])
    >>> data.dd(dir_type='to')
    array([0., 0., 0., 0.])

    >>> data.current_dir()
    array([190., 190., 190., 190.])
    >>> data.current_dir(dir_type='from')
    array([10., 10., 10., 10.])
    >>> data.current_dir(dir_type='to')
    array([190., 190., 190., 190.])

    >>> data.x_current()
    array([-0.03472964, -0.03472964, -0.03472964, -0.03472964])
    >>> data.x_wind()
    array([0, 0, 0, 0])

As you can see, you don't need to know if the class stores components or not.

The directional type is in the standard name, but you can also check it directly:

.. code-block:: python

    >>> data.core.meta_parameter('dd')
    <class 'geo_parameters.wind.WindDir'>
    >>> data.core.meta_parameter('dd').standard_name()
    'wind_from_direction'
    >>> data.core.meta_parameter('dd').dir_type()
    'from'

Using spectral conventions
++++++++++++++++++++++++++++

In practice, you just need to define which spectral convention you are reading in readers, and which one you want when you export. Dnora will take care of the transformations.

One example is the class that reads SWAN netcdf spectra:

.. code-block:: python

    class SWAN_Nc(PointNetcdf):
        def convention(self) -> str:
            return SpectralConvention.MET

        def __call__(self, *args, **kwargs):
            ds = super().__call__(*args, **kwargs)
            theta = ds.dirs.values
            dd = np.mod(np.rad2deg(theta) + 360, 360)
            ds["dirs"] = np.round(dd, 2)
            ds = ds.sortby("dirs")
            return ds

Here the ``convention`` method tells dnora that the spectral data returned by this reader is in meteorological convention. If you have more complex reader that inherits from ``ProductReader`` to take care of parising files and handling multiple data sources, it might look like this:

.. code-block:: python

    class NORA3(SpectralProductReader):
        product_configuration = ProductConfiguration(
            filename="SPC%Y%m%d00.nc",
            default_folders={
                DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m",
                DataSource.INTERNAL: "WINDSURFER/mw3hindcast/spectra/%Y/%m",
            },
            ds_creator_function=partial(basic_xarray_read, inds_var="x"),
            convention=SpectralConvention.OCEAN,
            default_data_source=DataSource.REMOTE,
            ds_aliases={"SPEC": gp.wave.Efth},
        )

        file_structure = FileStructure(
            stride=24,
            hours_per_file=24,
        )


To make sure that we get the right type of data when exporting, we just make a similar trick in the exporer class. E.g. the SWAN writer class starts like:

.. code-block:: python

    class SWAN(SpectraWriter):
        def convention(self) -> SpectralConvention:
            """Convention of spectra"""
            return SpectralConvention.MET


When the conventions have been defined in the reading and writing stages and any combination of readers and writers can be used. Dnora automatically transforms the spectra under the hood without the user needing to know what convention is used where.


List of spectral conventions
++++++++++++++++++++++++++++++++++++++

The spectral density in direction (2D) spectra doesn't have a directional convention in itself. Dnora has therefore implemented an additional directional convention attribute to keep track of this. They are stored in an ``Enum`` at ``dnora.type_manager.spectral_conventions.SpectralConvention``

Note, that all conventions use degrees, and they only describe how the wave direction is rotated and the to/from convention.

The available spectral conventions are:

OCEAN
---------

Oceanic convention: 

1. Direction TO
2. North = 0, East = 90
3. directional vector is monotonically increasing, e.g. [0,10,.....,340,350]
   
Value 180 means wave are propagation from North towards South.


MET
---------

Oceanic convention: 

1. Direction FROM
2. North = 0, East = 90
3. directional vector is monotonically increasing, e.g. [0,10,.....,340,350]
   
Value 180 means wave are propagation from South towards North.


WW3
---------

WAVEWATCH III convention: 

1. Direction TO
2. North = 0, East = 90
3. directional vector ordered according to the mathematical convention e.g. [90, 80, 70, ..., 0, 350, ..., 110, 100]
   
Value 180 means wave are propagation from North towards South.

In other words, all the values are in oceanographic convention, but the directional vector is reordered.

MATH
---------

Mathematical convention:

1. Direction TO
2. North = 90, East = 0
3. directional vector is monotonically increasing, e.g. [0,10,.....,340,350]

Value 180 means wave are propagation from East towards West.


MATHVEC
---------

Mathematical convention with reordered directional vector:

1. Direction TO
2. North = 90, East = 0
3. directional vector ordered according to the mathematical convention e.g. [90, 80, 70, ..., 0, 350, ..., 110, 100]

Value 180 means wave are propagation from East towards West.

UNDEFINED
------------

Included for purposes when se directional convention has not been set. This should only happen in experimental situatios.


