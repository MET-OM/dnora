.. dnora documentation master file, created by
   sphinx-quickstart on Sun Nov 15 14:18:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dnora's documentation!
=====================================
**DNORA** is a Python package for dynamical downscaling of NORA wave hindcast using the spectral wave models SWAN or WAVEWATCH III and wave-flow model SWASH.

This documentation is for the dnora version 2. 

The package contains functions that:
  * create a high-resolution grid using open-access bathymetry/topography datasets,
  * prepare the boundary spectra in the right directional convention from different sources
  * prepare forcing files (wind, current, ice, etc.) for the spectral models
  * create input parameter files based on the defined grid and files that have been prepared
  * run the wave models


Installing **DNORA**
=============================================
1. Install anaconda3 or miniconda3
2. Clone dnora:

.. code-block:: bash

   $ git clone https://github.com/MET-OM/dnora.git
   $ cd dnora/

3. Create environment with the required dependencies and install dnora

.. code-block:: bash

  $ conda config --add channels conda-forge
  $ conda env create -f environment.yml
  $ conda activate dnora
  $ pip install -e .
  
To update the enviroment using a new environment.yml, run:

.. code-block:: bash

   $ conda env update --file environment.yml --prune

Basic example
=============================================

DNORA was originally greated to easily downscale the NORA3 hindcast. While it can now do much more, here is how you do that. The code runs a 500 m SWAN model for one day using EMODNET topography, NORA3 wind forcing, and NORA3 boundary spectra::

   import dnora as dn
   
   grid = dn.grid.EMODNET(lon=(5.35, 5.6), lat=(59.00, 59.17), name="Boknafjorden")
   grid.import_topo()
   grid.set_spacing(dm=500)
   grid.mesh_grid()
   grid.set_boundary_points(dn.grid.mask.Edges(["N", "W", "S"]))
   
   model = dn.modelrun.NORA3(grid, year=2020, month=2, day=15)
   model.import_wind()
   model.import_spectra()
   
   exp = dn.export.SWAN(model)
   exp.export_grid()
   exp.export_wind()
   exp.export_spectra()
   
   exe = dn.executer.SWAN(model)
   exe.write_input_file()
   exe.run_model()

.. code-block:: python

Before running the example you have to have SWAN installed, since it is not a part of the DNORA package.

The basic workflow in DNORA scripts follow the same logic
  * Define an area you are working with by creating a Grid-object
  * Define a time period you are working with by creating a ModelRun-object
  * Import the data from the source you want
  * Define an exporter to export the data in the format you want
  * Define an executer to write input files and run the model

Specifically, the import from different sources (e.g. MET Norway, ECMWF) and the writing data in differen formats (e.g. SWAN, WAVEWATCH III) are separated, and you can always use any combination you want, while dnora takes care of making sure the data is in the right format for the model (e.g. spectral conventions).


The Grid-object
=====================================

+++++++
Setting area and resolution
+++++++

The grid object is initialized with the following command::

   import dnora as dn
   grid = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.set_spacing(dm=250) # Set spacing to around 250 metres
   
.. code-block:: rst

A desired grid spacing can be set either by providing a desired grid spacing in metres (approximate) or defining the amounts of grid points (exact)::

   grid.set_spacing(dm=250) # Set spacing to around 250 metres
   grid.set_spacing(nx=291, ny=249) # Create 291 (lon) x 249 (lat) grid points

.. code-block:: rst

Both  of  these  options  will  convert  the  input  to  the  native  resolution  in  longitude  and latitude.  These can, of course, also be set directly by::

   grid.set_spacing(dlon=0.0048780, dlat=0.0022573)

.. code-block:: rst

In this case ``dlon`` and ``dlat`` are not exact.  If an exact resolution needs to be forced, the ``floating_edge``-option can be used, e.g.,::

   grid.set_spacing(dlon=1/205, dlat=1/443, floating_edge=True)

.. code-block:: rst

This will enforce the resolution and instead change the initially set area slightly (if needed). The main properties of the grid can be accessed by methods::

   grid.lon() # Longitude vector
   grid.lat() # Latitude vector
   grid.name # Name given at initialization
   grid.nx() # Amount of point in longitude direction
   grid.ny() # Amount of point in latitude direction
   grid.size() # Tuple (nx, ny)

.. code-block:: rst

+++++++
Importing and meshing a topography
+++++++

The grid created above is empty. To get bathymetrical data we need to import it. Trivially we can just test it by "importing" a bathymetry with constant 50 meter depth::

   grid.import_topo(dn.read.generic.ConstantData(), topo=50)

.. code-block:: rst

The "imported" raw bathymetry can be found at ``grid.raw()`` if you want to look at it (the corresponding xarray Dataset is ``grid.raw().ds()``), but what you really want to do is mesh this raw topography to the grid spacing you defined::

   grid.mesh_grid()

.. code-block:: rst

This interpolates the "raw" topography, which can be either gridded or non-gridded, to the spacing you defined for your grid object. The meshing is a wrapper around ``scipy.interpolate.griddata`` using nearest neighbour interpolation. Change this by providing e.g. ``method='linear'``.

**NB!** DNORA has a system for keeping track of folders and possible data that has been downloaded automatically so it can be reused between projects. See section "Folders, filenames and URL's in DNORA" for more information.

If you want to modicy the meshed data, there are some functions to do that. E.g. to set everything less than 1 m deep to land, and everything less than 2 m deep to 2 m::

   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=1, to_land=True)
   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=2)

.. code-block:: rst

The possible data sources are:

+++++++
EMODNET (``dnora.read.grid.EMODNET``)
+++++++

Covers Europe with approximately 115 m resolution. To use this reader you can also just set the grid a subclass::

   grid = dn.grid.EMODNET(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo() # No reader needs to be specified

.. code-block:: rst

The EMODNET reader automatically identifies the tiles needed to cover the area, and downloads them if they do not exist in the folder.

Reference: https://emodnet.ec.europa.eu/en/bathymetry

+++++++
GEBCO (``dnora.read.grid.GEBCO``)
+++++++

Global coverage approximately 500 m resolution. To use this reader you can also just set the grid a subclass::

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(filename='name_of_downloaded_tile.nc') # No reader needs to be specified

.. code-block:: rst

The GEBCO reader is not as automated as the EMODNET reader, and basically wants to read custom tiles that have already been downloaded. You can download tiles here: https://download.gebco.net/ 

However, it does make use of an API provided by (https://github.com/cywhale/gebco). This means that if no filename is provided, then the relevant data is automatically fetched using the API. However, this doesn't work for very large areas. Recommended practice is therefore to download a faiirly large area once, and then make use of that several times.

The functionality will be improved in future versions.

Refereces: https://gebco.net, https://github.com/cywhale/gebco

+++++++
Karverket-50m (``dnora.read.grid.KartverketNo50m``)
+++++++

High-resolution coverage for the Norwegian coast at a 50 m resolution provided by the Norwegian Mapping Authority (Kartverket). To use this reader you can also just set the grid a subclass::

   grid = dn.grid.Kartverket(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(tile='B1008')
.. code-block:: rst

This reader is not automated. This means that you need to manually check which tile(s) you need and download them from Mapping Authorities web service, extract it and provide the tile to the reader.

References: https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e

+++++++
Netcdf (``dnora.read.generic.Netcdf``)
+++++++

Reads generic Netcdf files with structured data. E.g.::

   grid = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(dn.read.generic.Netcdf(), filename='some_random_topo.nc', folder='')

.. code-block:: rst

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.

+++++++
PointNetcdf (``dnora.read.generic.PointNetcdf``)
+++++++

Reads generic Netcdf files with unstructured data. This can be handy e.g. when you have donwloaded e.g. GEBCO data using the automatic API and wan't to save it for later::

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo()
   grid.raw().ds().to_netcdf('GEBCO_Bergen.nc')

   # Re-import from file
   grid.import_topo(dn.read.generic.PointNetcdf(), filename='GEBCO_Bergen.nc', folder='')

.. code-block:: rst

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.

+++++++
Boundary points
+++++++

Setting boundary points is now only important for being able to write the grid-files, but are also of consequence when importing boundary spectra. E.g. to set the norther, western and southern edges of the grid as boundaries::

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

.. code-block:: rst

Information about the boundary points that are set in the grid can be accessed using methods::

   grid.boundary_mask() # Boolean array where True is a boundary point
   grid.boundary_points() # Array containing a longitude, latitude list of the boundary points

.. code-block:: rst

The possible logics to define bounaries are:

+++++++
Whole edges (``dnora.grid.mask.Edges``)
+++++++

Set all wet points along a certain edge as boundary points::

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

.. code-block:: rst

This is probably what you want for spectral models.

+++++++
Mid points (``dnora.grid.mask.MidPoint``)
+++++++

Same as Edges, but only set the middle point::

   grid.set_boundary_points(dn.grid.mask.MidPoint(edges=['N', 'W'])

.. code-block:: rst

This might be needed in some phase resolving models.

+++++++
A list of points (``dnora.grid.mask.LonLat`` or ``dnora.grid.mask.XY``)
+++++++

Give a list of longitudes and latitudes (or x- and y-coordinates if grid is cartesian) and set the nearest grid point as a boundary point::

   grid.set_boundary_points(dn.grid.mask.LonLat(lon=(5.0, 5.1), lat=(60.3, 60.3)))

.. code-block:: rst

This is probably not used for setting bounadries, but this same function can be used to set a mask of e.g. spectral output points.

+++++++
All or none of the grid (``dnora.grid.mask.All`` or ``dnora.grid.mask.Clear``)
+++++++

Provided for completeness and used internally e.g. on one point grids::

   grid.set_boundary_points(dn.grid.mask.Clear())

.. code-block:: rst

The ModelRun-object
=====================================

++++++
Defining a time period
++++++

The ``ModelRun``-object is the second central object and contain all forcing and boundary data. This object is always defined for a certain grid and a certain time::

   model = dn.modelrun.ModelRun(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

.. code-block:: rst

If you need to prepare forcing files for a long hindcast, it might pay of to use e.g. monthly files. Then a good possibility can be to define the time with the alternative format::

   model = dn.modelrun.ModelRun(grid, year=2020, month=5)

.. code-block:: rst

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-05-01T23:00'``. If you need the first hour of the next month for e.g. hotstart-purposes, then::

   model = dn.modelrun.ModelRun(grid, year=2020, month=5, hotstart_hour=True)

.. code-block:: rst

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-06-01T00:00'``

+++++++++++++++++++++
Importing forcing data
+++++++++++++++++++++

The bathymetrical data has already been imported to the ``Grid``-object. All other data will be imported at this stage. To e.g. import NORA3 wind forcing that covers the area and time period::

   model.import_forcing(dn.read.wind.metno.NORA3())

.. code-block:: rst

In practice, you can use often use subclasses (defined in ``dnora.modelrun.templates``) of the ``ModelRun``-object and then you don't have to specify the reader. For example, if you just want to work with NORA3 data, then the above code simplifies to::

   model = dn.modelrun.NORA3(grid, year=2020, month=5)
   model.import_forcing()

.. code-block:: rst

You can also now import boundary spectra and ice simply by calling those import methods::

   model.import_spectra()
   model.import_ice()

.. code-block:: rst

+++++++++++++++++++++
Import of gridded data
+++++++++++++++++++++

You typically want to import two types of data: gridded or non-gridded. Examples of gridded data are wind, ice and currents. An example of non-gridded data is boundary spectra.

For gridded data the import-method automaticall gets all data that covers the grid you provided, but expands it by 20% to make sure the edges are covered (needed for interpolation etc.). If the imported area is too small or too large, this can be controlled by a keyword::

   model.import_forcing(expansion_factor=1.3)

.. code-block:: rst

will expand the are by 30%. The default value here is 1.2.

+++++++++++++++++++++
Import of non-gridded data
+++++++++++++++++++++

For non-gridded data we have some more options, since we might want different ways to choose e.g. boundary spectra. This is technically done by separating the logic of choosing the points to import and the actual import of data. Usually DNORA takes care of this automatically, and the logic is the following::

1) If no boundary points are set, then we will import all spectra that fall within the area of the grid that has been expanded by 50% (``expansion_factor = 1.5``)
2) If boundary points are set, then we will go through all boundary points in order, and find the nearest spectra in the database. No duplicates are given, so if all boundary points have the same nearest spectra, then only one spectum is imported.
3) If the grid consists of a single point, then we find the nearest spectrum to that point, regardless of if it has been set to a boundary point.

**Case 1)** the points are technically (and automatically) chosen by the ``PointPicker`` class ``dnora.pick.Area``. The user can change the area if e.g. the grid is small and the resolution of the spectra are low. To just take all NORA3 spectra the covers an area expanded to double size::

   area = dn.grid.Grid(lat=(60.1, 60.2), lon=(4.2, 4.3), name='TestArea')
   model = dn.modelrun.NORA3(area, year=2020, month=5)
   model.import_spectra(expansion_factor=2)

.. code-block:: rst

**Case 2)** we set the western edge to boundary points and then grab the nearest spectrum (technically implemented in ``dnora.pick.NearestGridPoint``) to each point::

   grid = dn.grid.Grid(lat=(60.0, 60.4), lon=(4.0, 4.3), name='TestGrid')
   grid.set_spacing(dlon=0.1, dlat=0.1)
   grid.set_boundary_points(dn.grid.mask.Edges(['W'])

   model = dn.modelrun.NORA3(grid, year=2020, month=5)
   model.import_spectra()

.. code-block:: rst

To speed up the process, DNORA actually limits the search for the nearest spectrum to an area expanded by 3 degrees in latitude and 6 degrees in longitude. This is mainly done to avoid having to calculate distances to points very far away, since fast cartesian calculations are not available if the database contains points above 84 degrees latitude. This is not visible to the user. The user can, however, choose to discard spectra that are too far from the boundary point. For example, if we rather not import any spectrum at all if the closest one is over 7 km away, then::

   model.import_spectra(max_dist=7)

.. code-block:: rst

By default, no max distance is applied.

*If* the user wants to import spectra for the entire area even though boundary points have been defined, then the default behaviour can be overridden by explicitly providing a ``PointPicker``. If, in addition, we want to make sure that no spectra from outside the grid area is imported, we can set the ``expansion_factor``::

   model.import_spectra(point_picker=dn.pick.Area(), expansion_factor=1)

.. code-block:: rst

**Case 3)** is implemented so that the user can easily just get the spectrum that is nearest to a single point from a database without having to set a technical boundary point. Simply::

   point = dn.grid.Grid(lat=60.0, lon=4.0, name='TestPoint')
   
   model = dn.modelrun.NORA3(point, year=2020, month=5)
   model.import_spectra()

.. code-block:: rst

Also here the ``max_dist`` keyword can be used. Note, that since DNORA uses a standard format for 2D-spectra, it can automatically integrate the spectra to omni-directional spectra and to some basic integrated parameters. E.g. To get the omnidirectional spectra of the NORA3 2D-spectra from above::

   model.spectra_to_1d()

.. code-block:: rst

You can now find the 1D-spectra at ``model.spectra1d()`` and the corresponding xarray Dataset at ``model.spectra1d().ds()``. To calculate integrated parameters from the 1D-spectra, use::
   
   model.spectra_to_waveseries()

.. code-block:: rst

You can now find the wave parameters ``model.waveseries()``. The method automatically calculates the 1D-spectra as an intermediate step if you have only imported 2D-spectra. In other words, the following code works::

   model.import_spectra()
   model.spectra_to_waveseries()

.. code-block:: rst

And after that you have 2D-spectra, 1D-spectra and wave parameters in ``model.spectra()``, ``model.spectra1d()`` and ``model.waveseries()`` respectively.

Exporting data
=====================================

The imported data is now in a model agnostic format inside the ``ModelRun`` object. We now need to export them to files so that the wave model can use them as forcing files. This is done by choosing a model-specific exporter. For example, if we want to run SWAN::

   exp = dn.export.SWAN(model)
   exp.export_grid()
   exp.export_wind()
   exp.export_spectra()

.. code-block:: rst

This exports the data to ASCII files that follows the SWAN format. It also makes sure that the directional convention of the spectra matches that that SWAN is expecting (in this case coming from). The filenames are DNORA defaults. If you need to change this, it can be done. See the section on Filenames and folders. 

To prepare forcing files for WAVEWATCH III (including converting to the proper spectral convention), just use a different exporter::

   exp = dn.export.WW3(model)
   
.. code-block:: rst

The available exporters are:

+++++++++++++++++++++
SWAN (``dn.export.SWAN``)
+++++++++++++++++++++

Can write grid, wind, spectra, current, waterlevel and ice. Writes ASCII files and all data into one file.

+++++++++++++++++++++
WW3 (``dn.export.WW3``)
+++++++++++++++++++++

Can write grid, wind, spectra and triangular grids. Writes netcdf's (wind files into monthly files). Since there are different verisions of the WAVEWATCH III code floating around, the spectral writer has some options:

``one_file=False`` (default: ``True``): Writes every point into a separate file.

``squeeze_lonlat=True`` (default: ``False``): Write the longitudes and latitudes without a time dimension.

``convention='ww3'`` (default: ``'ww3'``): Force a different convention for the spectra (see section on DNORA Spectral conventions). Use only if you know what you are doing.

``var_names={'efth': 'SPEC'}``: Change the names of the variables in the netcdf-files. E.g. change the spectra to ``'SPEC'`` (default: ``'efth'``).

+++++++++++++++++++++
SWASH (``dn.export.SWASH``)
+++++++++++++++++++++

Can write grid, wind, and spectra. Uses same ASCII format as SWAN.

+++++++++++++++++++++
REEF3D (``dn.export.REEF3D``)
+++++++++++++++++++++

Can write spectra1d.

+++++++++++++++++++++
Netcdf (``dn.export.Netcdf``)
+++++++++++++++++++++

Can write anything and is model-agnistic. Default writes all data into one file, but you can specify ``daily_files=True`` or ``monthly_files=True`` in the export method.

+++++++++++++++++++++
Cacher (``dn.export.Cacher``)
+++++++++++++++++++++

Used by the automatic caching function in DNORA (see section on Caching), and is probably not interesting for a user.


Running the model
=====================================

The process of running the model includes writing the input files and possibly running the executables. The needs and process is a bit different for differnt models. All models need to be compiled, installed and set up separately. Please see the section on Dependencies for some additional information.

If the model can output wave spectra (namely SWAN and WAVEWATCH III), then output can be requested by setting output-points in the grid. E.g. to request the model to output spectra at a single location::

   grid.set_output_points(dn.grid.mask.LonLat(lon=5.4, lat=59.1))

.. code-block:: rst

+++++++++++
SWAN
+++++++++++

If the files have been exported, SWAN can be ran with::

   exe = dn.executer.SWAN(model)
   exe.write_input_file()
   exe.run_model()

.. code-block:: rst

The first function writes the input file (``.swn``) for the model. In writing the input file it makes use that it already knows the resolution of the grid and forcing files, the names to where the data has been exported, the wanted spectral output points etc. There are options to control the writing of the input file:

``calibrate: dict`` default values: ``{"wind": 1, "wcap": 0.5000e-04, "waterlevel": 1, "current": 1, "ice": 1}``

   ``'wcap'``: the white capping coefficient
   
   ``'wind'``: factor that the wind is multiplied with, e.g. 1.1 increases the wind speed by 10%.
   
   ``'current'``: factor that the currents are multiplied with, e.g. 1.1 increases the current speed by 10%.
   
   ``'waterlevel'``: factor added to the waterlevels, e.g. 2 means 2 meters are added to the waterlevel.
   
   ``'ice'``: factor for ice (spefic meaning?)

``use_wind, use_waterlevel, use_spectra, use_current, use_ice: bool`` (default: ``True``): Set to false if you have maybe exported data that you don't want to use. By default forcing data is used if present (i.e. has been exported).

``structures: list[dict]``: Can be used to set breakwaters of a given shape, and with given transparency and reflection. From the doc string of the method::

     """structures given in format:
     E.g. One triangle and one line
     [{'lon': (5,5.1,4.9), 'lat': (60,59.9,59,9), 'trans': 0.3, 'refl': 0.0, 'closed': True, 'name': 'closed triangle'},
     {'lon': (4,4.1), 'lat': (60.1,60.1,60.1),'name': 'breakwater'}
     ]
   
     'trans' is the transparency (0 completely blocked, 1 completely open)
     'refl' is the reflection (default 0.0)
   
     if 'trans' or 'refl' is not given for an object, then the previous object's values are used.
     Especially, if values are only given for the first object, then they are used for all objects.
   
     'closed': True [default False] means that the first point is repeated in the input file, thus effectively closing the structure
   
     'name' is optional, but will be printed to the swn file if provided to give a better overview
     """

.. code-block:: rst

The physics are currently set to a given standard that has been found useful in the Norwegian fjords (see e.g. Christakos et al (202x))::

   GEN3 WESTH cds2=5e-05 AGROW
   FRICTION JON 0.067 
   PROP BSBT 
   NUM ACCUR NONST 1 

.. code-block:: rst

If the user needs to change something else than the whitecapping coefficient, the input file needs to be edited manually before running the model.

The second function simply executes the model. It runs with OpenMP parallelization. To change the amount of processes, ute the keyword ``nproc=2`` (default: 4).



File names and directories
=====================================

The default file names and directories used in dnora are defined in the ``defaults.py``-file. Different default can be generated for different models, but the styles are not inherently linked to a certain models (see example below).

The default file names and folders used by the different writers are set within the writers, and they convey their preference to the ``ModelRun``-object. These defaults are used if the user doesn't provide anything explicitly. For example, the default file name for writing wind forcing for the SWAN model is::

   wind#Forcing#Grid#T0_#T1.asc


.. code-block:: rst

where ``#Forcing`` and ``#Grid`` will be replaced by the name of the forcing and grid of the ``ModelRun``-object, and ``#T0`` and ``#T1`` will be replaced by the start and end times, formatted according to the default format: ``%Y%m%d``.

Let's say we want to run an operation version of the SWAN model for the grid name "Sula", and want to have the forcing file name in the format: ``WIND_Sula_2018010106Z.asc``, where the time is the start time. The first way to do this is to provide this information to the method doing the writing, e.g.::

   model.export_forcing(wnd.write.SWAN(), filestring="WIND_#Grid_#T0Z.asc", datestring="%Y%m%d%H")

.. code-block:: rst

If this is used often, then these values can be added to ``defaults.py`` under the name "SWAN_oper". Then we can simply set the preference format upon initialization of the ``ForcingWriter``::

   model.export_forcing(wnd.write.SWAN(out_format='SWAN_oper'))

.. code-block:: rst

The third level is to actualy create a new template for this type of model runs, which can be done (for example) as a subclass of the ``SWAN``-template::

   class SWAN_oper(SWAN):
       def _get_forcing_writer(self):
           return wnd.write.SWAN(out_format='SWAN_oper')

.. code-block:: rst

Download data
=============================================
You can easily execute only a part of the workflow. For example, say you only want to download NORA3 wind and wave data around the area of Bergen, Norway for January 2020, but don't really want to worry about any specific model or model grid::

   import dnora as dn
   area = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   
   model = dn.modelrun.NORA3(area, year=2020, month=1)
   model.import_wind()
   model.import_spectra()

.. code-block:: rst

Dnora automatically expands the download area with 20% for wind and 50% for wave spectra to make sure you get data covering the whole area. This might be important when the wave model interpolates the data to its own grid. If you want exactly the area you defined, use the keyword ``expansion_factor=1`` in both import methods.

Say you wan't your data in netcdf files, but you want to make sure that your wave spectra have the convention 'coming from'. You can then first set the convention, and use the Netcdf-exporter to write the data::

   model.spectra().set_convention('met')

   exp = dn.exporter.Netcdf(model)
   exp.export_wind()
   exp.export_spectra()

.. code-block:: rst

In this case the NORA3 wave spectra actually had the directional convention 'going to', but you don't neew to know that. If the convention would have been right from the start, the ``set_convention``-method would have done nothing. If you would have wanted ERA5 data instead, just change the line to::

   model = dn.modelrun.ERA5(area, year=2020, month=1)

.. code-block:: rst





Dependencies
=====================================
1. Installation of **SWAN**. The latest SWAN version can be downloaded from https://sourceforge.net/projects/swanmodel/files/swan/ . The installation procedure can be found in: https://swanmodel.sourceforge.io/online_doc/swanimp/node4.html

2. Installation of **WAVEWATCH III**. The latest model version can be downloaded from https://github.com/NOAA-EMC/WW3 . The installation procedure can be found in: https://github.com/NOAA-EMC/WW3/wiki/Quick-Start

3. Installation of **SWASH**. The latest model version can be downloaded from https://sourceforge.net/projects/swash/ . The installation procedure can be found in: https://swash.sourceforge.io/online_doc/swashimp/node4.html

4. Installation of **HOS-ocean**. The latest HOS-ocean version can be downloaded from https://github.com/LHEEA/HOS-ocean . The installation procedure can be found in: https://github.com/LHEEA/HOS-ocean/wiki/Installation

5. Installation of **REEF3D**. The latest REEF3D and DIVEMesh versions can be downloaded from https://github.com/REEF3D . The installation procedure can be found in: https://reef3d.wordpress.com/installation/ . Note: Only REEF3D::FNPF has been tested in dnora.

To run the models within dnora, the paths, where the models are installed, need to be defined in .bashrc, e.g., ::

   export PATH=${PATH}:/home/user/Programs/swan

   export PATH=${PATH}:/home/user/Programs/swash

   export PATH=${PATH}:/home/user/Programs/HOS-ocean/bin

   export PATH=${PATH}:/home/user/Programs/REEF3D_xx/DIVEMesh/bin
   export PATH=${PATH}:/home/user/Programs/REEF3D_xx/REEF3D/bin

   source ~/.bashrc


.. code-block:: rst

3. Fimex is a the File Interpolation, Manipulation and EXtraction library for gridded geospatial data (more info in \url{httpshttps://wiki.met.no/fimex/start}). Fimex is used in DNORA for the preparation of forcing fields (wind). **In case of running the spectral model without wind forcing, the fimex installation can be neglected**.  A detailed installation procedure can be find in https://wiki.met.no/fimex/install. For a quick installation in Linux-Ubuntu, you can follow the steps: open the Synaptic Package Manager (to install it: sudo apt-get install synaptic
) --> Settings --> Repositories --> Other Software --> Add:  **ppa:met-norway/fimex** to your system's Software Sources (see https://launchpad.net/~met-norway/+archive/ubuntu/fimex). Then, search and mark for installation a latest version of fimex, and press apply installation. Check the installation (usually it is installed in /usr/bin/) by typing in command line: fimex or fimex-xxx where xxx is the version number. In case that only fimex-xxx works then add a symbolic link::

      cd /usr/bin
      sudo ln -s fimex-xxx fimex

.. code-block:: rst

Folders, filenames and URL's in DNORA
=====================================

[To be added]

.. toctree::
   :maxdepth: 2
   :caption: Contents:
