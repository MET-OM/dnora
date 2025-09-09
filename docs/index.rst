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

Quick Installation: 
The easiest way to install dnora is to use pip:

.. code-block:: bash

   $ pip install dnora

Alternatively, you can clone the GitHub repository:

1. Install anaconda3 or miniconda3
2. Clone dnora:

.. code-block:: bash

   $ git clone https://github.com/MET-OM/dnora.git
   $ cd dnora/

3. Create an environment with the required dependencies and install dnora

.. code-block:: bash

  $ conda config --add channels conda-forge
  $ conda env create -f environment.yml
  $ conda activate dnora2
  $ pip install .
  
To update the environment using a new environment.yml, run:

.. code-block:: bash

   $ conda env update --file environment.yml --prune

For JupyterHub platforms such as WEKEO, please use:

.. code-block:: bash

   $ conda create -n dnora python=3.10 fimex=1.8.1
   $ conda activate dnora
   $ python -m pip install dnora

Basic example
=============================================

DNORA was originally greated to easily downscale the NORA3 hindcast. While it can now do much more, here is how you do that. The code runs a 500 m SWAN model for one day using EMODNET topography, NORA3 wind forcing, and NORA3 boundary spectra:

.. code-block:: python

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

The grid object is initialized with the following command:

.. code-block:: python

   import dnora as dn
   grid = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.set_spacing(dm=250) # Set spacing to around 250 metres
   


A desired grid spacing can be set either by providing a desired grid spacing in metres (approximate) or defining the amounts of grid points (exact):

.. code-block:: python

   grid.set_spacing(dm=250) # Set spacing to around 250 metres
   grid.set_spacing(nx=291, ny=249) # Create 291 (lon) x 249 (lat) grid points

Both  of  these  options  will  convert  the  input  to  the  native  resolution  in  longitude  and latitude.  These can, of course, also be set directly by:

.. code-block:: python

   grid.set_spacing(dlon=0.0048780, dlat=0.0022573)


In this case ``dlon`` and ``dlat`` are not exact.  If an exact resolution needs to be forced, the ``floating_edge``-option can be used, e.g.,:

.. code-block:: python

   grid.set_spacing(dlon=1/205, dlat=1/443, floating_edge=True)


This will enforce the resolution and instead change the initially set area slightly (if needed). The main properties of the grid can be accessed by methods:

.. code-block:: python

   grid.lon() # Longitude vector
   grid.lat() # Latitude vector
   grid.name # Name given at initialization
   grid.nx() # Amount of point in longitude direction
   grid.ny() # Amount of point in latitude direction
   grid.size() # Tuple (nx, ny)


+++++++
Importing and meshing a topography
+++++++

The grid created above is empty. To get bathymetrical data we need to import it. Trivially we can just test it by "importing" a bathymetry with constant 50 meter depth:

.. code-block:: python

   grid.import_topo(dn.read.generic.ConstantData(), topo=50)

The "imported" raw bathymetry can be found at ``grid.raw()`` if you want to look at it (the corresponding xarray Dataset is ``grid.raw().ds()``), but what you really want to do is mesh this raw topography to the grid spacing you defined:

.. code-block:: python

   grid.mesh_grid()

This interpolates the "raw" topography, which can be either gridded or non-gridded, to the spacing you defined for your grid object. The meshing is a wrapper around ``scipy.interpolate.griddata`` using nearest neighbour interpolation. Change this by providing e.g. ``method='linear'``.

**NB!** DNORA has a system for keeping track of folders and possible data that has been downloaded automatically so it can be reused between projects. See section "Folders, filenames and URL's in DNORA" for more information.

If you want to modicy the meshed data, there are some functions to do that. E.g. to set everything less than 1 m deep to land, and everything less than 2 m deep to 2 m:

.. code-block:: python

   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=1, to_land=True)
   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=2)


The possible data sources are:

+++++++
EMODNET (``dnora.read.grid.EMODNET``)
+++++++

Covers Europe with approximately 115 m resolution. To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.EMODNET(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo() # No reader needs to be specified

The EMODNET reader automatically identifies the tiles needed to cover the area, and downloads them if they do not exist in the folder.

Reference: https://emodnet.ec.europa.eu/en/bathymetry

+++++++
GEBCO (``dnora.read.grid.GEBCO``)
+++++++

Global coverage approximately 500 m resolution. To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(filename='name_of_downloaded_tile.nc') # No reader needs to be specified

The GEBCO reader is not as automated as the EMODNET reader, and basically wants to read custom tiles that have already been downloaded. You can download tiles here: https://download.gebco.net/ 

However, it does make use of an API provided by (https://github.com/cywhale/gebco). This means that if no filename is provided, then the relevant data is automatically fetched using the API. However, this doesn't work for very large areas. Recommended practice is therefore to download a faiirly large area once, and then make use of that several times.

The functionality will be improved in future versions.

References: https://gebco.net, https://github.com/cywhale/gebco

+++++++
Karverket-50m (``dnora.read.grid.KartverketNo50m``)
+++++++

High-resolution coverage for the Norwegian coast at a 50 m resolution provided by the Norwegian Mapping Authority (Kartverket). To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.Kartverket(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(tile='B1008')

This reader is not automated. This means that you need to manually check which tile(s) you need and download them from Mapping Authorities web service, extract it and provide the tile to the reader.

References: https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e

+++++++
Netcdf (``dnora.read.generic.Netcdf``)
+++++++

Reads generic Netcdf files with structured data. E.g.:

.. code-block:: python

   grid = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(dn.read.generic.Netcdf(), filename='some_random_topo.nc', folder='')

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.

+++++++
PointNetcdf (``dnora.read.generic.PointNetcdf``)
+++++++

Reads generic Netcdf files with unstructured data. This can be handy e.g. when you have donwloaded e.g. GEBCO data using the automatic API and wan't to save it for later:

.. code-block:: python

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo()
   grid.raw().ds().to_netcdf('GEBCO_Bergen.nc')

   # Re-import from file
   grid.import_topo(dn.read.generic.PointNetcdf(), filename='GEBCO_Bergen.nc', folder='')

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.

+++++++
Boundary points
+++++++

Setting boundary points is now only important for being able to write the grid-files, but are also of consequence when importing boundary spectra. E.g. to set the norther, western and southern edges of the grid as boundaries:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

Information about the boundary points that are set in the grid can be accessed using methods:

.. code-block:: python

   grid.boundary_mask() # Boolean array where True is a boundary point
   grid.boundary_points() # Array containing a longitude, latitude list of the boundary points

The possible logics to define bounaries are:

+++++++
Whole edges (``dnora.grid.mask.Edges``)
+++++++

Set all wet points along a certain edge as boundary points:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

This is probably what you want for spectral models.

+++++++
Mid points (``dnora.grid.mask.MidPoint``)
+++++++

Same as Edges, but only set the middle point:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.MidPoint(edges=['N', 'W'])

This might be needed in some phase resolving models.

+++++++
A list of points (``dnora.grid.mask.LonLat`` or ``dnora.grid.mask.XY``)
+++++++

Give a list of longitudes and latitudes (or x- and y-coordinates if grid is cartesian) and set the nearest grid point as a boundary point:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.LonLat(lon=(5.0, 5.1), lat=(60.3, 60.3)))

This is probably not used for setting bounadries, but this same function can be used to set a mask of e.g. spectral output points.

+++++++
All or none of the grid (``dnora.grid.mask.All`` or ``dnora.grid.mask.Clear``)
+++++++

Provided for completeness and used internally e.g. on one point grids:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Clear())

The ModelRun-object
=====================================

++++++
Defining a time period
++++++

The ``ModelRun``-object is the second central object and contain all forcing and boundary data. This object is always defined for a certain grid and a certain time:

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

If you need to prepare forcing files for a long hindcast, it might pay of to use e.g. monthly files. Then a good possibility can be to define the time with the alternative format:

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, year=2020, month=5)

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-05-01T23:00'``. If you need the first hour of the next month for e.g. hotstart-purposes, then:

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, year=2020, month=5, hotstart_hour=True)

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-06-01T00:00'``

+++++++++++++++++++++
Importing forcing data
+++++++++++++++++++++

The bathymetrical data has already been imported to the ``Grid``-object. All other data will be imported at this stage. To e.g. import NORA3 wind forcing that covers the area and time period:

.. code-block:: python

   model.import_forcing(dn.read.wind.metno.NORA3())

In practice, you can use often use subclasses (defined in ``dnora.modelrun.templates``) of the ``ModelRun``-object and then you don't have to specify the reader. For example, if you just want to work with NORA3 data, then the above code simplifies to:

.. code-block:: python

   model = dn.modelrun.NORA3(grid, year=2020, month=5)
   model.import_forcing()

You can also now import boundary spectra and ice simply by calling those import methods:

.. code-block:: python

   model.import_spectra()
   model.import_ice()

You typically want to import two types of data: gridded or non-gridded. Examples of gridded data are wind, ice and currents. An example of non-gridded data is boundary spectra.

+++++++++++++++++++++
Import of gridded data
+++++++++++++++++++++

For gridded data the import-method automaticall gets all data that covers the grid you provided, but expands it by 20% to make sure the edges are covered (needed for interpolation etc.). If the imported area is too small or too large, this can be controlled by a keyword:

.. code-block:: python

   model.import_forcing(expansion_factor=1.3)

will expand the are by 30%. The default value here is 1.2.

+++++++++++++++++++++
Import of non-gridded data
+++++++++++++++++++++

For non-gridded data we have some more options, since we might want different ways to choose e.g. boundary spectra. This is technically done by separating the logic of choosing the points to import and the actual import of data. Usually DNORA takes care of this automatically, and the logic is the following::

1) If no boundary points are set, then we will import all spectra that fall within the area of the grid that has been expanded by 50% (``expansion_factor = 1.5``)
2) If boundary points are set, then we will go through all boundary points in order, and find the nearest spectra in the database. No duplicates are given, so if all boundary points have the same nearest spectra, then only one spectum is imported.
3) If the grid consists of a single point, then we find the nearest spectrum to that point, regardless of if it has been set to a boundary point.

**Case 1)** the points are technically (and automatically) chosen by the ``PointPicker`` class ``dnora.pick.Area``. The user can change the area if e.g. the grid is small and the resolution of the spectra are low. To just take all NORA3 spectra the covers an area expanded to double size:

.. code-block:: python

   area = dn.grid.Grid(lat=(60.1, 60.2), lon=(4.2, 4.3), name='TestArea')
   model = dn.modelrun.NORA3(area, year=2020, month=5)
   model.import_spectra(expansion_factor=2)

**Case 2)** we set the western edge to boundary points and then grab the nearest spectrum (technically implemented in ``dnora.pick.NearestGridPoint``) to each point:

.. code-block:: python

   grid = dn.grid.Grid(lat=(60.0, 60.4), lon=(4.0, 4.3), name='TestGrid')
   grid.set_spacing(dlon=0.1, dlat=0.1)
   grid.set_boundary_points(dn.grid.mask.Edges(['W'])

   model = dn.modelrun.NORA3(grid, year=2020, month=5)
   model.import_spectra()

To speed up the process, DNORA actually limits the search for the nearest spectrum to an area expanded by 3 degrees in latitude and 6 degrees in longitude. This is mainly done to avoid having to calculate distances to points very far away, since fast cartesian calculations are not available if the database contains points above 84 degrees latitude. This is not visible to the user. The user can, however, choose to discard spectra that are too far from the boundary point. For example, if we rather not import any spectrum at all if the closest one is over 7 km away, then:

.. code-block:: python

   model.import_spectra(max_dist=7)

By default, no max distance is applied.

*If* the user wants to import spectra for the entire area even though boundary points have been defined, then the default behaviour can be overridden by explicitly providing a ``PointPicker``. If, in addition, we want to make sure that no spectra from outside the grid area is imported, we can set the ``expansion_factor``:

.. code-block:: python

   model.import_spectra(point_picker=dn.pick.Area(), expansion_factor=1)

**Case 3)** is implemented so that the user can easily just get the spectrum that is nearest to a single point from a database without having to set a technical boundary point. Simply:

.. code-block:: python

   point = dn.grid.Grid(lat=60.0, lon=4.0, name='TestPoint')
   
   model = dn.modelrun.NORA3(point, year=2020, month=5)
   model.import_spectra()

Also here the ``max_dist`` keyword can be used. Note, that since DNORA uses a standard format for 2D-spectra, it can automatically integrate the spectra to omni-directional spectra and to some basic integrated parameters. E.g. To get the omnidirectional spectra of the NORA3 2D-spectra from above:

.. code-block:: python

   model.spectra_to_1d()


You can now find the 1D-spectra at ``model.spectra1d()`` and the corresponding xarray Dataset at ``model.spectra1d().ds()``. To calculate integrated parameters from the 1D-spectra, use:

.. code-block:: python
   
   model.spectra_to_waveseries()

You can now find the wave parameters ``model.waveseries()``. The method automatically calculates the 1D-spectra as an intermediate step if you have only imported 2D-spectra. In other words, the following code works:

.. code-block:: python

   model.import_spectra()
   model.spectra_to_waveseries()

And after that you have 2D-spectra, 1D-spectra and wave parameters in ``model.spectra()``, ``model.spectra1d()`` and ``model.waveseries()`` respectively.

Exporting data
=====================================

The imported data is now in a model agnostic format inside the ``ModelRun`` object. We now need to export them to files so that the wave model can use them as forcing files. This is done by choosing a model-specific exporter. For example, if we want to run SWAN:

.. code-block:: python

   exp = dn.export.SWAN(model)
   exp.export_grid()
   exp.export_wind()
   exp.export_spectra()

This exports the data to ASCII files that follows the SWAN format. It also makes sure that the directional convention of the spectra matches that that SWAN is expecting (in this case coming from). The filenames are DNORA defaults. If you need to change this, it can be done. See the section on Filenames and folders. 

To prepare forcing files for WAVEWATCH III (including converting to the proper spectral convention), just use a different exporter:

.. code-block:: python

   exp = dn.export.WW3(model)
   
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

If the model can output wave spectra (namely SWAN and WAVEWATCH III), then output can be requested by setting output-points in the grid. E.g. to request the model to output spectra at a single location:

.. code-block:: python

   grid.set_output_points(dn.grid.mask.LonLat(lon=5.4, lat=59.1))

+++++++++++
SWAN
+++++++++++

If the files have been exported, SWAN can be ran with:

.. code-block:: python

   exe = dn.executer.SWAN(model)
   exe.write_input_file()
   exe.run_model()

The first function writes the input file (``.swn``) for the model. In writing the input file it makes use that it already knows the resolution of the grid and forcing files, the names to where the data has been exported, the wanted spectral output points etc. There are options to control the writing of the input file:

``calibrate: dict`` default values: ``{"wind": 1, "wcap": 0.5000e-04, "waterlevel": 1, "current": 1, "ice": 1}``

   ``'wcap'``: the white capping coefficient
   
   ``'wind'``: factor that the wind is multiplied with, e.g. 1.1 increases the wind speed by 10%.
   
   ``'current'``: factor that the currents are multiplied with, e.g. 1.1 increases the current speed by 10%.
   
   ``'waterlevel'``: factor added to the waterlevels, e.g. 2 means 2 meters are added to the waterlevel.
   
   ``'ice'``: factor for ice (spefic meaning?)

``use_wind, use_waterlevel, use_spectra, use_current, use_ice: bool`` (default: ``True``): Set to false if you have maybe exported data that you don't want to use. By default forcing data is used if present (i.e. has been exported).

``structures: list[dict]``: Can be used to set breakwaters of a given shape, and with given transparency and reflection. From the doc string of the method:

.. code-block:: python

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

The physics are currently set to a given standard that has been found useful in the Norwegian fjords (see e.g. Christakos et al (202x)):

.. code-block:: python

   GEN3 WESTH cds2=5e-05 AGROW
   FRICTION JON 0.067 
   PROP BSBT 
   NUM ACCUR NONST 1 

If the user needs to change something else than the whitecapping coefficient, the input file needs to be edited manually before running the model.

The second function simply executes the model. It runs with OpenMP parallelization. To change the amount of processes, ute the keyword ``nproc=2`` (default: 4).

+++++++
WAVEWATCH III
+++++++

WAVEWATCH III is a bit more complicated to run compared to SWAN, since you need to pre-process the grid, forcing etc. The workflow here still needs to be improved, especially if using a super computer, but running a test locally for now works like this (after exporting all the data etc.):

.. code-block:: python

   exe = dn.executer.WW3(model)
   exe.write_grid_file()
   exe.write_wind_file()
   exe.write_spectra_file()
   exe.write_input_file()

   exe.run_grid()
   exe.run_wind()
   exe.run_spectra()
   exe.run_model()

This assumes that ``DNORA_WW3_PATH`` has been set (see section on environmental variables). Alternatively the path can be provided using the keyword ``model_folder``. Note, that the ``run_model`` method also automatically runs the post-processing programs to get gridded and spectral files from binary to netcdf files (``ounf`` and ``ounp``).

If you are just preparing the forcing an plan on copying the forcing files to a server and running WAVEWATCH III there, then you can speficy where the forcing files will be on the server when writing the input files for the pre-processors:

.. code-block:: python

   exe.write_grid_file(folder_on_server='/server/data/grids')
   exe.write_wind_file(folder_on_server='/server/data/winds')
   exe.write_spectra_file(folder_on_server='/server/data/spectra')

For the input file written for the spectral pre-processor (``write_spectra_file``) you can also specify the method WAVEWATCH III will use to interpolate the boundary spectra: ``method='nearest'`` (default) or ``'linear'``.

WAVEWATCH III can also be run with homogeneous input. You could creaty dummy constant files, but you can also set this when writing the input file:

.. code-block:: python

   exe.write_input_file(homog={"wind": (10, 0)})

Here, ``(10, 0)`` referres to the ``u`` and ``v`` components respectively. You can also specify a constant current and waterlevel: 

``homog={"wind": (10, 0), "current": (-1,0), "waterlevel": 3}``

The physics used by WAVEWATCH III (i.e. ST4, ST6) is determined at compile time, and specifying these are therefore not a part of DNORA. If you want to change other setting (such as BETAMAX), this is done in the ``namelists.nml`` file. It is empty by default, so this needs to be done manually for now.


Main data sets
=====================================

DNORA can read data from many different sources, but the main ones are given below, with a breif description and guide to usage.


++++++++
NORA3 (``dnora.modelrun.NORA3``)
++++++++

**Description**: NORA3 is an atmospheric and wave hindcasts downscaled from ERA5 to a 3 km resolution. The athmospheric hindcast covers North Sea, the Norwegian Sea, and the Barents Sea. The wave hindcast covers a larger area, and uses ERA5 winds outside the atmospheric domain. 

**Availablity in DNORA**: You can download the wind (3 km resolution), wave spectra (stored at about a 30 km resolution) and the ice used to force the wave model (3 km resolution). 

**Steps needed to access data**: None. Available on the MET Norway thredds server. 

**Minimal example**: A minimal example to get NORA3 data for an area and save it as netcdf is:

.. code-block:: python

   import dnora as dn
   
   area = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   model = dn.modelrun.NORA3(area, year=2020, month=2, day=15)
   model.import_wind()
   model.import_spectra()
   model.import_ice()
   
   exp = dn.export.Netcdf(model)
   exp.export_wind()
   exp.export_spectra()
   exp.export_ice()
   

**References**: 

Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., & Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of the North Sea, the Norwegian Sea, and the Barents Sea, Journal of Applied Meteorology and Climatology, 60(10), 1443-1464, https://doi.org/10.1175/JAMC-D-21-0029.1

Breivik, Ø., Carrasco, A., Haakenstad, H., Aarnes, O. J., Behrens, A., Bidlot, J.-R., et al. (2022). The impact of a reduced high-wind Charnock parameter on wave growth with application to the North Sea, the Norwegian Sea, and the Arctic Ocean. Journal of Geophysical Research: Oceans, 127, e2021JC018196. https://doi.org/10.1029/2021JC018196


++++++++
Operational WW3-4km and MEPS (``dnora.modelrun.WW3_4km``)
++++++++

**Description**: This is the operational wave and wind data of MET Norway that is run four times per day (00, 06, 12 and 18 UTC) for 67 (MEPS) and 73 (WW3) hours. The wave model is a 4 km WAVEWATCH III implementation forced by one deterministic member of the MEPS atmospheric forecasta, which is a 2.5 km Arome model.

**Availablity in DNORA**: You can download the wind (2.5 km resolution) and wave spectra (irregular spacing). Since this is an operational product the files overlap. DNORA automatically patches the data to use the first 6 hours, but uses longer lead times if files are missing to get a continues data set. See more in section "Working with forecast data"

**Steps needed to access data**: None. Available on the MET Norway thredds server. 

**Minimal example**: A minimal example to get WW3 wave spectra and MEPS winds for an area and save it as netcdf is:

.. code-block:: python

   import dnora as dn
   
   area = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   model = dn.modelrun.WW3_4km(area, year=2020, month=2, day=15)
   model.import_wind()
   model.import_spectra()
   
   exp = dn.export.Netcdf(model)
   exp.export_wind()
   exp.export_spectra()

**References**: TBA

Working with forecast data
=====================================

Folders, filenames and URL's in DNORA
=====================================

+++++
Using filenames and folders
+++++

The default file names and directories used in DNORA are defined in the ``export_default.yml``-file. These are names that we deemed to make most sense for specific models, but you can specify other names on export. The names are created based on tags, which can be used both in default names and user specified names.

Say we have the following setup:

.. code-block:: python

   grid = dn.grid.Grid(lon=(5,10), lat=(55,58), name='ExampleGrid')
   model = dn.modelrun.NORA3(grid, year=2020, month=3)
   model.import_spectra()
   model.import_wind()

   exp = dn.export.SWAN(model)

If we now run ``exp.export_wind()``, then the SWAN defaults are

   filename: ``'wind#WIND#GRID#T0_#T1'``

   extension: ``'asc'``

   folder: ``#GRID_SWAN``

This means that a folder will be created and the file will be written to: ``ExampleGrid_SWAN/windNORA3ExampleGrid20200301_20200330.asc`` (we use dense naming because SWAN has character limits in the input file)

If you want to export this to a folder that has a nested monthly structure, drop the grid name from the filename and change the extension to ``'.txt'``:

.. code-block:: python

   exp.export_wind(folder='#GRID_SWAN/#FT0', dateformat_folder='%Y/%m', filename='wind#WIND_#T0_#T1.txt')


This creates the nested folders and writes to: ``ExampleGrid_SWAN/2020/03/windNORA3_20200301_20200330.txt``

These tags are also used by the function that read data. E.g. the url to the NORA3 spectral files is a combination of ``'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m'`` and ``'SPC%Y%m%d00.nc'``-


+++++++
Explanation of tags
+++++++


The following tags referring to objects now mean (example where possible):

   ``'#GRID'``: Name of the grid object ``model.grid()`` (e.g. ``'ExampleGrid'``)

   ``'#WIND'``: Name of the wind forcing object ``model.wind()`` (e.g. ``'NORA3'``)

   ``'#SPECTRA'``: Name of the boundary spectra object ``model.spectra()`` (e.g. ``'NORA3'``)

   ``'#SPECTRA1D'``: Name of the omnidirectional spectra object ``model.spectra1d()`` 

   ``'#ICE'``: Name of the ice object ``model.ice()``

   ``'#WATERLEVEL'``: Name of the waterlevel object ``model.waterlevel()``

   ``'#CURRENT'``: Name of the current object ``model.current()``


The following tags refer to the properties of the setup:


   ``'#T0'``: Start time of the model run (e.g. ``'20200301T00'`` if given in combination with ``dateformat='%Y%m%dT%H'``) 

   ``'#T1'``: End time of the model run (e.g. ``'20200330T23'`` if given in combination with ``dateformat='%Y%m%dT%H'``)

   ``'#FT0'``: Same as ``'#T0'``, but for folders (e.g. ``'202003'`` if given in combination with ``dateformat_folder='%Y%m'``)

   ``'#FT1'``: Same as ``'#T1'``, but for folders (e.g. ``'202003'`` if given in combination with ``dateformat_folder='%Y%m'``)

   ``'#LON0'``: Eastern edge of object (e.g. ``'05.00'`` if given in combination with ``coord_format='02.2f'`` when exporting ``model.grid()``)

   ``'#LON1'``: Western edge of object (e.g. ``'10.00'`` if given in combination with ``coord_format='02.2f'`` when exporting ``model.grid()``)

   ``'#LAT0'``: Southern edge of object (e.g. ``'05.00'`` if given in combination with ``coord_format='02.2f'`` when exporting ``model.grid()``)

   ``'#LAT1'``: Northern edge of object (e.g. ``'10.00'`` if given in combination with ``coord_format='02.2f'`` when exporting ``model.grid()``)

   ``'#X0'``: Left edge of cartesian object. Format specifiel by ``cartesian_coord_format``

   ``'#X1'``: Right edge of cartesian object. Format specifiel by ``cartesian_coord_format``

   ``'#Y0'``: Lower edge of cartesian object. Format specifiel by ``cartesian_coord_format``

   ``'#Y1'``: Upper edge of cartesian object. Format specifiel by ``cartesian_coord_format``


There are a couple more that are used mostly for caching purposes. These are all forced to lowercase:

   ``'#DataType'``: The name of the datatype being exported, e.g. ``'wind'`` when exporting ``model.wind()``

   ``'#ObjectName'``: The name of the object being exported, e.g. ``'nora3'`` when exporting ``model.wind()`` in our example.

   ``'#ModelRun'``: The name of the ModelRun-object containing the data. Default to ``'dnoramodelrun'`` in our example.


Data sources and environmental variables
=============================================

+++++++
Data sources
+++++++

Since the same data might exist in several places, DNORA has formalized different sources to get data from. They are defined in ``dnora.type_manager.data_sources.DataSource``, but the user can refer to them using string if needed:

   ``DataSource.LOCAL``: A local file on e.g. your own computer

   ``DataSource.REMOTE``: A remote server, e.g. MET Norway's thredds-server or the copernicus database

   ``DataSource.CREATION``: Data that is being created on the fly, e.g. a "reader" that creates constant data

   ``DataSource.INTERAL``: Used internally in MET Norway, but can be used to have an alternative pre-defined path to ``LOCAL`` or ``REMOTE``

   ``DataSource.IMMUTABLE``: As ``INTERNAL``

   ``DataSource.UNDEFINED``: Used as a default source in abstract classes etc.

+++++++
Environmental variables
+++++++

DNORA can use enironmental variables to point to fixed locations. That way a folder name that is being used constantly doesn't need to be provided. The environmental variables can be valid for an enitre ``DataSource``, or specific for a ``DataSource`` and data type. A data type specific variable has priority over a general one. For example, to read boundary spectra from an internal source, DNORA first searches for the environmental variable ``DNORA_INTERAL_SPECTRA_PATH``. If that is not found, it uses ``DNORA_INTERNAL_PATH``. 

The enivronmental variables can also be used to specify where the different wave models are installed, e.g. by specifying ``DNORA_SWAN_PATH`` or ``DNORA_WW3_PATH``.

You can set the enivronmental variable in three ways:

1) ``export DNORA_INTERNAL_PATH=/path/to/internal/source``, preferrable in your ``.bashrc``
2) Set the environmental variable for a project folder only by creating a ``.env`` file with ``DNORA_INTERNAL_PATH=/path/to/internal/source``
3) Set the environmental variable for a single script: ``os.environ['DNORA_INTERNAL_PATH'] = '/path/to/internal/source'``


++++++++
Data sources for storing bathymetry
++++++++

Specifically we might want to store bathymetric data in one place, since they are typically being used across several project. For example if we want to use EMODNET data for the first time:

.. code-block:: python

   area = dn.grid.EMODNET(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   area.import_topo()

The needed tiles will be download into a folder ``'EMODNET/2022'`` in the folder where you run your script. If you then want to use it in another place, the reader won't find it. You can tell the reader to store it in a central place by specifying a folder:

.. code-block:: python

   area.import_topo(folder='/folder/for/bathymetric_data')

Then this needs to be given every time. The better way is to set an environmental variable:

.. code-block:: python

   export DNORA_LOCAL_GRID_PATH=/folder/for/bathymetric_data

Then all bethymetric data that is automatically be downloaded will be stored there. If you download data manually and place it there (e.g. Kartverket50m), then it will also be found automatically between projects.

If you have an internal common storage that is separate from your own computer, then you can also define a path to that:

.. code-block:: python

   export DNORA_INTERNAL_GRID_PATH=/lustre/projects/wave_modelling/bathymetry

You can access that storage by providing the folder explicitly, but by setting the environmental variable you can just switch the source:

.. code-block:: python

   # These two are now equivalent
   area.import_topo(folder='/lustre/projects/wave_modelling/bathymetry')
   area.import_topo(source='internal')


++++++++
Running models
++++++++

To run a model in dnora, use the ``executer`` object:

.. code-block:: python

   exe = dn.executer.SWAN(model)
   exe.write_input_file()
   exe.run_model()

Or for WAVEWATCH III that requires several input files and pre-processors:

.. code-block:: python

   exe = dn.executer.WW3(model)
   exe.write_grid_file()
   exe.write_wind_file()
   exe.write_spectra_file()
   exe.write_input_file()

   exe.run_grid()
   exe.run_wind()
   exe.run_spectra()
   exe.run_model()

To run the models within dnora, the paths where the models are installed, need to be defined somehow. Alternatives are:

i) In the ``~/.bashrc`` file:

.. code-block::

   export PATH=${PATH}:/home/user/Programs/swan
   export PATH=${PATH}:/home/user/Programs/swash
   export PATH=${PATH}:/home/user/Programs/HOS-ocean/bin
   export PATH=${PATH}:/home/user/Programs/REEF3D_xx/DIVEMesh/bin
   export PATH=${PATH}:/home/user/Programs/REEF3D_xx/REEF3D/bin

   source ~/.bashrc

ii) As an environmental variable:

.. code-block::

   DNORA_SWAN_PATH=/home/user/Programs/swan
   DNORA_SWASH_PATH=/home/user/Programs/swash

or alternatively in python:

.. code-block:: python

   import os
   os.environ["DNORA_SWAN_PATH"] = "/home/user/Programs/swan"

iii) As a keyword in the run method:

.. code-block:: python

   exe.run_model(model_folder='/home/user/Programs/swan')

++++++++
Dependencies
++++++++

1. Installation of **SWAN**. The latest SWAN version can be downloaded from https://sourceforge.net/projects/swanmodel/files/swan/ . The installation procedure can be found in: https://swanmodel.sourceforge.io/online_doc/swanimp/node4.html

2. Installation of **WAVEWATCH III**. The latest model version can be downloaded from https://github.com/NOAA-EMC/WW3 . The installation procedure can be found in: https://github.com/NOAA-EMC/WW3/wiki/Quick-Start

3. Installation of **SWASH**. The latest model version can be downloaded from https://sourceforge.net/projects/swash/ . The installation procedure can be found in: https://swash.sourceforge.io/online_doc/swashimp/node4.html

4. Installation of **HOS-ocean**. The latest HOS-ocean version can be downloaded from https://github.com/LHEEA/HOS-ocean . The installation procedure can be found in: https://github.com/LHEEA/HOS-ocean/wiki/Installation

5. Installation of **REEF3D**. The latest REEF3D and DIVEMesh versions can be downloaded from https://github.com/REEF3D . The installation procedure can be found in: https://reef3d.wordpress.com/installation/ . Note: Only REEF3D::FNPF has been tested in dnora.

6. Fimex is a the File Interpolation, Manipulation and EXtraction library for gridded geospatial data (more info in \url{httpshttps://wiki.met.no/fimex/start}). Fimex is used in DNORA for the preparation of forcing fields (wind). **In case of running the spectral model without wind forcing, the fimex installation can be neglected**.  A detailed installation procedure can be find in https://wiki.met.no/fimex/install. For a quick installation in Linux-Ubuntu, you can follow the steps: open the Synaptic Package Manager (to install it: sudo apt-get install synaptic
) --> Settings --> Repositories --> Other Software --> Add:  **ppa:met-norway/fimex** to your system's Software Sources (see https://launchpad.net/~met-norway/+archive/ubuntu/fimex). Then, search and mark for installation a latest version of fimex, and press apply installation. Check the installation (usually it is installed in /usr/bin/) by typing in command line: fimex or fimex-xxx where xxx is the version number. In case that only fimex-xxx works then add a symbolic link::

      cd /usr/bin
      sudo ln -s fimex-xxx fimex

[To be added]

.. toctree::
   :maxdepth: 2
   :caption: Contents:
