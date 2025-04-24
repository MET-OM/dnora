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

   $ git clone https://github.com/KonstantinChri/dnora.git
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

Basic workflow
=============================================
The basic workflow in dnora scripts follow the same logic
  * Define an area you are working with by creating a Grid-object
  * Define a time period you are working with by creating a ModelRun-object
  * Import the data from the source you want
  * Define an exporter to export the data in the format you want
  * Define an executer to write input files and run the model

Specifically, the import from different sources (e.g. MET Norway, ECMWF) and the writing data in differen formats (e.g. SWAN, WAVEWATCH III) are separated, and you can always use any combination you want, while dnora takes care of making sure the data is in the right format for the model (e.g. spectral conventions).

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


Creating a Grid-object
=====================================
In the previous example we only used the Grid-object to define an area. We will now create a proper grid. The grid object is initialized with the following command::

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

Importing and meshing a topography
=====================================

The grid created above is empty. To get bathymetrical data we need to import it. Trivially we can just test it by "importing" a bathymetry with constant 50 meter depth::

   grid.import_topo(dn.read.generic.ConstantData(), topo=50)

.. code-block:: rst

The "imported" raw bathymetry can be found at ``grid.raw()`` if you want to look at it (the corresponding xarray Dataset is ``grid.raw().ds()``), but what you really want to do is mesh this raw topography to the grid spacing you defined::

   grid.mesh_grid()

.. code-block:: rst

This interpolates the "raw" topography, which can be either gridded or non-gridded, to the spacing you defined for your grid object. The meshing is a wrapper around ``scipy.interpolate.griddata`` using nearest neighbour interpolation. Change this by providing e.g. ``method='linear'``.

*NB!* DNORA has a system for keeping track of folders and possible data that has been downloaded automatically so it can be reused between projects. See section "Folders, filenames and URL's in DNORA" for more information.

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

Setting boundary points in a grid
=====================================

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

Creating a ModelRun-object
=====================================

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



The grid data can now be exported in a certain format using a ``GridWriter``. To export in WAVEWATCH III format::

   model.export_grid(grd.write.WW3())

.. code-block:: rst

Boundary and Forcing data can be read using ``BoundaryReaders`` and ``ForcingReaders``. To read in boundary spectra and wind forcing from the MET Norway NORA3 hindcast, use::

   model.import_boundary(bnd.read_metno.NORA3(), point_picker=bnd.pick.Area())
   model.import_forcing(wnd.read_metno.NORA3())

.. code-block:: rst

Here, the ``PointPicker`` object defines how spectra are chosen from the database. In WW3, we take all spectra around the grid area, and let WW3 interpolate. For SWAN, we would want to use::

   model.import_boundary(bnd.read_metno.NORA3(), point_picker=bnd.pick.NearestGridPoint())

.. code-block:: rst

to connect each defined boundary point to a certain spectra (even though we can get duplicates).

to write the boundary spectra in WAVEWATCH III format and wind forcing in SWAN format, use::

   model.export_boundary(bnd.write.WW3())
   model.export_forcing(wnd.write.SWAN())

.. code-block:: rst

The spectral convention is defined in the ``BoundaryReader`` and ``BoundaryWriter``, and the ``ModelRun``-object automatically takes care of convention changes (if needed).

**NB!** The WW3 convention here is thath of the WW3 *output* files, i.e. directional vector looks like a mathematical convention, but it is actually oceanic. This is in line with the bounc.ftn file used in the develop-branch of WAVEWATCH III.

Creating templates
=====================================

Several features that are typically used together can be packaged as a "template" by creating a subclass of the ``ModelRun`` object. These are defined in the ``mdl/models.py``-file. For example, a ``WW3``-template is defined as::

   class WW3(ModelRun):
       def _get_boundary_writer(self):
           return bnd.write.WW3()

       def _get_forcing_writer(self):
           return wnd.write.WW3()

       def _get_point_picker(self):
           return bnd.pick.Area()

       def _get_grid_writer(self):
           return grd.write.WW3()

.. code-block:: rst

These defaults can be used by::

   model = mdl.WW3(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

   model.import_boundary(bnd.read_metno.NORA3()) # PointPicker defined in template
   model.export_boundary() # BoundaryWriter defined in template

.. code-block:: rst

Further subclasses can also be defied. For example to have a ``ModelRun``-object that uses WW3 conventions and gets the forcing data from the NORA3-hindcast, a ``WW3_NORA3``-template is defined using the above ``WW3``-template::

   class WW3_NORA3(WW3):
       def _get_boundary_reader(self):
           return bnd.read_metno.NORA3()

       def _get_forcing_reader(self):
           return wnd.read_metno.NORA3()

.. code-block:: rst

The above importing and exporting of NORA3 boundary is now simplified to::

   model = mdl.WW3_NORA3(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

   model.import_boundary() # BoundaryReader and PointPicker defined in template
   model.export_boundary() # BoundaryWriter defined in template
.. code-block:: rst

The defaults of the templates can always be overridden by explicitly providing an object to the method. For example, the following code import WAM 4km boundary spectra, not NORA3 spectra::

   model = mdl.WW3_NORA3(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

   model.import_boundary(bnd.read_metno.WAM4km()) # Override BoundaryReader but use template PointPicker
   model.export_boundary() # BoundaryWriter defined in template
.. code-block:: rst



Run the spectral model
=====================================

This fucntionality is at the moment only available for SWAN and SWASH. To run the model automatically we first need to generate an input file::

   model.write_input_file(inp.SWAN())

.. code-block:: rst

This generates an input file based on the grid, boundary and forcing data that is available in the object. After that, the model can be automatically ran by::


   model.run_model(run.SWAN())

.. code-block:: rst

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
