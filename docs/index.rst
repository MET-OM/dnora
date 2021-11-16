.. dnora documentation master file, created by
   sphinx-quickstart on Sun Nov 15 14:18:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dnora's documentation!
=====================================
**dnora** is a Python package for dynamical downscaling of NORA wave hindcast using the spectral wave models SWAN or WAVEWATCH III. 

The package contains functions that: 
  * create a high-resolution grid using open-access bathymetry/topography datasets (e.g., Emodnet 2018 bathymetry)
  * prepare the boundary conditions (NORA3, WAM4-operational wave model at MET Norway, available in https://thredds.met.no/thredds/catalog.html) for the spectral model
  * prepare the wind (NORA3, WAM4) forcing for the spectral model 
  * create input parameter files (e.g., .swn, .inp) for the spectral model
  * run the spectral model


Grid generation
=====================================
This section document the grd-module.  The main idea is that the aGrid-object is created,and  a  fixed  set  of  methods  are  used  to  import  a  topography,  mesh  it  down  to  a  grid,  orfilter the data.  The functionality of these methods are controlled by passing them an object(which is a callable function).Adding e.g.  a topography source thus means adding a newTopoFetcher-class that canthen me passed to the Grid-object’s.importtopo()-method.  A similar principle goes foradding meshing of filtering functionalities.

The grid object is initialized with the following command::

   example_grid = grd.Grid(lon_min, lon_max, lat_min, lat_max, name = ’GridName’)

.. code-block:: rst

Use ``print(grid)`` to print out the status of te object.

A desired grid spacing can be set either by providing a desired grid spacing in metres (ap-proximate) or defining the amounts of grid points (exact)::

   grid.set_spacing(dm = 250) # Set spacing to around 250 metres
   grid.set_spacing(nx = 291, ny = 249) # Create 291 (lon) x 249 (lat) grid points

.. code-block:: rst

Both  of  these  options  will  convert  the  input  to  the  native  resolution  in  longitude  andlatitude.  These can, of course, also be set directly by::
   
   grid.set_spacing(dlon = 0.0048780, dlat = 0.0022573)
   
.. code-block:: rst

In this case ``dlon`` and ``dlat`` are not exact.  If an exact resolution needs to be forced, the ``floatingedge-option`` can be used, e.g.,::

   grid.set_spacing(dlon = 1/205, dlat = 1/443, floating_edge = True)
   
.. code-block:: rst
   
This will enforce the resolution and instead change the initially set ``lon_max`` and ``lat_max`` slightly (if needed).The longitude and latitude vectors, and the name, of the grid can be accessed by::

   grid.lon()
   grid.lat()
   grid.name()
   
.. code-block:: rst

Run the spectral model
=====================================

Run SWAN::

   run.run_SWAN(input_file_name,swan_directory=swan_directory)

.. code-block:: rst


.. toctree::
   :maxdepth: 2
   :caption: Contents:



