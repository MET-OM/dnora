Setting a grid and importing bathymetrical data
====================================================


Setting area and resolution
+++++++++++++++++++++++++++++++++++

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



Available bathymetrical data
++++++++++++++++++++++++++++

EMODNET (``dnora.read.grid.EMODNET``)
---------------------------------------------

Covers Europe with approximately 115 m resolution. To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.EMODNET(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo() # No reader needs to be specified

The EMODNET reader automatically identifies the tiles needed to cover the area, and downloads them if they do not exist in the folder.

Reference: https://emodnet.ec.europa.eu/en/bathymetry


GEBCO (``dnora.read.grid.GEBCO``)
---------------------------------------------

Global coverage approximately 500 m resolution. To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(filename='name_of_downloaded_tile.nc') # No reader needs to be specified

The GEBCO reader is not as automated as the EMODNET reader, and basically wants to read custom tiles that have already been downloaded. You can download tiles here: https://download.gebco.net/ 

However, it does make use of an API provided by (https://github.com/cywhale/gebco). This means that if no filename is provided, then the relevant data is automatically fetched using the API. However, this doesn't work for very large areas. Recommended practice is therefore to download a faiirly large area once, and then make use of that several times.

The functionality will be improved in future versions.

References: https://gebco.net, https://github.com/cywhale/gebco


Karverket-50m (``dnora.read.grid.KartverketNo50m``)
---------------------------------------------------

High-resolution coverage for the Norwegian coast at a 50 m resolution provided by the Norwegian Mapping Authority (Kartverket). To use this reader you can also just set the grid a subclass:

.. code-block:: python

   grid = dn.grid.Kartverket(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(tile='B1008')

This reader is not automated. This means that you need to manually check which tile(s) you need and download them from Mapping Authorities web service, extract it and provide the tile to the reader.

References: https://kartkatalog.geonorge.no/metadata/dybdedata-terrengmodeller-50-meters-grid-landsdekkende/bbd687d0-d34f-4d95-9e60-27e330e0f76e


Netcdf (``dnora.read.generic.Netcdf``)
--------------------------------------

Reads generic Netcdf files with structured data. E.g.:

.. code-block:: python

   grid = dn.grid.Grid(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(dn.read.generic.Netcdf(), filename='some_random_topo.nc', folder='')

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.


PointNetcdf (``dnora.read.generic.PointNetcdf``)
-------------------------------------------------------------------------

Reads generic Netcdf files with unstructured data. This can be handy e.g. when you have donwloaded e.g. GEBCO data using the automatic API and wan't to save it for later:

.. code-block:: python

   grid = dn.grid.GEBCO(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo()
   grid.raw().ds().to_netcdf('GEBCO_Bergen.nc')

   # Re-import from file
   grid.import_topo(dn.read.generic.PointNetcdf(), filename='GEBCO_Bergen.nc', folder='')

Assumes that the file has a variable named ``'topo'`` or that it has the standard name ``'sea_floor_depth_below_sea_surface'``. Assumes depths are positive.

For any more regular use, it is recommended to build a reader for that specific product.


ConstantData (``dn.read.generic.ConstantData``)
-------------------------------------------------------------------------

To import a constant value for a grid you can use the Constant class:

.. code-block:: python

   grid = dn.grid.Constant(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo(topo=50)


Meshing and processing bathymetrical data
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

A bathymetry that has been imported from a source can be in any format (gridded, non-gridded, cartesian coordinates etc.). We therefore need to mesh the available data down to our grid specifications. 

This original data is in itself a grid object and is stored at ``grid.raw()``. For normal use cases, however, you don't need to access this object, since the meshing (wrapper around ``scipy.interpolate.griddata``) is done automatically by interpolation (either nearest neighbour or bi-linear). 

As an example, import EMODENET data for the Bergen area and mesh it down to a 500 m regualar grid:

.. code-block:: python

   grid = dn.grid.EMODNET(lat=(60.0, 60.6), lon=(4.4, 5.9), name="Bergen")
   grid.import_topo()
   grid.set_spacing(dm=500) 
   grid.mesh_grid(method='nearest') # method = 'nearest' [default] or 'linear'

The grid now has automatically defined a land-sea mask based on the depth values (``grid.land_mask()`` and ``grid.sea_mask()``).

**NB!** dnora has a system for keeping track of folders and possible data that has been downloaded automatically so it can be reused between projects. See section "Folders, filenames and URL's in DNORA" for more information.

If you want to modify the meshed data, there are some functions to do that. E.g. to set everything less than 1 m deep to land, and everything less than 2 m deep to 2 m:

.. code-block:: python

   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=1, to_land=True)
   grid.process_grid(dn.grid.process.SetMinDepth(), min_depth=2)



