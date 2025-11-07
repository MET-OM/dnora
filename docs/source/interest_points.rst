Setting interest (e.g. boundary) points for a grid
====================================================

Meaning of "interest points"
++++++++++++++++++++++++++++++++++++++

The interest points are logical masks that can be used to define points of particular interest in the grid. They serve two main purposes:

a. Determine import of data
b. Write information for model grid and definition files

There are three different interest points that can be set for the grid.

1. Boundary points
-------------------

These points represent the boundary points in a wave model grid. They are used for:

a. Determining for which points to import wave spectra when a ``NearestGridPoint`` ``PointPicker`` is used.
b. Used to write the proper grid and/or input file for the model so the boundary conditions propagate properly into the grid.


2. WaveSeries points
---------------------

These points are much like the boundary points, but for when importing waveseries (Hs, Tp etc.) data instead of wave spectra. Note, that if you want to set the boundary conditions by importing integrated parameters and converting them to spectra, you need to set the waveseries points for import, and the boundary points for writing the grid etc. files.

3. Ouput points
---------------------

These have no consequency for importing data, but they can be used to determine the points where the wave model should output wave spectra.

Accessing and setting interest points
++++++++++++++++++++++++++++++++++++++

Interest points are stored in boolean masks. To get the mask of boundary points for a grid:

.. code-block:: python

   mask = grid.boundary_mask()

To get the corresponding longitude and latitude values for the set points:

.. code-block:: python

   bnd_lon, bnd_lat = grid.boundary_points()

Similar methods exist for ``waveseries_mask`` and ``output_mask``

Interest points can be set just like the land-sea mask:

.. code-block:: python

   grid.set_boundary_mask(mask)

However, this is not very convenient, since creating the mask is the difficult part. You can therefore use another method which takes a class that will constuct the mask for you using different logics. The most common case is setting boundary points for the edges of a grid:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

This defines the North, West and South edges of the grids as boundary points.


Available logics for setting points
+++++++++++++++++++++++++++++++++++

These same classes can be used to set all of the different types of interest points.

Whole edges (``dnora.grid.mask.Edges``)
----------------------------------------

Set all wet points along a certain edge as boundary points:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Edges(edges=['N', 'W', 'S'])

This is probably what you want for spectral models.

Mid points (``dnora.grid.mask.MidPoint``)
----------------------------------------

Same as Edges, but only set the middle point:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.MidPoint(edges=['N', 'W'])

This might be needed in some phase resolving models.

A list of points (``dnora.grid.mask.LonLat`` or ``dnora.grid.mask.XY``)
-------------------------------------------------------------------------

Give a list of longitudes and latitudes (or x- and y-coordinates if grid is cartesian) and set the nearest grid point as a boundary point:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.LonLat(lon=(5.0, 5.1), lat=(60.3, 60.3)))

This is probably not used for setting bounadries, but this same function can be used to set a mask of e.g. spectral output points.


All points (``dnora.grid.mask.All``)
----------------------------------------

To set all points as boundary points:

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.All())

No points (``dnora.grid.mask.Clear``)
----------------------------------------

To clear all boundary points

.. code-block:: python

   grid.set_boundary_points(dn.grid.mask.Clear())