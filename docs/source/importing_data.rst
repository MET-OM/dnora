Importing data
=====================================

To import data ou need to create a ``ModelRun``-object. This object is always defined for a certain grid and a certain time and contains all forcing and boundary data.

Setting a time period
+++++++++++++++++++++++

The ``Grid``-object is used to define the area, so what is left is to define the time period.

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, start_time='2018-08-25T00:00', end_time='2018-08-25T01:00')

If you need to prepare forcing files for a long hindcast, it might pay of to use e.g. monthly files. Then a good possibility can be to define the time with the alternative format:

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, year=2020, month=5)

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-05-01T23:00'``. If you need the first hour of the next month for e.g. hotstart-purposes, then:

.. code-block:: python

   model = dn.modelrun.ModelRun(grid, year=2020, month=5, hotstart_hour=True)

This is equivalent to using a ``start_time='2020-05-01T00:00'`` and ``end_time='2020-06-01T00:00'``



Importing gridded data
+++++++++++++++++++++++

The bathymetrical data has already been imported to the ``Grid``-object. Other gridded data will be imported to the ``ModelRun``-object. To e.g. import NORA3 wind forcing that covers the area and time period:

.. code-block:: python

   model.import_wind(dn.read.wind.metno.NORA3())

In practice, you can use often use subclasses (defined in ``dnora.modelrun.templates``) of the ``ModelRun``-object and then you don't have to specify the reader. For example, if you just want to work with NORA3 data, then the above code simplifies to:

.. code-block:: python

   model = dn.modelrun.NORA3(grid, year=2020, month=5)
   model.import_forcing()


For gridded data (such as wind) the import-method automatically gets all data that covers the grid you provided, but expands it by 20% to make sure the edges are covered (needed for interpolation etc.). If the imported area is too small or too large, this can be controlled by a keyword:

.. code-block:: python

   model.import_forcing(expansion_factor=1.3)

will expand the are by 30%. The default value here is 1.2.

Other exampeles of gridded data area ice, current, ocean, and athmosphere. For an overview of what the different data types consists of, please see the section on data types.

Importing non-gridded data
+++++++++++++++++++++++++++

For non-gridded data we have some more options, since we might want different ways to choose e.g. boundary spectra. This is technically done by separating the logic of choosing the points to import and the actual import of data. Usually dnora takes care of this automatically, and the logic is the following:

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

To speed up the process, DNORA actually limits the search for the nearest spectrum to an area expanded by 3 degrees in latitude and 6 degrees in longitude. This is mainly done to avoid having to calculate distances to points very far away, since fast cartesian calculations are not available if the database contains points above 84 degrees latitude. This is not visible to the user. 
The user can, however, choose to discard spectra that are too far from the boundary point. For example, if we rather not import any spectrum at all if the closest one is over 7 km away, then:

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

Also here the ``max_dist`` keyword can be used. 

Note, that since dnora uses a standard format for 2D-spectra, it can automatically integrate the spectra to omni-directional spectra and to some basic integrated parameters. E.g. To get the omnidirectional spectra of the NORA3 2D-spectra from above:

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

Non-gridded data that are imported using this logic are spectra, spectra1d and waveseries. For an overview of what the different data types consists of, please see the section on data types.