




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
