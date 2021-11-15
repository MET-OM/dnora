.. dnora documentation master file, created by
   sphinx-quickstart on Sun Nov 15 14:18:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dnora's documentation!
=====================================
What is dnora? 
dnora is a software for dynamical downscaling of NORA3 wave hindcast.

dnora can:

* create high resolution grid from open access Emodnet 2018 bathymetry within the NORA3 domain
* create wind forcing and 2d spectra at the boundaries based on NORA3 (hindcast) and WAM4 (forecast)
* create parameter file(s) for the spectral model (SWAN and WW3)
* run the spectral model

Run SWAN::

   run.run_SWAN(input_file_name,swan_directory=swan_directory)

.. code-block:: rst


.. toctree::
   :maxdepth: 2
   :caption: Contents:


