.. dnora documentation master file, created by
   sphinx-quickstart on Sun Nov 15 14:18:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. toctree::
   :hidden:

   self

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

Contents
===============

.. toctree::
   :maxdepth: 3
   :glob:

   install
   basic_example
   grid
   interest_points
   importing_data
   rest
   technical_info




