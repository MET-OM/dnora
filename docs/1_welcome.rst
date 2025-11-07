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
