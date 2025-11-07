Basic example
=============================================

dnora was originally created to easily downscale the NORA3 hindcast. While it can now do much more, here is how you do that. The code runs a 200 m SWAN model nested inside a 500 m SWAN model for one day using EMODNET topography, NORA3 wind forcing, and NORA3 boundary spectra:

.. code-block:: python

   import dnora as dn
   
   sula500 = dn.grid.EMODNET(lon=(5.21, 6.66), lat=(62.25, 62.89), name="sula500")
   sula500.set_spacing(dm=500)
   sula500.import_topo()
   sula500.mesh_grid()
   # Propagate boundary spectra for North and West edges
   sula500.set_boundary_points(dn.grid.mask.Edges(edges=["N", "W"]))

   sula200 = dn.grid.EMODNET(lon=(5.5, 6.66), lat=(62.25, 62.6), name="sula200")
   sula200.set_spacing(dm=200)
   sula200.import_topo()
   sula200.mesh_grid()
   
   model = dn.modelrun.NORA3(sula500, year=2020, month=2, day=15)
   # sula500 will output boundary conditions for sula200
   model.set_nested_grid(sula200) 
   # use same wind for both runs
   # model.nest().import_wind(...) to import another wind for sula200
   model.import_wind() 
   # Only for sula500, since sula200 uses results from sula500
   model.import_spectra() 
   
   # Exporter will automatically export files and keep track of all nesting
   exp = dn.export.SWAN(model)
   exp.export_grid()
   exp.export_wind()
   exp.export_spectra()
   
   # Executer will execute all the nested models in the right order
   exe = dn.executer.SWAN(model)
   exe.write_input_file()
   exe.run_model()


Before running the example you have to have SWAN installed, since it is not a part of the DNORA package.

The basic workflow in dnora scripts follow the same logic
  * Define an area you are working with by creating a Grid-object
  * Define a time period you are working with by creating a ModelRun-object
  * Import the data from the source you want
  * Define an exporter to export the data in the format you want
  * Define an executer to write input files and run the model

Specifically, the import from different sources (e.g. MET Norway, ECMWF) and the writing data in differen formats (e.g. SWAN, WAVEWATCH III) are separated, and you can always use any combination you want, while dnora takes care of making sure the data is in the right format for the model (e.g. spectral conventions).


