# dnora
# ![logo](https://user-images.githubusercontent.com/67804784/145466261-f50dbc27-f242-4db0-8d99-e23d0bd0dbbc.png)

What is dnora? 

dnora is a software for dynamical downscaling of wave products i.e., NORA3 wave hindcast and WW3 wave forecast from [Norwegian Meteorological Institute](https://www.met.no/) and ERA5 from [ECMWF](https://www.ecmwf.int/).

[Documentation and installation instructions](https://dnora.readthedocs.io/en/latest/)

[Quick Installation Procedure for DNORA and wave models](https://docs.google.com/document/d/1VzRDkQk5pqq2rkNKa_woC3ehJ9Dn2pcRNy4LtrcOQXY/edit?usp=sharing)

[Excercises using DNORA/SWAN](https://docs.google.com/document/d/1FlIE6ByyXRF7QOSXOdjFaLnX4q3P6YzZq_4YdVEjATw/edit?usp=sharing)

Examples in [youtube](https://youtu.be/pTmjBnsXNz8) 

Example of downscaling NORA3 using dnora/SWAN offshore Tromsø, Norway:
![dnora](https://user-images.githubusercontent.com/67804784/147151236-b9ef920c-34a2-4da0-9877-6241723eff80.gif)

Example of long-crested wave propagation at Stad, Norway, using dnora/SWASH: 
![eta](https://user-images.githubusercontent.com/67804784/160290851-ca743601-2ac7-48b5-be52-da3ec8c31e13.gif)


dnora contains the code of the function read_sms_mesh from the PyFVCOM package (distributed under the MIT License). The original code can be found at https://github.com/pwcazenave/pyfvcom

