# dnora
# ![logo](https://user-images.githubusercontent.com/67804784/145466261-f50dbc27-f242-4db0-8d99-e23d0bd0dbbc.png)


Welcome to dnora version 2!


# Quick Installation

To install the latest version of dnora with pip:

```shell
$ pip install dnora 
```

The latest version of the old 'v1' dnora is 1.4.3, and can be installed with:

```shell
$ pip install dnora==1.4.3 
```

You can find the code of the origin dnora v1 in the branch 'v1' on GitHub. Please not that it is no longer mainained.


# What is dnora? 

dnora is a software for dynamical downscaling of wave products i.e., NORA3 wave hindcast and WW3 wave forecast from [Norwegian Meteorological Institute](https://www.met.no/) and ERA5 from [ECMWF](https://www.ecmwf.int/).

# Documentation

Please  [documentation](https://dnora.readthedocs.io/en/latest/) for dnora is still ongoing, but can be found at:

Old examples for dnora-v1 on [youtube](https://youtu.be/pTmjBnsXNz8) 

# Example output

Example of downscaling NORA3 using dnora/SWAN offshore Tromsø, Norway:
![dnora](https://user-images.githubusercontent.com/67804784/147151236-b9ef920c-34a2-4da0-9877-6241723eff80.gif)

Example of long-crested wave propagation at Stad, Norway, using dnora/SWASH: 
![eta](https://user-images.githubusercontent.com/67804784/160290851-ca743601-2ac7-48b5-be52-da3ec8c31e13.gif)


dnora contains the code of the function read_sms_mesh from the PyFVCOM package (distributed under the MIT License). The original code can be found at https://github.com/pwcazenave/pyfvcom

