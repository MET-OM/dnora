#!/usr/bin/env python

import os
import setuptools

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'dnora/version.py')).read())

setuptools.setup(
    name        = 'dnora',
    description = 'dnora - a software for dynamical downscaling of wave products i.e., NORA3 wave hindcast and WAM4 wave forecast from MET Norway',
    author      = 'Konstantinos Christakos & Jan-Victor Björkqvist / MET Norway',
    url         = 'https://github.com/KonstantinChri/dnora',
    download_url = 'https://github.com/KonstantinChri/dnora',
    version = __version__,
    license = 'GPLv2',
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'netCDF4'
    ],
    packages = setuptools.find_packages(),
    include_package_data = True,
    setup_requires = ['setuptools_scm'],
    tests_require = ['pytest'],
    scripts = []
)
