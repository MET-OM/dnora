import numpy as np
import pandas as pd
import xarray as xr
from abc import ABC, abstractmethod
import netCDF4
from dnora2 import msg, spec, wnd
from dnora2.aux import distance_2points, min_distance, day_list, month_list, create_time_stamps
from copy import copy





