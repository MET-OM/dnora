from __future__ import annotations
from pathlib import Path
import yaml
from dotenv import load_dotenv
import os
import sys
from typing import TYPE_CHECKING
from dnora import msg
if TYPE_CHECKING:
    from type_manager.dnora_types import DnoraDataType, DataSource


def read_defaults(filename: str, from_module: bool = False):
    if from_module:
        defaults_file = Path(__file__).parent.joinpath(Path(filename))
    else:
        defaults_file = Path(filename)
    with open(defaults_file, "r") as file:
        default_values = yaml.safe_load(file)

    return default_values


def read_environment_variable(obj_type: DnoraDataType, data_source: DataSource) -> str:
    """Reads an environmental variable:
    1) Loads possible .env-file
    2) Reads e.g. DNORA_LOCAL_GRID_PATH
    3) If 2) fails, reads e.g. DNORA_LOCAL_PATH"""
    load_dotenv(f"{sys.path[0]}/.env")
    long_var = f"DNORA_{data_source.name}_{obj_type.name}_PATH"
    short_var = f"DNORA_{data_source.name}_PATH"
    
    used_var = long_var
    value = os.getenv(long_var)

    if value is None:
        used_var = short_var
        value = os.getenv(short_var)
    
    if value is not None:
        value = os.path.expanduser(value)
        msg.plain(f"Reading {used_var}={value}")        
    return value
