from __future__ import annotations
from pathlib import Path
import yaml
from dotenv import load_dotenv
import os
import sys
from typing import TYPE_CHECKING

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

    load_dotenv(f"{sys.path[0]}/.env")
    value = os.getenv(f"DNORA_{data_source.name}_{obj_type.name}_PATH")
    if value is None:
        value = os.getenv(f"DNORA_{data_source.name}_PATH")
    return value
