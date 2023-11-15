from pathlib import Path
import yaml
from dotenv import load_dotenv
import os


def read_defaults(filename: str, from_module: bool = False):
    if from_module:
        defaults_file = Path(__file__).parent.joinpath(Path(filename))
    else:
        defaults_file = Path(filename)
    with open(defaults_file, "r") as file:
        default_values = yaml.safe_load(file)

    return default_values


def data_sources(data_source: str) -> str:
    data_source = data_source.lower()
    if data_source.lower() not in ["local", "internal"]:
        raise KeyError(
            f"data_source should be 'local' or 'internal', not {data_source}!"
        )
    load_dotenv()
    defaults = read_defaults("data_sources.yml", from_module=True)

    folder = os.getenv(f"DNORA_{data_source.upper()}") or defaults[data_source]

    if data_source == "internal" and not folder:
        raise Exception(
            "You haven't set an internal location. Set DNORA_INTERNAL in your Linux environment or .env file!"
        )

    return folder
