from pathlib import Path
import yaml
def read_defaults(filename: str, from_module: bool=False):
    if from_module:
        defaults_file = Path(__file__).parent.joinpath(Path(filename))
    else:
        defaults_file = Path(filename)
    with open(defaults_file, 'r') as file:
        default_values = yaml.safe_load(file)

    return default_values