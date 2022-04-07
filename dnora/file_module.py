from pathlib import Path
import yaml
from dataclasses import dataclass
import pandas as pd
import re

def add_prefix(filename: str, prefix: str) -> str:
    """Adds a prefix to a filename, e.g. FileName.txt -> new_FileName.txt"""
    if (not prefix == '') and (not prefix[-1] == '_'):
        return f"{prefix}_{filename}"
    else:
        return f"{prefix}{filename}"

def add_suffix(filename: str, suffix: str) -> str:
    """Adds a suffix to a filename, e.g. FileName.txt -> FileName_new.txt"""
    if (not suffix == '') and (not suffix[0] == '_'):
        suffix = f"_{suffix}"

    filename_list = filename.split('.')

    if len(filename_list) == 1:
        return filename+suffix
    else:
        return '.'.join(filename_list[0:-1]) + suffix + '.' + filename_list[-1]

def replace_times(filename: str, dateformat: str, times: list) -> str:
    """Substitutes the strings #T0, #T1, #T2... etc. in filname with time
    stamps from a list of times, using the format in dateformat.

    e.g. #T0_file.txt, ['2020-05-04 18:00'], %Y%m%d%H%M -> 202005041800_file.txt
    """


    for ct, t in enumerate(times):
        filename = re.sub(f"#T{ct}", pd.Timestamp(t).strftime(dateformat), filename)

    return filename

def replace_lonlat(filename: str, lon: float, lat: float) -> str:
    """Substitutes the strings #Lon, #Lat in filename with values of lon and
    lat.

    e.g. #Lon_#Lat_file.txt, 8.0, 60.05 -> 08.0000000_60.05000000_file.txt
    """

    filename = re.sub("#Lon", f"{lon:010.7f}", filename)
    filename = re.sub("#Lat", f"{lat:010.7f}", filename)

    return filestring

def replace_objects(filename: str, objects: list) -> str:
    """Substitutes the strings #{Object} in filename with the name given to
    the object.

    e.g. #Grid_#Forcing_file.txt, [Grid(..., name="Sula"), Forcing(..., name='NORA3')]
        -> Sula_NORA3_file.txt
    """

    for object in objects:
        if object is not None:
            obj_str = type(object).__name__
            obj_name = object.name()
            filename = re.sub(f"#{obj_str}", obj_name, filename)

    return filename

def clean(filename: str, list_of_placeholders: list[str]) -> str:
    """ Cleans out the file name from possible used placeholders, e.g. #Grid
    as given in the list.

    Also removes multiple underscores '___'
    """

    for s in list_of_placeholders:
            filename = re.sub(s, '', filename)

    filename = re.sub("_{2,10}", '_', filename)

    return filename

def get_default_value(key: str, module:str, primary: dict, fallback: dict):
    """Get a key (e.g. folder) from the defaults list.

    1) Tries Model+module specific value (e.g. SWAN-wnd-folder)

    2) Tries Model specifc values (e.g. SWAN-folder)

    3) Returns ModelRun defaults (e.g. ModulRun-wnd-folder)
    """

    if module not in fallback.keys():
        raise ValueError(f'Default values not defined for module {module}!')
    fallback_filename = fallback[module].get(key)

    # Try module specific filename is module settings defined
    module_filename = None
    if primary.get(module) is not None:
        module_filename = primary[module].get(key) or module_filename

    # If filename not defined for specific module, try Model specific name
    module_filename = module_filename or primary.get(key)

    return module_filename or fallback_filename

def add_folder_to_filename(filename: str, folder: str) -> str:
    return str(Path(folder).joinpath(filename))

def split_filepath(filepath: str) -> tuple[str, str]:
    folder = str(Path(filepath).parent)
    filename = Path(filepath).name
    return filename, folder

@dataclass
class FileNames:
    format: str
    clean_names: bool
    list_of_objects: list
    _filename: str
    _folder: str
    _dateformat: str
    extension: str
    module: str = None
    dnora_obj: str = None

    def __post_init__(self):
        defaults_file = Path(__file__).parent.joinpath(Path('defaults.yml'))
        with open(defaults_file, 'r') as file:
          self._defaults = yaml.safe_load(file)
        self.fallback = self._defaults['ModelRun']
        self.primary = self._defaults[self.format]
        self.placeholders = self._defaults['list_of_placeholders']
        if self.module is None and self.dnora_obj is None:
            raise ValueError('Provide either module or object!')
        if self.module is None:
            self.module = self._defaults['list_of_modules'][self.dnora_obj]

    def dateformat(self) -> str:
        return self._dateformat or get_default_value('dateformat', self.module, self.primary, self.fallback)

    def filename(self) -> str:
        filename = self._filename or get_default_value('filename', self.module, self.primary, self.fallback)
        filename = self.replace_placeholders(filename, self.dateformat())
        return Path(filename).with_suffix(f'.{self.extension}')

    def folder(self) -> str:
        folder = self._folder or get_default_value('folder', self.module, self.primary, self.fallback)
        return Path(self.replace_placeholders(folder, self.dateformat()))

    def filepath(self) -> str:
        return add_folder_to_filename(self.filename(), self.folder())

    def create_folder(self) -> None:
        folder = Path(self.folder())
        if not folder.is_dir():
            msg.plain(f"Creating folder {str(folder)}")
            folder.mkdir()

    def replace_placeholders(self, unclean_string: str, dateformat: str) -> str:
        unclean_string = replace_objects(unclean_string, self.list_of_objects)

        start_time = self.list_of_objects[0].start_time
        end_time = self.list_of_objects[0].end_time
        clean_string = replace_times(unclean_string, dateformat, [start_time, end_time])

        if self.clean_names:
            clean_string = clean(clean_string, list_of_placeholders=self.placeholders)

        return clean_string
