from dnora.utils.time import create_time_stamps
from dataclasses import dataclass
from dnora.utils.time import create_monthly_stamps, create_yearly_stamps

from typing import Union

@dataclass
class FileStructure:
    """stride:
    int in hours (e.g. 24 for daily files)
    'month' for monthly files
    'year' for yearly files
    None if all data is in one file (but points can be spread across different files)
    """

    stride: Union[int, str, None] = None
    hours_per_file: int = None
    lead_time: int = 0
    last_file: str = ""
    offset: int = 0
    tile: str = None

    def create_time_stamps(self, start_time: str, end_time: str, last_file: str = ''):
        if self.stride is None:
            start_times, end_times = [start_time], [end_time]
            file_times = start_times
        elif self.stride == "month":
            start_times, end_times = create_monthly_stamps(start_time, end_time)
            file_times = start_times
        elif self.stride == "year":
            start_times, end_times = create_yearly_stamps(start_time, end_time)
            file_times = start_times
        else:
            start_times, end_times, file_times = create_time_stamps(
                start_time,
                end_time,
                stride=self.stride,
                hours_per_file=self.hours_per_file,
                last_file=(last_file or self.last_file),
                lead_time=self.lead_time,
                offset=self.offset,
            )

        return start_times, end_times, file_times
