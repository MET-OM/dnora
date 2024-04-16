from dnora import msg
from collections.abc import Iterable
import pandas as pd
import numpy as np


def print_constant_values(
    data_dict: dict[str, np.ndarray], obj_type, time=None
) -> None:
    """Prints out what constant values have been set to the data variables"""
    for key, value in data_dict.items():
        if time is not None:
            val_vec, time_vec = decode_constant_array(value, time)
            for v, t in zip(val_vec, time_vec):
                msg.plain(f"Setting {key}={v} for {obj_type.name} from {t}")
        else:
            msg.plain(f"Setting {key}={value} for {obj_type.name}")


def create_constant_array(val, time, obj_size, time_vec):
    """Creates a constnat array based on an iterable value and iterable start times"""
    if not isinstance(val, Iterable):
        return np.full(obj_size, val)
    else:

        val_array = np.full(obj_size, 0)
        assert len(val) == len(time)
    time = pd.to_datetime(time)
    for n, t in enumerate(time):
        if t < time_vec[0] or t > time_vec[-1]:
            msg.warning(
                f"Given time {t} is outside ModelRun time range {time_vec[0].strftime('%Y-%m-%d %H:%M')} - {time_vec[-1].strftime('%Y-%m-%d %H:%M')}!"
            )
        mask = time_vec >= t
        val_array[mask, ...] = val[n]
    return val_array


def decode_constant_array(val_array, time_vec):
    # Assume time is first

    mean_var = np.mean(val_array, axis=tuple(range(1, len(val_array.shape))))
    inds = np.where(np.diff(mean_var) > 0.0001)[0]

    if inds.size == 0:  # All constant values
        return np.atleast_1d(mean_var[0]), pd.to_datetime(np.atleast_1d(time_vec[0]))

    val = np.zeros(len(inds) + 1)
    time = []
    time.append(time_vec[0])

    for n, ind in enumerate(inds):
        val[n] = mean_var[ind]
        time.append(time_vec[ind + 1])
    val[n + 1] = mean_var[ind + 1]

    return val, pd.to_datetime(time)
