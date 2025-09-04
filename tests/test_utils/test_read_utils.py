from dnora.read.ds_read_functions import data_left_to_try_with


def test_data_left_to_try_with():
    hours_per_file = 12
    file_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 06:00:00",
        "2020-01-01 12:00:00",
        "2020-01-01 18:00:00",
    ]
    start_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 06:00:00",
        "2020-01-01 12:00:00",
        "2020-01-01 18:00:00",
    ]
    end_times = [
        "2020-01-01 05:00:00",
        "2020-01-01 11:00:00",
        "2020-01-01 17:00:00",
        "2020-01-01 23:00:00",
    ]
    n = 0
    ct = 0
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    # First file needs to exist!
    n = 0
    ct = 1
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    # Second file is missing byt we can go back one file since they contain 12 hours
    n = 1
    ct = 0
    assert data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    # Third and second file is missing
    # We can NOT go back TWO files since they contain only 12 hours
    n = 2
    ct = 2
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])


def test_data_left_to_try_with_no_overlap():
    """If we have no overlap (e.g. in a hindcast) then we should always get false if a file is missing"""
    hours_per_file = 6
    file_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 06:00:00",
        "2020-01-01 12:00:00",
        "2020-01-01 18:00:00",
    ]
    start_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 06:00:00",
        "2020-01-01 12:00:00",
        "2020-01-01 18:00:00",
    ]
    end_times = [
        "2020-01-01 05:00:00",
        "2020-01-01 11:00:00",
        "2020-01-01 17:00:00",
        "2020-01-01 23:00:00",
    ]
    n = 0
    ct = 0
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 0
    ct = 1
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 1
    ct = 0
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 2
    ct = 2
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])


def test_data_left_to_try_with_no_overlap_1h():
    """If we have no overlap (e.g. in a hindcast) then we should always get false if a file is missing"""
    hours_per_file = 1
    file_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 01:00:00",
        "2020-01-01 02:00:00",
        "2020-01-01 03:00:00",
    ]
    start_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 01:00:00",
        "2020-01-01 02:00:00",
        "2020-01-01 03:00:00",
    ]
    end_times = [
        "2020-01-01 00:00:00",
        "2020-01-01 01:00:00",
        "2020-01-01 02:00:00",
        "2020-01-01 03:00:00",
    ]
    n = 0
    ct = 0
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 0
    ct = 1
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 1
    ct = 0
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])

    n = 2
    ct = 2
    assert not data_left_to_try_with(hours_per_file, n, ct, file_times, end_times[n])
