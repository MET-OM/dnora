from dnora.utils.time import create_time_stamps
import pandas as pd


def test_wam800():
    """Files 00 and 12"""
    stride = 12
    hours_per_file = 73
    for lead_time in [0, 12, 24]:
        for hh in [00, 12]:
            start_time = pd.to_datetime(f"2024-03-21 {hh:02.0f}:00:00")
            end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
            last_file = start_time - pd.Timedelta(hours=lead_time)

            start_times, end_times, file_times = create_time_stamps(
                start_time,
                end_time,
                stride=stride,
                hours_per_file=hours_per_file,
                last_file=last_file,
                lead_time=lead_time,
            )

            assert len(start_times) == 1
            assert len(end_times) == 1
            assert len(file_times) == 1
            assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == start_time.strftime(
                "%Y-%m-%d %H:%M:00"
            )
            assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == (
                end_time - pd.Timedelta(hours=lead_time)
            ).strftime("%Y-%m-%d %H:%M:00")
            assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == (
                start_time - pd.Timedelta(hours=lead_time)
            ).strftime("%Y-%m-%d %H:%M:00")


def test_wam800_earlier_start_date_00():
    """Files 00 and 12"""
    stride = 12
    hours_per_file = 73
    start_time = pd.to_datetime(f"2024-03-21 00:00:00")
    end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
    last_file = pd.to_datetime(f"2024-03-21 12:00:00")

    start_times, end_times, file_times = create_time_stamps(
        start_time,
        end_time,
        stride=stride,
        hours_per_file=hours_per_file,
        last_file=last_file,
    )

    assert len(start_times) == 2
    assert len(end_times) == 2
    assert len(file_times) == 2
    assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 00:00:00"
    assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == end_time.strftime(
        "%Y-%m-%d %H:%M:00"
    )
    assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 11:00:00"

    assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 00:00:00"
    assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 12:00:00"


def test_wam800_earlier_start_date_18():
    """Files 00 and 12"""
    stride = 12
    hours_per_file = 73
    start_time = pd.to_datetime(f"2024-03-20 18:00:00")
    end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
    last_file = pd.to_datetime(f"2024-03-21 00:00:00")

    start_times, end_times, file_times = create_time_stamps(
        start_time,
        end_time,
        stride=stride,
        hours_per_file=hours_per_file,
        last_file=last_file,
    )

    assert len(start_times) == 2
    assert len(end_times) == 2
    assert len(file_times) == 2
    assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-20 18:00:00"
    assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == end_time.strftime(
        "%Y-%m-%d %H:%M:00"
    )
    assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-20 23:00:00"

    assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-20 12:00:00"
    assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 00:00:00"


def test_wam3():
    """Files 06 and 18"""
    stride = 12
    hours_per_file = 121
    offset = 6
    for lead_time in [0, 12, 24]:
        for hh in [6, 18]:
            start_time = pd.to_datetime(f"2024-03-21 {hh:02.0f}:00:00")
            end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
            last_file = start_time - pd.Timedelta(hours=lead_time)

            start_times, end_times, file_times = create_time_stamps(
                start_time,
                end_time,
                stride=stride,
                hours_per_file=hours_per_file,
                last_file=last_file,
                lead_time=lead_time,
                offset=offset,
            )

            assert len(start_times) == 1
            assert len(end_times) == 1
            assert len(file_times) == 1
            assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == start_time.strftime(
                "%Y-%m-%d %H:%M:00"
            )
            assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == (
                end_time - pd.Timedelta(hours=lead_time)
            ).strftime("%Y-%m-%d %H:%M:00")
            assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == (
                start_time - pd.Timedelta(hours=lead_time)
            ).strftime("%Y-%m-%d %H:%M:00")


def test_wam3_earlier_start_date_00():
    """Files 06 and 18"""
    stride = 12
    hours_per_file = 121
    offset = 6
    start_time = pd.to_datetime(f"2024-03-21 00:00:00")
    end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
    last_file = pd.to_datetime(f"2024-03-21 06:00:00")

    start_times, end_times, file_times = create_time_stamps(
        start_time,
        end_time,
        stride=stride,
        hours_per_file=hours_per_file,
        last_file=last_file,
        offset=offset,
    )

    assert len(start_times) == 2
    assert len(end_times) == 2
    assert len(file_times) == 2
    assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 00:00:00"
    assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == end_time.strftime(
        "%Y-%m-%d %H:%M:00"
    )
    assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 05:00:00"

    assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-20 18:00:00"
    assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 06:00:00"


def test_wam3_earlier_start_date_12():
    """Files 06 and 18"""
    stride = 12
    hours_per_file = 121
    offset = 6
    start_time = pd.to_datetime(f"2024-03-21 12:00:00")
    end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
    last_file = pd.to_datetime(f"2024-03-21 18:00:00")

    start_times, end_times, file_times = create_time_stamps(
        start_time,
        end_time,
        stride=stride,
        hours_per_file=hours_per_file,
        last_file=last_file,
        offset=offset,
    )

    assert len(start_times) == 2
    assert len(end_times) == 2
    assert len(file_times) == 2
    assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 12:00:00"
    assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == end_time.strftime(
        "%Y-%m-%d %H:%M:00"
    )
    assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 17:00:00"

    assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 06:00:00"
    assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 18:00:00"


def test_wam3_loop_start():
    """Files 06 and 18"""
    stride = 12
    hours_per_file = 121
    offset = 6

    file0 = (
        ["2024-03-20 18:00:00"] * 6
        + ["2024-03-21 06:00:00"] * 12
        + ["2024-03-21 18:00:00"] * 6
    )

    lens0 = [3] * 6 + [2] * 12 + [1] * 6

    t1 = ["2024-03-21 05:00:00"] * 6 + ["2024-03-21 17:00:00"] * 12

    for h0 in range(0, 24):
        start_time = pd.to_datetime(f"2024-03-21 00:00:00") + pd.Timedelta(hours=h0)
        end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
        last_file = pd.to_datetime(f"2024-03-21 18:00:00")

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            offset=offset,
        )

        assert len(start_times) == lens0[h0]
        assert len(end_times) == lens0[h0]
        assert len(file_times) == lens0[h0]
        assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == start_time.strftime(
            "%Y-%m-%d %H:%M:00"
        )

        if h0 < 18:
            assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == t1[h0]
        assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == min(
            end_time, last_file + pd.Timedelta(hours=hours_per_file - 1)
        ).strftime("%Y-%m-%d %H:%M:00")
        assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == file0[h0]
        assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-21 18:00:00"


def test_wam800_loop_start():
    """Files 00 and 12"""
    stride = 12
    hours_per_file = 73

    file0 = (
        ["2024-03-21 00:00:00"] * 6
        + ["2024-03-21 12:00:00"] * 12
        + ["2024-03-22 00:00:00"] * 6
    )

    lens0 = [3] * 6 + [2] * 12 + [1] * 6

    t1 = ["2024-03-21 11:00:00"] * 6 + ["2024-03-21 23:00:00"] * 12

    for h0 in range(0, 24):
        start_time = pd.to_datetime(f"2024-03-21 06:00:00") + pd.Timedelta(hours=h0)
        end_time = start_time + pd.Timedelta(hours=hours_per_file - 1)
        last_file = pd.to_datetime(f"2024-03-22 00:00:00")

        start_times, end_times, file_times = create_time_stamps(
            start_time,
            end_time,
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
        )

        assert len(start_times) == lens0[h0]
        assert len(end_times) == lens0[h0]
        assert len(file_times) == lens0[h0]
        assert start_times[0].strftime("%Y-%m-%d %H:%M:00") == start_time.strftime(
            "%Y-%m-%d %H:%M:00"
        )
        if h0 < 18:
            assert end_times[0].strftime("%Y-%m-%d %H:%M:00") == t1[h0]

        assert end_times[-1].strftime("%Y-%m-%d %H:%M:00") == min(
            end_time, last_file + pd.Timedelta(hours=hours_per_file - 1)
        ).strftime("%Y-%m-%d %H:%M:00")
        assert file_times[0].strftime("%Y-%m-%d %H:%M:00") == file0[h0]
        assert file_times[-1].strftime("%Y-%m-%d %H:%M:00") == "2024-03-22 00:00:00"
