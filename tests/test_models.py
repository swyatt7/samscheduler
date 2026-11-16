"""Tests for the validated public data models."""

import pytest
from pydantic import ValidationError

from samscheduler import DEFAULT_TARGETS, Exposure, Filter, ScheduleEntry, Target


def test_default_catalog_has_visible_coordinates_and_ubvri_exposures():
    assert len(DEFAULT_TARGETS) == 5
    expected_filters = set(Filter)
    for target in DEFAULT_TARGETS:
        assert 0 <= target.ra_deg < 360
        assert -90 <= target.dec_deg <= 90
        assert {exposure.filter for exposure in target.exposures} == expected_filters
        assert all(200 <= exposure.seconds <= 500 for exposure in target.exposures)


@pytest.mark.parametrize(
    "values",
    [
        {"filter": "X", "seconds": 300},
        {"filter": "V", "seconds": 199},
        {"filter": "R", "seconds": 501},
        {"filter": "B", "seconds": 300, "count": 0},
    ],
)
def test_exposure_rejects_invalid_filter_or_range(values):
    with pytest.raises(ValidationError):
        Exposure(**values)


@pytest.mark.parametrize(
    "field,value",
    [
        ("ra_deg", 360),
        ("ra_deg", -0.1),
        ("dec_deg", 90.1),
        ("priority", -1),
        ("max_airmass", 0.9),
    ],
)
def test_target_rejects_invalid_coordinates_and_constraints(field, value):
    values = {
        "name": "Test",
        "ra_deg": 20,
        "dec_deg": 10,
        "exposures": [{"filter": "V", "seconds": 200}],
        field: value,
    }
    with pytest.raises(ValidationError):
        Target(**values)


def test_target_coordinate_and_schedule_duration():
    target = Target(
        name="Test",
        ra_deg=20,
        dec_deg=10,
        exposures=(Exposure(filter=Filter.V, seconds=200),),
    )
    assert target.skycoord.ra.deg == pytest.approx(20)

    entry = ScheduleEntry(
        target=target,
        exposures=(Exposure(filter=Filter.V, seconds=200),),
        start_utc="2026-01-01T00:00:00+00:00",
        end_utc="2026-01-01T00:03:20+00:00",
        airmass=1.2,
    )
    assert entry.duration_seconds == pytest.approx(200)
