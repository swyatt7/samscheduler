"""Tests for scheduling and plotting behavior."""

from datetime import date

import pytest

from samscheduler import DEFAULT_TARGETS, Exposure, Filter, Scheduler, Target


def make_target(name="Test", priority=1, count=1):
    return Target(
        name=name,
        ra_deg=255,
        dec_deg=32,
        priority=priority,
        min_moon_separation_deg=0,
        exposures=(Exposure(filter=Filter.V, seconds=200, count=count),),
    )


def make_scheduler(targets, **kwargs):
    options = {
        "observing_date": date(2025, 6, 21),
        "local_utc_offset_hours": -7,
        "time_step_minutes": 5,
    }
    options.update(kwargs)
    return Scheduler(targets, **options)


def test_run_populates_schedule_with_visible_targets_and_repeats():
    target = make_target(count=2)
    scheduler = make_scheduler([target])

    schedule = scheduler.run()

    assert len(schedule) == 1
    assert schedule == scheduler.schedule
    assert all(entry.target.name == "Test" for entry in schedule)
    assert all(entry.airmass <= 2.5 for entry in schedule)
    assert schedule[0].exposures == target.exposures
    assert schedule[0].duration_seconds == pytest.approx(410)


def test_target_with_multiple_filters_is_scheduled_as_one_block():
    target = Target(
        name="Multi-filter",
        ra_deg=255,
        dec_deg=32,
        min_moon_separation_deg=0,
        exposures=(
            Exposure(filter=Filter.U, seconds=200),
            Exposure(filter=Filter.R, seconds=300),
        ),
    )

    schedule = make_scheduler([target]).run()

    assert len(schedule) == 1
    assert schedule[0].target == target
    assert schedule[0].exposures == target.exposures
    assert schedule[0].duration_seconds == pytest.approx(510)


def test_run_respects_priority_and_is_repeatable():
    high_priority = make_target(name="Higher", priority=0)
    lower_priority = make_target(name="Lower", priority=2)
    scheduler = make_scheduler([lower_priority, high_priority])

    first = scheduler.run()
    second = scheduler.run()

    assert [entry.target.name for entry in first] == ["Higher", "Lower"]
    assert [entry.target.name for entry in second] == ["Higher", "Lower"]


def test_default_catalog_produces_a_schedule_on_a_sample_night():
    schedule = Scheduler(
        DEFAULT_TARGETS,
        observing_date=date(2025, 6, 21),
    ).run()

    assert schedule
    assert all(entry.airmass <= entry.target.max_airmass for entry in schedule)
    assert len({entry.target.name for entry in schedule}) == len(schedule)


def test_run_returns_empty_when_targets_are_not_visible():
    below_horizon = Target(
        name="South",
        ra_deg=255,
        dec_deg=-80,
        min_moon_separation_deg=0,
        exposures=(Exposure(filter=Filter.R, seconds=200),),
    )
    assert make_scheduler([below_horizon]).run() == []


def test_run_skips_daylight_and_rejects_invalid_configuration():
    scheduler = make_scheduler(
        [make_target()],
        sun_altitude_limit_deg=-90,
    )
    assert scheduler.run() == []

    with pytest.raises(ValueError, match="duration_hours"):
        make_scheduler([make_target()], duration_hours=0)
    with pytest.raises(ValueError, match="time_step_minutes"):
        make_scheduler([make_target()], time_step_minutes=0)
    with pytest.raises(ValueError, match="local_utc_offset_hours"):
        make_scheduler([make_target()], local_utc_offset_hours=15)
    with pytest.raises(ValueError, match="start_hour"):
        make_scheduler([make_target()], start_hour=24)
    with pytest.raises(ValueError, match="cannot be negative"):
        make_scheduler([make_target()], readout_seconds=-1)


def test_run_does_not_schedule_an_exposure_that_will_not_fit():
    scheduler = make_scheduler(
        [make_target()],
        start_hour=23,
        duration_hours=0.05,
    )
    assert scheduler.run() == []


def test_plot_returns_airmass_figure_with_and_without_schedule():
    scheduler = make_scheduler([make_target()])

    figure = scheduler.plot()
    assert len(figure.axes) == 1
    assert figure.axes[0].get_ylabel() == "Airmass"
    assert figure.axes[0].lines

    scheduler.run()
    scheduled_figure = scheduler.plot()
    assert len(scheduled_figure.axes[0].patches) == 1
