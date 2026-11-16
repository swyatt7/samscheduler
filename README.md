# samscheduler

A small, importable astronomy observing scheduler. The library uses Astropy for
coordinate transforms and Pydantic models for target, exposure, and schedule
data. It has no database, telescope-queue, or network-service integration.

## Install

```shell
python -m pip install -r requirements.txt
python -m pip install -e .
```

For development and tests:

```shell
python -m pip install -r requirements-dev.txt
pytest
```

## Use

The package ships with a small example catalog. Coordinates are ICRS degrees,
and each example target has one exposure in each of the UBVRI filters between
200 and 500 seconds.

```python
from datetime import date

from samscheduler import DEFAULT_TARGETS, Scheduler

scheduler = Scheduler(DEFAULT_TARGETS, observing_date=date(2026, 11, 16))
schedule = scheduler.run()
figure = scheduler.plot()
figure.savefig("schedule.png")

for target_block in schedule:
    print(
        target_block.target.name,
        [(exposure.filter.value, exposure.seconds, exposure.count)
         for exposure in target_block.exposures],
        target_block.start_utc,
        target_block.airmass,
    )
```

`observing_date` is the local calendar date on which the night begins.
`local_utc_offset_hours` defaults to `-7` (Arizona standard time); the observing
window starts at 18:00 local time and runs for 12 hours. The scheduler greedily
selects the highest-priority observable target, then the lowest airmass, and
schedules that target's complete exposure list as one block. A smaller numeric
priority is higher priority. `run()` replaces and returns one schedule entry
per scheduled target; `plot()` returns a Matplotlib figure with target airmass
tracks and scheduled target blocks.
