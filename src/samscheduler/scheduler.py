"""Greedy, visibility-constrained observing scheduler."""

from datetime import date, datetime, time, timedelta, timezone
from typing import Iterable, List, Optional, Tuple

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, get_body, get_sun
from astropy.time import Time
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from samscheduler.catalog import DEFAULT_TARGETS, OBSERVATORY
from samscheduler.schedule import ScheduleEntry
from samscheduler.target import Target


class Scheduler:
    """Build and visualize a one-night observing schedule.

    ``observing_date`` is the local calendar date on which the night begins.
    The fixed UTC offset defaults to Arizona standard time (UTC-7).
    """

    def __init__(
        self,
        targets: Iterable[Target] = DEFAULT_TARGETS,
        observing_date: Optional[date] = None,
        location=OBSERVATORY,
        local_utc_offset_hours: float = -7,
        start_hour: int = 18,
        duration_hours: float = 12,
        time_step_minutes: int = 5,
        readout_seconds: int = 10,
        slew_seconds: int = 60,
        sun_altitude_limit_deg: float = -12,
    ) -> None:
        if not -12 <= local_utc_offset_hours <= 14:
            raise ValueError("local_utc_offset_hours must be between -12 and 14")
        if not 0 <= start_hour <= 23:
            raise ValueError("start_hour must be between 0 and 23")
        if duration_hours <= 0:
            raise ValueError("duration_hours must be positive")
        if time_step_minutes <= 0:
            raise ValueError("time_step_minutes must be positive")
        if readout_seconds < 0 or slew_seconds < 0:
            raise ValueError("readout_seconds and slew_seconds cannot be negative")

        self.targets = tuple(targets)
        self.observing_date = observing_date or date.today()
        self.location = location
        self.local_utc_offset_hours = local_utc_offset_hours
        self.start_hour = start_hour
        self.duration_hours = duration_hours
        self.time_step_minutes = time_step_minutes
        self.readout_seconds = readout_seconds
        self.slew_seconds = slew_seconds
        self.sun_altitude_limit_deg = sun_altitude_limit_deg
        self.schedule: List[ScheduleEntry] = []

    def _night_bounds(self) -> Tuple[Time, Time]:
        local_start = datetime.combine(self.observing_date, time(self.start_hour))
        utc_start = local_start - timedelta(hours=self.local_utc_offset_hours)
        start = Time(utc_start, scale="utc")
        return start, start + self.duration_hours * u.hour

    def _frame(self, at: Time) -> AltAz:
        return AltAz(obstime=at, location=self.location, pressure=0 * u.hPa)

    def _airmass(self, target: Target, at: Time) -> Optional[float]:
        altitude = target.skycoord.transform_to(self._frame(at))
        if altitude.alt.deg <= 0:
            return None
        value = float(altitude.secz.value)
        if not np.isfinite(value) or value > target.max_airmass:
            return None
        return value

    def _is_dark(self, at: Time) -> bool:
        sun_altitude = get_sun(at).transform_to(self._frame(at)).alt.deg
        return bool(sun_altitude < self.sun_altitude_limit_deg)

    def _moon_is_clear(self, target: Target, at: Time) -> bool:
        if target.min_moon_separation_deg == 0:
            return True
        moon = get_body("moon", at, location=self.location)
        return bool(target.skycoord.separation(moon).deg >= target.min_moon_separation_deg)

    def _target_block_end(self, target: Target, start: Time, night_end: Time) -> Optional[Time]:
        """Return the block end if all of a target's exposures remain observable."""
        exposure_sequence = [
            exposure
            for exposure in target.exposures
            for _ in range(exposure.count)
        ]
        current = start
        block_end = start

        for index, exposure in enumerate(exposure_sequence):
            exposure_end = current + exposure.seconds * u.second
            if exposure_end > night_end:
                return None
            for at in (current, exposure_end):
                if (
                    self._airmass(target, at) is None
                    or not self._is_dark(at)
                    or not self._moon_is_clear(target, at)
                ):
                    return None
            block_end = exposure_end
            if index < len(exposure_sequence) - 1:
                current = exposure_end + self.readout_seconds * u.second

        return block_end

    def run(self) -> List[ScheduleEntry]:
        """Populate the schedule with one entry per observable target."""
        start, end = self._night_bounds()
        pending = list(self.targets)
        self.schedule = []
        current = start
        step = self.time_step_minutes * u.minute

        while pending and current < end:
            if not self._is_dark(current):
                current += step
                continue

            visible = []
            for target in pending:
                airmass = self._airmass(target, current)
                block_end = self._target_block_end(target, current, end)
                if airmass is not None and block_end is not None:
                    visible.append((target, airmass, block_end))

            if not visible:
                current += step
                continue

            target, airmass, block_end = min(
                visible,
                key=lambda item: (item[0].priority, item[1], item[0].name),
            )

            self.schedule.append(
                ScheduleEntry(
                    target=target,
                    exposures=target.exposures,
                    start_utc=current.to_datetime(timezone=timezone.utc),
                    end_utc=block_end.to_datetime(timezone=timezone.utc),
                    airmass=airmass,
                )
            )
            pending.remove(target)
            current = block_end + self.slew_seconds * u.second

        return list(self.schedule)

    def plot(self) -> Figure:
        """Return a Matplotlib figure of target airmass tracks and exposures."""
        start, end = self._night_bounds()
        sample_count = max(2, int(self.duration_hours * 60 / self.time_step_minutes) + 1)
        sample_times = start + np.linspace(0, self.duration_hours, sample_count) * u.hour
        frame = self._frame(sample_times)

        figure = Figure(figsize=(10, 5))
        FigureCanvasAgg(figure)
        axis = figure.add_subplot(1, 1, 1)
        hours = np.linspace(0, self.duration_hours, sample_count)

        for target in self.targets:
            altaz = target.skycoord.transform_to(frame)
            airmass = np.asarray(altaz.secz.value, dtype=float)
            airmass[(altaz.alt.deg <= 0) | (airmass > target.max_airmass)] = np.nan
            axis.plot(hours, airmass, label=target.name)

        for entry in self.schedule:
            offset_start = (
                Time(entry.start_utc) - start
            ).to_value(u.hour)
            offset_end = (
                Time(entry.end_utc) - start
            ).to_value(u.hour)
            axis.axvspan(offset_start, offset_end, color="tab:red", alpha=0.15)
            filters = "/".join(
                "{}x{}".format(exposure.filter.value, exposure.count)
                if exposure.count > 1
                else exposure.filter.value
                for exposure in entry.exposures
            )
            axis.annotate(
                filters,
                ((offset_start + offset_end) / 2, entry.airmass),
                fontsize=7,
                ha="center",
            )

        axis.set_xlim(0, (end - start).to_value(u.hour))
        axis.set_ylim(3.5, 1)
        axis.set_xlabel("Hours from observing-window start (local)")
        axis.set_ylabel("Airmass")
        axis.set_title("Target airmass and scheduled exposures")
        axis.grid(True, alpha=0.3)
        if self.targets:
            axis.legend(loc="best")
        figure.tight_layout()
        return figure
