"""Validated scheduled-observation models."""

from datetime import datetime
from typing import Tuple

from pydantic import BaseModel, ConfigDict, Field

from samscheduler.exposure import Exposure
from samscheduler.target import Target


class ScheduleEntry(BaseModel):
    """A target and its complete exposure sequence on the observing timeline."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    target: Target
    exposures: Tuple[Exposure, ...] = Field(min_length=1)
    start_utc: datetime
    end_utc: datetime
    airmass: float = Field(ge=1)

    @property
    def duration_seconds(self) -> float:
        """Return the scheduled exposure duration."""
        return (self.end_utc - self.start_utc).total_seconds()
