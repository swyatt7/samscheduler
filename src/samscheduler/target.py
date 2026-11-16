"""Validated observing target models."""

from typing import Tuple

from astropy.coordinates import SkyCoord
from astropy import units as u
from pydantic import BaseModel, ConfigDict, Field

from samscheduler.exposure import Exposure


class Target(BaseModel):
    """An ICRS target and its observing constraints."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str = Field(min_length=1)
    ra_deg: float = Field(ge=0, lt=360)
    dec_deg: float = Field(ge=-90, le=90)
    priority: int = Field(default=3, ge=0)
    exposures: Tuple[Exposure, ...] = Field(min_length=1)
    max_airmass: float = Field(default=2.5, ge=1)
    min_moon_separation_deg: float = Field(default=10, ge=0, le=180)

    @property
    def skycoord(self) -> SkyCoord:
        """Return this target's ICRS coordinate."""
        return SkyCoord(ra=self.ra_deg * u.deg, dec=self.dec_deg * u.deg, frame="icrs")
