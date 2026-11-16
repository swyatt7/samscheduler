"""Example observatory and static target catalog."""

from astropy import units as u
from astropy.coordinates import EarthLocation

from samscheduler.exposure import Exposure, Filter
from samscheduler.target import Target

OBSERVATORY = EarthLocation(
    lat=32.4167 * u.deg,
    lon=-110.7344 * u.deg,
    height=2510 * u.m,
)

DEFAULT_TARGETS = (
    Target(
        name="Andromeda Galaxy",
        ra_deg=10.6847,
        dec_deg=41.2692,
        priority=5,
        exposures=(
            Exposure(filter=Filter.U, seconds=400),
            Exposure(filter=Filter.B, seconds=350),
            Exposure(filter=Filter.V, seconds=300),
            Exposure(filter=Filter.R, seconds=300),
            Exposure(filter=Filter.I, seconds=350),
        ),
    ),
    Target(
        name="Orion Nebula",
        ra_deg=83.8221,
        dec_deg=-5.3911,
        priority=1,
        exposures=(
            Exposure(filter=Filter.U, seconds=300),
            Exposure(filter=Filter.B, seconds=250),
            Exposure(filter=Filter.V, seconds=200),
            Exposure(filter=Filter.R, seconds=250),
            Exposure(filter=Filter.I, seconds=300),
        ),
    ),
    Target(
        name="Whirlpool Galaxy",
        ra_deg=202.4696,
        dec_deg=47.1952,
        priority=2,
        exposures=(
            Exposure(filter=Filter.U, seconds=450),
            Exposure(filter=Filter.B, seconds=400),
            Exposure(filter=Filter.V, seconds=350),
            Exposure(filter=Filter.R, seconds=350),
            Exposure(filter=Filter.I, seconds=400),
        ),
    ),
    Target(
        name="North America Nebula",
        ra_deg=314.7000,
        dec_deg=44.5000,
        priority=2,
        exposures=(
            Exposure(filter=Filter.U, seconds=350),
            Exposure(filter=Filter.B, seconds=300),
            Exposure(filter=Filter.V, seconds=250),
            Exposure(filter=Filter.R, seconds=300),
            Exposure(filter=Filter.I, seconds=350),
        ),
    ),
    Target(
        name="Bode's Galaxy",
        ra_deg=148.8882,
        dec_deg=69.0653,
        priority=1,
        exposures=(
            Exposure(filter=Filter.U, seconds=500),
            Exposure(filter=Filter.B, seconds=450),
            Exposure(filter=Filter.V, seconds=400),
            Exposure(filter=Filter.R, seconds=400),
            Exposure(filter=Filter.I, seconds=450),
        ),
    ),
)
