"""Validated observing exposure models."""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field


class Filter(str, Enum):
    """Standard Johnson-Cousins broadband filters."""

    U = "U"
    B = "B"
    V = "V"
    R = "R"
    I = "I"


class Exposure(BaseModel):
    """One filter's exposure setup for a target."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    filter: Filter
    seconds: int = Field(ge=200, le=500)
    count: int = Field(default=1, ge=1)
