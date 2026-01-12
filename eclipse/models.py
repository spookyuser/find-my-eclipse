"""
Data models for eclipse calculations.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
Data source: https://eclipse.gsfc.nasa.gov/eclipse_besselian_from_mysqldump2.csv
"""

from __future__ import annotations

import math
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

EclipseTypeMain = Literal["P", "A", "T", "H"]
EclipseKind = Literal["none", "partial", "annular", "total"]


class Location(BaseModel):
    """A geographic location on Earth."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    latitude_deg: float = Field(..., ge=-90.0, le=90.0)
    longitude_deg_east: float = Field(..., ge=-180.0, le=180.0)
    altitude_m: float = Field(default=0.0)


class BesselianCoefficients(BaseModel):
    """Besselian elements for an eclipse."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    t0_tdt_hours: float

    x0: float
    x1: float
    x2: float
    x3: float

    y0: float
    y1: float
    y2: float
    y3: float

    d0: float
    d1: float
    d2: float

    mu0: float
    mu1: float
    mu2: float

    l10: float
    l11: float
    l12: float

    l20: float
    l21: float
    l22: float

    tan_f1: float
    tan_f2: float

    tmin: float
    tmax: float


class EclipseRecord(BaseModel):
    """A solar eclipse record from the NASA catalog."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    year: int
    month: int
    day: int

    eclipse_type_raw: str = Field(..., alias="eclipse_type")
    main_type: EclipseTypeMain

    delta_t_seconds: float = Field(..., alias="dt")
    julian_date_ge_tdt: float = Field(..., alias="julian_date")
    td_ge: str

    cat_no: int
    bessel: BesselianCoefficients

    lat_ge: Optional[float] = None
    lon_ge: Optional[float] = None

    @staticmethod
    def parse_main_type(raw: str) -> EclipseTypeMain:
        if not raw:
            raise ValueError("Empty eclipse_type")
        c = raw[0].upper()
        if c not in ("P", "A", "T", "H"):
            raise ValueError(f"Unknown eclipse_type {raw!r}")
        return c  # type: ignore[return-value]

    @property
    def date_iso(self) -> str:
        """Return the eclipse date as YYYY-MM-DD."""
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"


class FundamentalArgs(BaseModel):
    """Intermediate values computed during eclipse calculations."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    t: float

    X: float
    Y: float
    d: float
    mu: float

    Xp: float
    Yp: float
    d_dot: float
    mu_dot: float

    L1: float
    L2: float

    H: float

    rho_sin_phi_prime: float
    rho_cos_phi_prime: float

    xi: float
    eta: float
    zeta: float
    xi_p: float
    eta_p: float

    L1p: float
    L2p: float

    u: float
    v: float
    a: float
    b: float
    n: float


class LocalEclipseResult(BaseModel):
    """Result of computing a local eclipse for a specific location."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    eclipse: EclipseRecord
    kind: EclipseKind

    t_max: float
    jd_ut_max: float

    m_er: float
    L1p: float
    L2p: float

    jd_ut_c2: Optional[float] = None
    jd_ut_c3: Optional[float] = None

    @property
    def duration_seconds(self) -> Optional[float]:
        """Duration of totality in seconds, if applicable."""
        if self.jd_ut_c2 is not None and self.jd_ut_c3 is not None:
            return (self.jd_ut_c3 - self.jd_ut_c2) * 86400.0
        return None
