"""
Data source (CSV):
    https://eclipse.gsfc.nasa.gov/eclipse_besselian_from_mysqldump2.csv

Attribution requested by NASA:
    "Eclipse Predictions by Fred Espenak, NASA's GSFC"
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Literal, Optional, Tuple

from astropy.time import Time
from pydantic import BaseModel, ConfigDict, Field

EclipseTypeMain = Literal["P", "A", "T", "H"]
EclipseKind = Literal["none", "partial", "annular", "total"]

RAD = math.pi / 180.0

# Values match the Meeus-style examples commonly used with NASA polynomial elements.
EARTH_EQ_RADIUS_M = 6_378_140.0
EARTH_1ME2_SQRT = 0.99664719
DEG_PER_SEC = 0.00417807  # Earth's sidereal rotation rate in deg/s


class Location(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    latitude_deg: float = Field(..., ge=-90.0, le=90.0)
    longitude_deg_east: float = Field(..., ge=-180.0, le=180.0)  # east-positive
    altitude_m: float = Field(...)


class BesselianCoefficients(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    t0_tdt_hours: float  # t0 in TDT decimal hours

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
    model_config = ConfigDict(extra="forbid", frozen=True)

    year: int
    month: int
    day: int

    eclipse_type_raw: str = Field(..., alias="eclipse_type")
    main_type: EclipseTypeMain

    delta_t_seconds: float = Field(..., alias="dt")
    julian_date_ge_tdt: float = Field(..., alias="julian_date")
    td_ge: str  # HH:MM:SS in TDT

    cat_no: int
    bessel: BesselianCoefficients

    # Peak coordinates from NASA (only defined for central eclipses)
    lat_ge: Optional[float] = (
        None  # latitude of greatest eclipse (degrees, north positive)
    )
    lon_ge: Optional[float] = (
        None  # longitude of greatest eclipse (degrees, east positive)
    )

    @staticmethod
    def parse_main_type(raw: str) -> EclipseTypeMain:
        if not raw:
            raise ValueError("Empty eclipse_type")
        c = raw[0].upper()
        if c not in ("P", "A", "T", "H"):
            raise ValueError(f"Unknown eclipse_type {raw!r}")
        return c  # type: ignore[return-value]


class FundamentalArgs(BaseModel):
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

    H: float  # local hour angle (deg)

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
    model_config = ConfigDict(extra="forbid", frozen=True)

    eclipse: EclipseRecord
    kind: EclipseKind

    t_max: float
    jd_ut_max: float

    # NOTE: this is a distance on the fundamental plane in Earth radii (not meters).
    m_er: float
    L1p: float
    L2p: float

    jd_ut_c2: Optional[float] = None  # start of totality (only if kind == "total")
    jd_ut_c3: Optional[float] = None  # end of totality (only if kind == "total")


def parse_hms_to_hours(hms: str) -> float:
    parts = hms.strip().split(":")
    if len(parts) != 3:
        raise ValueError(f"Expected HH:MM:SS, got {hms!r}")
    h = int(parts[0])
    m = int(parts[1])
    s = float(parts[2])
    return h + (m / 60.0) + (s / 3600.0)


def poly3(c0: float, c1: float, c2: float, c3: float, t: float) -> float:
    return ((c3 * t + c2) * t + c1) * t + c0


def poly2(c0: float, c1: float, c2: float, t: float) -> float:
    return (c2 * t + c1) * t + c0


def fundamental_args(e: EclipseRecord, loc: Location, t: float) -> FundamentalArgs:
    b = e.bessel

    X = poly3(b.x0, b.x1, b.x2, b.x3, t)
    Y = poly3(b.y0, b.y1, b.y2, b.y3, t)
    d = poly2(b.d0, b.d1, b.d2, t)
    mu = poly2(b.mu0, b.mu1, b.mu2, t)

    Xp = b.x1 + 2.0 * b.x2 * t + 3.0 * b.x3 * t * t
    Yp = b.y1 + 2.0 * b.y2 * t + 3.0 * b.y3 * t * t
    d_dot = b.d1 + 2.0 * b.d2 * t
    mu_dot = b.mu1 + 2.0 * b.mu2 * t

    L1 = poly2(b.l10, b.l11, b.l12, t)
    L2 = poly2(b.l20, b.l21, b.l22, t)

    # Longitude convention:
    # - loc.longitude_deg_east is east-positive.
    # - The Besselian elements are in TDT/TT; convert to UT by subtracting Earth's rotation during ΔT.
    H = mu + loc.longitude_deg_east - DEG_PER_SEC * e.delta_t_seconds

    phi = loc.latitude_deg
    h_m = loc.altitude_m

    u1 = math.atan(EARTH_1ME2_SQRT * math.tan(phi * RAD)) / RAD
    rho_sin_phi_prime = EARTH_1ME2_SQRT * math.sin(u1 * RAD) + (
        h_m / EARTH_EQ_RADIUS_M
    ) * math.sin(phi * RAD)
    rho_cos_phi_prime = math.cos(u1 * RAD) + (h_m / EARTH_EQ_RADIUS_M) * math.cos(
        phi * RAD
    )

    xi = rho_cos_phi_prime * math.sin(H * RAD)
    eta = rho_sin_phi_prime * math.cos(d * RAD) - rho_cos_phi_prime * math.cos(
        H * RAD
    ) * math.sin(d * RAD)
    zeta = rho_sin_phi_prime * math.sin(d * RAD) + rho_cos_phi_prime * math.cos(
        H * RAD
    ) * math.cos(d * RAD)

    xi_p = RAD * mu_dot * rho_cos_phi_prime * math.cos(H * RAD)
    eta_p = RAD * (mu_dot * xi * math.sin(d * RAD) - zeta * d_dot)

    L1p = L1 - zeta * b.tan_f1
    L2p = L2 - zeta * b.tan_f2

    u = X - xi
    v = Y - eta
    a = Xp - xi_p
    bb = Yp - eta_p
    n = math.hypot(a, bb)

    return FundamentalArgs(
        t=t,
        X=X,
        Y=Y,
        d=d,
        mu=mu,
        Xp=Xp,
        Yp=Yp,
        d_dot=d_dot,
        mu_dot=mu_dot,
        L1=L1,
        L2=L2,
        H=H,
        rho_sin_phi_prime=rho_sin_phi_prime,
        rho_cos_phi_prime=rho_cos_phi_prime,
        xi=xi,
        eta=eta,
        zeta=zeta,
        xi_p=xi_p,
        eta_p=eta_p,
        L1p=L1p,
        L2p=L2p,
        u=u,
        v=v,
        a=a,
        b=bb,
        n=n,
    )


def solve_t_max(
    e: EclipseRecord, loc: Location, *, max_iter: int = 10, tol: float = 1e-8
) -> Tuple[float, FundamentalArgs]:
    t = 0.0
    for _ in range(max_iter):
        fa = fundamental_args(e, loc, t)
        tau = -((fa.u * fa.a) + (fa.v * fa.b)) / (fa.n * fa.n)
        t = t + tau
        if abs(tau) < tol:
            break
    fa = fundamental_args(e, loc, t)
    return t, fa


def classify_at_max(fa: FundamentalArgs) -> EclipseKind:
    # NOTE: do NOT reject on zeta here.
    # A sunrise/sunset total eclipse can have zeta cross 0 between C2 and C3.
    # We instead apply a dedicated horizon check after computing C2/C3 for total eclipses.
    m = math.hypot(fa.u, fa.v)
    if m > fa.L1p:
        return "none"
    if m < abs(fa.L2p):
        if fa.L2p < 0:
            return "total"
        return "annular"
    return "partial"


def jd_tdt_from_t(e: EclipseRecord, t: float) -> float:
    td_ge_hours = parse_hms_to_hours(e.td_ge)
    # Handle day wrap: t0 is the nearest integer hour to td_ge, so when td_ge
    # is close to midnight (e.g. 23:45), t0=0 refers to the NEXT day.
    delta_hours = e.bessel.t0_tdt_hours - td_ge_hours
    if delta_hours < -12:
        delta_hours += 24
    elif delta_hours > 12:
        delta_hours -= 24
    jd_t0 = e.julian_date_ge_tdt + delta_hours / 24.0
    return jd_t0 + t / 24.0


def jd_ut_from_t(e: EclipseRecord, t: float) -> float:
    return jd_tdt_from_t(e, t) - (e.delta_t_seconds / 86_400.0)


def refine_contact(
    e: EclipseRecord,
    loc: Location,
    t_guess: float,
    *,
    radius: Literal["L1p", "L2p"],
    which: Literal["minus", "plus"],
    max_iter: int = 10,
) -> float:
    t = t_guess
    for _ in range(max_iter):
        fa = fundamental_args(e, loc, t)
        L = fa.L1p if radius == "L1p" else abs(fa.L2p)

        S = (fa.a * fa.v - fa.u * fa.b) / (fa.n * L)
        root = math.sqrt(1.0 - S * S)

        base = -((fa.u * fa.a) + (fa.v * fa.b)) / (fa.n * fa.n)
        corr = (L / fa.n) * root

        tau = base - corr if which == "minus" else base + corr
        t = t + tau
    return t


def totality_visible(e: EclipseRecord, loc: Location, t_c2: float, t_c3: float) -> bool:
    # Reject the "night-side" solution:
    # If totality is entirely below the horizon (zeta <= 0 for the full C2..C3 interval),
    # then this location cannot observe totality.
    z2 = fundamental_args(e, loc, t_c2).zeta
    z3 = fundamental_args(e, loc, t_c3).zeta
    return max(z2, z3) > 0.0


def compute_local_eclipse(e: EclipseRecord, loc: Location) -> LocalEclipseResult:
    t_max, fa = solve_t_max(e, loc)
    jd_ut_max = jd_ut_from_t(e, t_max)
    m = math.hypot(fa.u, fa.v)

    def make_result(
        kind: EclipseKind,
        jd_ut_c2: Optional[float] = None,
        jd_ut_c3: Optional[float] = None,
    ) -> LocalEclipseResult:
        return LocalEclipseResult(
            eclipse=e,
            kind=kind,
            t_max=t_max,
            jd_ut_max=jd_ut_max,
            m_er=m,
            L1p=fa.L1p,
            L2p=fa.L2p,
            jd_ut_c2=jd_ut_c2,
            jd_ut_c3=jd_ut_c3,
        )

    # If we drift outside the validity interval of the polynomial fit, treat as not-visible.
    if not (e.bessel.tmin <= t_max <= e.bessel.tmax):
        return make_result("none")

    kind = classify_at_max(fa)

    if kind != "total":
        return make_result(kind)

    L = abs(fa.L2p)
    S = (fa.a * fa.v - fa.u * fa.b) / (fa.n * L)
    tau = (L / fa.n) * math.sqrt(1.0 - S * S)

    t_c2 = refine_contact(e, loc, t_max - tau, radius="L2p", which="minus")
    t_c3 = refine_contact(e, loc, t_max + tau, radius="L2p", which="plus")

    if not totality_visible(e, loc, t_c2, t_c3):
        return make_result("none")

    return make_result("total", jd_ut_from_t(e, t_c2), jd_ut_from_t(e, t_c3))


def load_catalog(csv_path: str) -> list[EclipseRecord]:
    p = Path(csv_path)
    if not p.exists():
        raise FileNotFoundError(p)

    records: list[EclipseRecord] = []
    with p.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)

        required = {
            "year",
            "month",
            "day",
            "eclipse_type",
            "dt",
            "julian_date",
            "td_ge",
            "cat_no",
            "t0",
            "x0",
            "x1",
            "x2",
            "x3",
            "y0",
            "y1",
            "y2",
            "y3",
            "d0",
            "d1",
            "d2",
            "mu0",
            "mu1",
            "mu2",
            "l10",
            "l11",
            "l12",
            "l20",
            "l21",
            "l22",
            "tan_f1",
            "tan_f2",
            "tmin",
            "tmax",
        }
        header = set(reader.fieldnames or [])
        missing = required - header
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")

        for row in reader:
            raw_type = row["eclipse_type"]
            f = {
                k: float(row[k])
                for k in row
                if k
                not in ("eclipse_type", "td_ge", "lat_ge", "lng_ge", "central_duration")
            }
            b = BesselianCoefficients(
                t0_tdt_hours=f["t0"],
                x0=f["x0"],
                x1=f["x1"],
                x2=f["x2"],
                x3=f["x3"],
                y0=f["y0"],
                y1=f["y1"],
                y2=f["y2"],
                y3=f["y3"],
                d0=f["d0"],
                d1=f["d1"],
                d2=f["d2"],
                mu0=f["mu0"],
                mu1=f["mu1"],
                mu2=f["mu2"],
                l10=f["l10"],
                l11=f["l11"],
                l12=f["l12"],
                l20=f["l20"],
                l21=f["l21"],
                l22=f["l22"],
                tan_f1=f["tan_f1"],
                tan_f2=f["tan_f2"],
                tmin=f["tmin"],
                tmax=f["tmax"],
            )

            # Parse peak coordinates (only defined for central eclipses)
            lat_ge = None
            lon_ge = None
            if row.get("lat_dd_ge") and row.get("lng_dd_ge"):
                try:
                    lat_ge = float(row["lat_dd_ge"])
                    lon_ge = float(row["lng_dd_ge"])
                except (ValueError, TypeError):
                    pass

            records.append(
                EclipseRecord(
                    year=int(f["year"]),
                    month=int(f["month"]),
                    day=int(f["day"]),
                    eclipse_type=raw_type,
                    main_type=EclipseRecord.parse_main_type(raw_type),
                    dt=f["dt"],
                    julian_date=f["julian_date"],
                    td_ge=row["td_ge"],
                    cat_no=int(f["cat_no"]),
                    bessel=b,
                    lat_ge=lat_ge,
                    lon_ge=lon_ge,
                )
            )

    return records


def find_last_next_total(
    catalog: list[EclipseRecord],
    loc: Location,
    ref_jd_ut: float,
) -> Tuple[Optional[LocalEclipseResult], Optional[LocalEclipseResult]]:
    prev_res: Optional[LocalEclipseResult] = None
    next_res: Optional[LocalEclipseResult] = None

    for e in catalog:
        if e.main_type not in ("T", "H"):
            continue

        res = compute_local_eclipse(e, loc)
        if res.kind != "total":
            continue

        jd = res.jd_ut_max
        if jd < ref_jd_ut:
            if prev_res is None or jd > prev_res.jd_ut_max:
                prev_res = res
        else:
            if next_res is None or jd < next_res.jd_ut_max:
                next_res = res

    return prev_res, next_res


def print_result(label: str, r: Optional[LocalEclipseResult]) -> None:
    if r is None:
        print(f"{label}: none found in catalog range")
        return

    e = r.eclipse
    print(
        f"{label}: cat_no={e.cat_no}  date={e.year:04d}-{e.month:02d}-{e.day:02d}  type={e.eclipse_type_raw}"
    )
    print(
        f"  max: {Time(r.jd_ut_max, format='jd', scale='utc').iso} UTC  (kind={r.kind})"
    )

    if r.jd_ut_c2 is not None and r.jd_ut_c3 is not None:
        dur_s = (r.jd_ut_c3 - r.jd_ut_c2) * 86400.0
        print(f"  C2 : {Time(r.jd_ut_c2, format='jd', scale='utc').iso} UTC")
        print(f"  C3 : {Time(r.jd_ut_c3, format='jd', scale='utc').iso} UTC")
        print(f"  duration: {dur_s:.1f} s")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv", required=True, help="Path to eclipse_besselian_from_mysqldump2.csv"
    )
    ap.add_argument(
        "--lat", required=True, type=float, help="Latitude in degrees (north positive)"
    )
    ap.add_argument(
        "--lon",
        required=True,
        type=float,
        help="Longitude in degrees east (east positive; west is negative)",
    )
    ap.add_argument(
        "--alt", type=float, default=0.0, help="Altitude in meters (default: 0)"
    )
    ap.add_argument(
        "--ref-utc",
        default=None,
        help="Reference time in ISO8601 (e.g. 2026-01-03T00:00:00Z). Default: now.",
    )
    args = ap.parse_args()

    if args.ref_utc is None:
        ref_time = Time.now()
        print(f"ref-utc not provided; using now: {ref_time.iso}")
    else:
        ref_time = Time(args.ref_utc, scale="utc")

    loc = Location(
        latitude_deg=args.lat, longitude_deg_east=args.lon, altitude_m=args.alt
    )
    catalog = load_catalog(args.csv)
    ref_jd = ref_time.jd

    prev_r, next_r = find_last_next_total(catalog, loc, ref_jd)

    print(
        f"Location: lat={loc.latitude_deg:.6f}, lon_east={loc.longitude_deg_east:.6f}, alt_m={loc.altitude_m:.1f}"
    )
    print(f"Reference: {ref_time.iso} UTC (JD_UT={ref_jd:.6f})")
    print_result("Previous total eclipse", prev_r)
    print_result("Next total eclipse", next_r)


if __name__ == "__main__":
    main()
