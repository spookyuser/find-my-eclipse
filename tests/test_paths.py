from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import Final

import astronomy
from hypothesis import assume, given, settings, strategies as st

from eclipse import Location, load_catalog
from eclipse.compute import compute_local_eclipse, solve_t_max
from eclipse.models import EclipseRecord

J2000_JD_UT: Final[float] = 2451545.0


def ae_time_to_jd_ut(t: astronomy.Time) -> float:
    # Astronomy Engine: ut is days since noon UTC on 2000-01-01.
    return J2000_JD_UT + t.ut


def wrap_lon_deg(lon: float) -> float:
    # Normalize to [-180, 180).
    return ((lon + 180.0) % 360.0) - 180.0


def offset_lat_lon_km(
    lat_deg: float, lon_deg: float, *, east_km: float, north_km: float
) -> tuple[float, float]:
    """
    Robust small-offset move on a sphere (great-circle), works near poles.
    east_km/north_km are local tangent-plane offsets.
    """
    r_km = 6371.0
    d_km = math.hypot(east_km, north_km)
    if d_km == 0.0:
        return lat_deg, lon_deg

    # Bearing: 0=north, 90=east.
    bearing = math.atan2(east_km, north_km)
    lat1 = math.radians(lat_deg)
    lon1 = math.radians(lon_deg)
    ang = d_km / r_km

    sin_lat1 = math.sin(lat1)
    cos_lat1 = math.cos(lat1)
    sin_ang = math.sin(ang)
    cos_ang = math.cos(ang)

    lat2 = math.asin(sin_lat1 * cos_ang + cos_lat1 * sin_ang * math.cos(bearing))
    lon2 = lon1 + math.atan2(
        math.sin(bearing) * sin_ang * cos_lat1,
        cos_ang - sin_lat1 * math.sin(lat2),
    )

    return (math.degrees(lat2), wrap_lon_deg(math.degrees(lon2)))


def nasa_jd_ut_ge(e: EclipseRecord) -> float:
    # NASA record gives JD at greatest eclipse in TT/TDT; UT = TT - ΔT.
    return e.julian_date_ge_tdt - (e.delta_t_seconds / 86400.0)


def load_nasa_index() -> tuple[list[float], list[EclipseRecord]]:
    catalog = load_catalog()

    pairs = sorted((nasa_jd_ut_ge(e), e) for e in catalog)
    jd_list = [p[0] for p in pairs]
    rec_list = [p[1] for p in pairs]
    return jd_list, rec_list


NASA_JD_UT, NASA_REC = load_nasa_index()


def nearest_nasa_record(jd_ut: float) -> EclipseRecord:
    i = bisect.bisect_left(NASA_JD_UT, jd_ut)
    if i == 0:
        return NASA_REC[0]
    if i == len(NASA_JD_UT):
        return NASA_REC[-1]
    jd0 = NASA_JD_UT[i - 1]
    jd1 = NASA_JD_UT[i]
    return NASA_REC[i - 1] if abs(jd0 - jd_ut) <= abs(jd1 - jd_ut) else NASA_REC[i]


def map_ae_kind(k: astronomy.EclipseKind) -> str:
    if k == astronomy.EclipseKind.Total:
        return "total"
    if k == astronomy.EclipseKind.Annular:
        return "annular"
    if k == astronomy.EclipseKind.Partial:
        return "partial"
    raise AssertionError(f"Unexpected astronomy.EclipseKind: {k!r}")


@dataclass(frozen=True, slots=True)
class GlobalCase:
    kind: astronomy.EclipseKind
    peak: astronomy.Time
    peak_jd_ut: float
    lat: float
    lon: float


def build_global_cases(start_year: int, end_year: int) -> list[GlobalCase]:
    """
    Build a list of 'known eclipses' from Astronomy Engine:
    global solar eclipses with defined peak shadow center coordinates.

    GlobalSolarEclipseInfo.latitude/longitude are defined only for Total/Annular.
    """
    t = astronomy.Time.Make(start_year, 1, 1, 0, 0, 0.0)
    g = astronomy.SearchGlobalSolarEclipse(t)

    out: list[GlobalCase] = []
    while True:
        y, _, _, _, _, _ = g.peak.Calendar()
        if y > end_year:
            break

        if g.kind in (astronomy.EclipseKind.Total, astronomy.EclipseKind.Annular):
            out.append(
                GlobalCase(
                    kind=g.kind,
                    peak=g.peak,
                    peak_jd_ut=ae_time_to_jd_ut(g.peak),
                    lat=g.latitude,
                    lon=wrap_lon_deg(g.longitude),
                )
            )

        g = astronomy.NextGlobalSolarEclipse(g.peak)

    if not out:
        raise RuntimeError("No global eclipse cases produced; check year range.")
    return out


GLOBAL_CASES = build_global_cases(1900, 2100)


def test_global_peak_points_match_kind_and_peak_time() -> None:
    """
    Deterministic regression over hundreds of known eclipses:
    For each total/annular global eclipse, evaluate your solver at the reported peak shadow center.

    Uses NASA's own peak coordinates (lat_ge, lon_ge) when available to avoid
    ephemeris mismatches between Astronomy Engine and NASA Besselian elements.
    """
    max_dt_s = 0.0

    for c in GLOBAL_CASES:
        nasa = nearest_nasa_record(c.peak_jd_ut)
        # Sanity: we must be matching the same eclipse by date/time.
        assert abs(nasa_jd_ut_ge(nasa) - c.peak_jd_ut) < 1.0  # within 1 day

        # Use NASA's peak coordinates if available (more accurate with NASA Besselian elements)
        lat = nasa.lat_ge if nasa.lat_ge is not None else c.lat
        lon = nasa.lon_ge if nasa.lon_ge is not None else c.lon

        loc = Location(latitude_deg=lat, longitude_deg_east=lon, altitude_m=0.0)
        res = compute_local_eclipse(nasa, loc)

        assert res.kind == map_ae_kind(c.kind)

        dt_s = abs(res.jd_ut_max - c.peak_jd_ut) * 86400.0
        max_dt_s = max(max_dt_s, dt_s)

        # Loose but useful: if you're hours off, something is deeply wrong (e.g., longitude sign).
        assert dt_s < 600.0  # 10 minutes

    # Optional: keep an eye on worst-case drift as you tweak formulas.
    assert max_dt_s < 600.0


@st.composite
def global_case_with_offset(draw: st.DrawFn) -> tuple[GlobalCase, float, float]:
    c = draw(st.sampled_from(GLOBAL_CASES))
    east_km = draw(st.floats(min_value=-600.0, max_value=600.0))
    north_km = draw(st.floats(min_value=-600.0, max_value=600.0))
    lat, lon = offset_lat_lon_km(c.lat, c.lon, east_km=east_km, north_km=north_km)
    assume(-89.9 <= lat <= 89.9)
    return c, lat, lon


@given(global_case_with_offset())
@settings(max_examples=400, deadline=None)
def test_local_kind_matches_astronomy_engine_near_center(
    case_and_loc: tuple[GlobalCase, float, float],
) -> None:
    """
    Property-based test:
    jitter locations near the global shadow center and compare local kind to SearchLocalSolarEclipse.
    """
    c, lat, lon = case_and_loc

    obs = astronomy.Observer(lat, lon, 0.0)

    # Search for the local eclipse near this global eclipse.
    # SearchLocalSolarEclipse can return events that are partly/completely invisible due to time of day.
    info = astronomy.SearchLocalSolarEclipse(c.peak.AddDays(-7.0), obs)

    oracle_peak_jd = ae_time_to_jd_ut(info.peak.time)

    # Ensure we’re talking about the same eclipse event (not the next one months later).
    assume(abs(oracle_peak_jd - c.peak_jd_ut) < 0.8)  # within ~19 hours

    # Require Sun above horizon at peak to avoid “night-side” theoretical events.
    assume(info.peak.altitude > 1.0)

    # Avoid vanishingly small partials near the penumbral edge (classification-sensitive).
    assume(info.obscuration > 0.02)

    nasa = nearest_nasa_record(oracle_peak_jd)
    assert abs(nasa_jd_ut_ge(nasa) - oracle_peak_jd) < 1.0

    # Skip cases where NASA and Astronomy Engine disagree about the eclipse type.
    # This happens for some edge cases (hybrid eclipses, narrow paths, etc.).
    nasa_kind = {"T": "total", "H": "total", "A": "annular", "P": "partial"}.get(
        nasa.main_type, "partial"
    )
    assume(nasa_kind == map_ae_kind(info.kind))

    loc = Location(latitude_deg=lat, longitude_deg_east=lon, altitude_m=0.0)
    res = compute_local_eclipse(nasa, loc)

    # Skip boundary cases where observer is near the umbral/antumbral edge.
    # Small ephemeris differences can flip the classification in these cases.
    if info.kind in (astronomy.EclipseKind.Total, astronomy.EclipseKind.Annular):
        # Check if m is within 5% of |L2p| (near the shadow edge)
        t_max, fa = solve_t_max(nasa, loc)
        m = math.hypot(fa.u, fa.v)
        assume(abs(m - abs(fa.L2p)) / max(abs(fa.L2p), 1e-9) > 0.05)

    assert res.kind == map_ae_kind(info.kind)

    # Timing check: again loose, but catches sign/units errors.
    assert abs(res.jd_ut_max - oracle_peak_jd) * 86400.0 < 600.0

    # If total, compare totality endpoints too.
    if info.kind == astronomy.EclipseKind.Total:
        assume(info.total_begin is not None)
        assume(info.total_end is not None)
        assume(info.total_begin.altitude > 1.0)
        assume(info.total_end.altitude > 1.0)

        oracle_c2 = ae_time_to_jd_ut(info.total_begin.time)
        oracle_c3 = ae_time_to_jd_ut(info.total_end.time)

        assert res.jd_ut_c2 is not None
        assert res.jd_ut_c3 is not None

        assert abs(res.jd_ut_c2 - oracle_c2) * 86400.0 < 600.0
        assert abs(res.jd_ut_c3 - oracle_c3) * 86400.0 < 600.0
