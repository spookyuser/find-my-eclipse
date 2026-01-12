"""
Command-line interface for eclipse calculations.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
"""

from __future__ import annotations

import argparse
from typing import Optional

from astropy.time import Time

from .catalog import load_catalog
from .compute import find_total_eclipses
from .models import LocalEclipseResult, Location


def print_result(label: str, r: Optional[LocalEclipseResult]) -> None:
    if r is None:
        print(f"{label}: none found in catalog range")
        return

    e = r.eclipse
    print(f"{label}: cat_no={e.cat_no}  date={e.date_iso}  type={e.eclipse_type_raw}")
    print(
        f"  max: {Time(r.jd_ut_max, format='jd', scale='utc').iso} UTC  (kind={r.kind})"
    )

    if r.jd_ut_c2 is not None and r.jd_ut_c3 is not None:
        print(f"  C2 : {Time(r.jd_ut_c2, format='jd', scale='utc').iso} UTC")
        print(f"  C3 : {Time(r.jd_ut_c3, format='jd', scale='utc').iso} UTC")
        print(f"  duration: {r.duration_seconds:.1f} s")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Find total solar eclipses visible from a location."
    )
    ap.add_argument(
        "--csv",
        default=None,
        help="Path to eclipse CSV catalog (default: bundled catalog)",
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

    prev_r, next_r = find_total_eclipses(catalog, loc, ref_jd)

    print(
        f"Location: lat={loc.latitude_deg:.6f}, lon_east={loc.longitude_deg_east:.6f}, alt_m={loc.altitude_m:.1f}"
    )
    print(f"Reference: {ref_time.iso} UTC (JD_UT={ref_jd:.6f})")
    print_result("Previous total eclipse", prev_r)
    print_result("Next total eclipse", next_r)


if __name__ == "__main__":
    main()
