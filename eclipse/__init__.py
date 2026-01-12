"""
Eclipse - Solar eclipse calculations using NASA Besselian elements.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
Data source: https://eclipse.gsfc.nasa.gov/eclipse_besselian_from_mysqldump2.csv

Example usage:
    from eclipse import Location, find_next_total_eclipse

    loc = Location(latitude_deg=40.7128, longitude_deg_east=-74.0060)
    result = find_next_total_eclipse(loc)
    if result:
        print(f"Next total eclipse: {result.eclipse.date_iso}")
"""

from __future__ import annotations

from datetime import datetime
from typing import Optional

from astropy.time import Time

from .catalog import load_catalog
from .compute import compute_local_eclipse, find_total_eclipses
from .models import (
    EclipseKind,
    EclipseRecord,
    LocalEclipseResult,
    Location,
)

__all__ = [
    "Location",
    "EclipseRecord",
    "LocalEclipseResult",
    "EclipseKind",
    "find_next_total_eclipse",
    "find_previous_total_eclipse",
    "find_total_eclipses_around",
    "compute_local_eclipse",
    "load_catalog",
]


def find_next_total_eclipse(
    location: Location,
    *,
    after: Optional[datetime] = None,
    catalog_path: Optional[str] = None,
) -> Optional[LocalEclipseResult]:
    """Find the next total solar eclipse visible from a location.

    Args:
        location: Geographic location to check.
        after: Reference datetime (default: now).
        catalog_path: Path to CSV catalog file (default: bundled catalog).

    Returns:
        LocalEclipseResult if found, None otherwise.
    """
    catalog = load_catalog(catalog_path)
    ref_time = Time(after) if after else Time.now()
    _, next_eclipse = find_total_eclipses(catalog, location, ref_time.jd)
    return next_eclipse


def find_previous_total_eclipse(
    location: Location,
    *,
    before: Optional[datetime] = None,
    catalog_path: Optional[str] = None,
) -> Optional[LocalEclipseResult]:
    """Find the most recent total solar eclipse visible from a location.

    Args:
        location: Geographic location to check.
        before: Reference datetime (default: now).
        catalog_path: Path to CSV catalog file (default: bundled catalog).

    Returns:
        LocalEclipseResult if found, None otherwise.
    """
    catalog = load_catalog(catalog_path)
    ref_time = Time(before) if before else Time.now()
    prev_eclipse, _ = find_total_eclipses(catalog, location, ref_time.jd)
    return prev_eclipse


def find_total_eclipses_around(
    location: Location,
    *,
    reference: Optional[datetime] = None,
    catalog_path: Optional[str] = None,
) -> tuple[Optional[LocalEclipseResult], Optional[LocalEclipseResult]]:
    """Find both previous and next total solar eclipses visible from a location.

    Args:
        location: Geographic location to check.
        reference: Reference datetime (default: now).
        catalog_path: Path to CSV catalog file (default: bundled catalog).

    Returns:
        Tuple of (previous_eclipse, next_eclipse). Either may be None.
    """
    catalog = load_catalog(catalog_path)
    ref_time = Time(reference) if reference else Time.now()
    return find_total_eclipses(catalog, location, ref_time.jd)
