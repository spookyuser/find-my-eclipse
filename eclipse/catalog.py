"""
Eclipse catalog loading and management.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
Data source: https://eclipse.gsfc.nasa.gov/eclipse_besselian_from_mysqldump2.csv
"""

from __future__ import annotations

import csv
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Optional

from .models import BesselianCoefficients, EclipseRecord


def load_catalog(csv_path: Optional[str] = None) -> list[EclipseRecord]:
    """Load eclipse catalog from CSV file.

    Args:
        csv_path: Path to CSV file. If None, uses bundled catalog data.

    Returns:
        List of EclipseRecord objects.
    """
    if csv_path is None:
        return _load_bundled_catalog()

    p = Path(csv_path)
    if not p.exists():
        raise FileNotFoundError(p)

    with p.open("r", encoding="utf-8", newline="") as f:
        return _parse_catalog(f)


@lru_cache(maxsize=1)
def _load_bundled_catalog() -> list[EclipseRecord]:
    """Load the bundled NASA eclipse catalog."""
    try:
        files = resources.files("eclipse.data")
        csv_file = files.joinpath("eclipse_besselian.csv")
        with resources.as_file(csv_file) as path:
            with path.open("r", encoding="utf-8", newline="") as f:
                return _parse_catalog(f)
    except (FileNotFoundError, TypeError):
        raise FileNotFoundError(
            "Bundled catalog not found. Please provide a csv_path argument."
        )


def _parse_catalog(f) -> list[EclipseRecord]:
    """Parse eclipse catalog from a file-like object."""
    records: list[EclipseRecord] = []
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
        f_vals = {
            k: float(row[k])
            for k in row
            if k
            not in ("eclipse_type", "td_ge", "lat_ge", "lng_ge", "central_duration")
        }
        b = BesselianCoefficients(
            t0_tdt_hours=f_vals["t0"],
            x0=f_vals["x0"],
            x1=f_vals["x1"],
            x2=f_vals["x2"],
            x3=f_vals["x3"],
            y0=f_vals["y0"],
            y1=f_vals["y1"],
            y2=f_vals["y2"],
            y3=f_vals["y3"],
            d0=f_vals["d0"],
            d1=f_vals["d1"],
            d2=f_vals["d2"],
            mu0=f_vals["mu0"],
            mu1=f_vals["mu1"],
            mu2=f_vals["mu2"],
            l10=f_vals["l10"],
            l11=f_vals["l11"],
            l12=f_vals["l12"],
            l20=f_vals["l20"],
            l21=f_vals["l21"],
            l22=f_vals["l22"],
            tan_f1=f_vals["tan_f1"],
            tan_f2=f_vals["tan_f2"],
            tmin=f_vals["tmin"],
            tmax=f_vals["tmax"],
        )

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
                year=int(f_vals["year"]),
                month=int(f_vals["month"]),
                day=int(f_vals["day"]),
                eclipse_type=raw_type,
                main_type=EclipseRecord.parse_main_type(raw_type),
                dt=f_vals["dt"],
                julian_date=f_vals["julian_date"],
                td_ge=row["td_ge"],
                cat_no=int(f_vals["cat_no"]),
                bessel=b,
                lat_ge=lat_ge,
                lon_ge=lon_ge,
            )
        )

    return records
