"""
Web API for eclipse calculations.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from astropy.time import Time
from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .catalog import load_catalog
from .compute import find_total_eclipses
from .models import Location

app = FastAPI(title="Eclipse Finder")

CATALOG = load_catalog()
STATIC_DIR = Path(__file__).parent / "static"


class EclipseInfo(BaseModel):
    cat_no: int
    date: str
    eclipse_type: str
    max_time_utc: str
    c2_time_utc: Optional[str] = None
    c3_time_utc: Optional[str] = None
    duration_seconds: Optional[float] = None


class EclipseResponse(BaseModel):
    lat: float
    lon: float
    ref_time_utc: str
    previous: Optional[EclipseInfo] = None
    next: Optional[EclipseInfo] = None


@app.get("/api/eclipses")
def get_eclipses(lat: float, lon: float, ref_utc: Optional[str] = None) -> EclipseResponse:
    if ref_utc:
        ref_time = Time(ref_utc, scale="utc")
    else:
        ref_time = Time.now()

    loc = Location(latitude_deg=lat, longitude_deg_east=lon, altitude_m=0)
    prev_r, next_r = find_total_eclipses(CATALOG, loc, ref_time.jd)

    def to_info(r) -> Optional[EclipseInfo]:
        if r is None:
            return None
        return EclipseInfo(
            cat_no=r.eclipse.cat_no,
            date=r.eclipse.date_iso,
            eclipse_type=r.eclipse.eclipse_type_raw,
            max_time_utc=Time(r.jd_ut_max, format="jd", scale="utc").iso,
            c2_time_utc=Time(r.jd_ut_c2, format="jd", scale="utc").iso if r.jd_ut_c2 else None,
            c3_time_utc=Time(r.jd_ut_c3, format="jd", scale="utc").iso if r.jd_ut_c3 else None,
            duration_seconds=r.duration_seconds,
        )

    return EclipseResponse(
        lat=lat,
        lon=lon,
        ref_time_utc=ref_time.iso,
        previous=to_info(prev_r),
        next=to_info(next_r),
    )


@app.get("/")
def index():
    return FileResponse(STATIC_DIR / "index.html")


if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


def run():
    import uvicorn
    uvicorn.run("eclipse.web:app", host="127.0.0.1", port=8000, reload=True)
