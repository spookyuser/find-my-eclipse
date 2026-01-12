"""
Eclipse calculation functions using Besselian elements.

Attribution: "Eclipse Predictions by Fred Espenak, NASA's GSFC"
"""

from __future__ import annotations

import math
from typing import Literal, Optional, Tuple

from .models import (
    EclipseKind,
    EclipseRecord,
    FundamentalArgs,
    LocalEclipseResult,
    Location,
)

RAD = math.pi / 180.0

EARTH_EQ_RADIUS_M = 6_378_140.0
EARTH_1ME2_SQRT = 0.99664719
DEG_PER_SEC = 0.00417807


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


def find_total_eclipses(
    catalog: list[EclipseRecord],
    loc: Location,
    ref_jd_ut: float,
) -> Tuple[Optional[LocalEclipseResult], Optional[LocalEclipseResult]]:
    """Find the previous and next total solar eclipses visible from a location.

    Args:
        catalog: List of eclipse records to search.
        loc: Geographic location.
        ref_jd_ut: Reference Julian Date (UT) to search around.

    Returns:
        Tuple of (previous_eclipse, next_eclipse). Either may be None if not found.
    """
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
