#!/usr/bin/env python3
"""
Speed-of-sound validation (EDMD or TIME batch outputs).

Reads the speed-of-sound batch CSVs produced by 00ALLINONE:
  wall_x_positions_L0_{L0x10}_wallmassfactor_{M}_run{r}.csv

For each run:
  - extract divider displacement vs time after release
  - estimate fundamental oscillation frequency via FFT peak

For each (L0, M):
  - aggregate across repeats (mean + stderr)

For each L0:
  - compute K (fundamental root of cot(K)=alpha*K with alpha=M/(2N_side))
  - fit nu vs x where x = K / (2*pi*L_eff)
  - slope = c_s with error bar from weighted regression

Outputs:
  - speed_of_sound_runs.csv
  - speed_of_sound_groups.csv
  - speed_of_sound_summary.csv
  - speed_of_sound_cs_vs_eta.pdf
  - speed_of_sound_per_L0_with_eta.pdf
  - FINAL speed_of_sound_on_packing_fracture.pdf (optional)
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np


FILENAME_RE = re.compile(r"wall_x_positions_L0_(\d+)_wallmassfactor_(\d+)_run(\d+)\.csv$")

# Roman et al. (2002) reference values (as used historically in wall_x_FFT.py).
# These are useful as a "literature points" overlay in the final cs(eta) figure.
ROMAN_2002_BY_L0 = {
    7.5: (5.99, 0.09),
    10.0: (3.78, 0.08),
    15.0: (2.61, 0.03),
    20.0: (2.20, 0.02),
    25.0: (2.01, 0.02),  # ##CHRIS 2026-09-13: Roman et al. 2002 Table I gives 2.01 +- 0.02 (was mistranscribed as 2.10)
    30.0: (1.89, 0.02),
    35.0: (1.81, 0.02),
}


@dataclass(frozen=True)
class RunKey:
    l0: float
    wall_mass_factor: int
    run: int


def find_latest_sim_dir(base_dir: Path, prefix: str) -> Optional[Path]:
    if not base_dir.exists():
        return None
    candidates = [p for p in base_dir.iterdir() if p.is_dir() and p.name.startswith(prefix)]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def parse_run_key(path: Path) -> Optional[RunKey]:
    m = FILENAME_RE.search(path.name)
    if not m:
        return None
    l0 = int(m.group(1)) / 10.0
    mf = int(m.group(2))
    run = int(m.group(3))
    return RunKey(l0=l0, wall_mass_factor=mf, run=run)


def read_wall_csv_5cols(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Robust CSV reader that only keeps the first 5 columns:
      Time, Wall_X, Displacement(σ), Left_Count, Right_Count
    Skips malformed rows.
    """
    t: list[float] = []
    wall_x: list[float] = []
    disp: list[float] = []
    left: list[int] = []
    right: list[int] = []

    with path.open("r", newline="") as f:
        rdr = csv.reader(f)
        # header
        try:
            next(rdr)
        except StopIteration:
            return np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for row in rdr:
            if len(row) < 5:
                continue
            try:
                tt = float(row[0].strip())
                wx = float(row[1].strip())
                dd = float(row[2].strip())
                ll = int(float(row[3].strip()))
                rr = int(float(row[4].strip()))
            except Exception:
                continue
            if not math.isfinite(tt) or not math.isfinite(wx) or not math.isfinite(dd):
                continue
            t.append(tt)
            wall_x.append(wx)
            disp.append(dd)
            left.append(ll)
            right.append(rr)

    return (
        np.asarray(t, dtype=float),
        np.asarray(wall_x, dtype=float),
        np.asarray(disp, dtype=float),
        np.asarray(left, dtype=int),
        np.asarray(right, dtype=int),
    )


def read_wall_csv_metadata(path: Path) -> dict[str, str]:
    """Read the first data row as metadata without disturbing legacy 5-column input."""
    try:
        with path.open("r", newline="") as f:
            row = next(csv.DictReader(f), None)
        return dict(row or {})
    except Exception:
        return {}


def dominant_frequency_fft(
    t: np.ndarray,
    disp: np.ndarray,
    *,
    drop_frac: float,
    f_min: float,
    f_max: Optional[float],
    peak_signal: str,
) -> Tuple[float, float, float, float, float, float, int, float]:
    """
    Returns (nu_peak, peak_power, dt, span, df, peak_snr,
             peak_on_search_boundary, raw_bin_frequency).
    span: time span used for FFT (after transient drop), in σ-time.
    df: frequency resolution ~ 1/span.
    """
    if len(t) < 8:
        return (float("nan"),) * 6 + (1, float("nan"))

    # Ensure increasing time
    order = np.argsort(t)
    t = t[order]
    disp = disp[order]

    # Drop initial transient fraction
    if drop_frac > 0:
        k0 = int(len(t) * drop_frac)
        if k0 >= len(t) - 8:
            k0 = max(0, len(t) - 8)
        t = t[k0:] - t[k0]
        disp = disp[k0:]

    dt = float(np.median(np.diff(t))) if len(t) > 1 else float("nan")
    if not math.isfinite(dt) or dt <= 0:
        return (float("nan"),) * 6 + (1, float("nan"))

    span = float(t[-1] - t[0]) if len(t) > 1 else float("nan")
    df = (1.0 / span) if (math.isfinite(span) and span > 0) else float("nan")

    # IMPORTANT:
    # For some configs (especially very light walls), displacement can have strong low-frequency
    # drift/random-walk components that dominate the PSD near f_min and cause "snapping" to f_min.
    # Using the velocity signal (d/dt displacement) suppresses that drift and yields a much cleaner peak.
    def _psd_peak(tt: np.ndarray, yy: np.ndarray) -> Tuple[float, float, float, int, float]:
        xx = np.asarray(tt, dtype=float)
        y0 = np.asarray(yy, dtype=float)
        if len(xx) >= 2:
            p = np.polyfit(xx, y0, deg=1)
            y1 = y0 - (p[0] * xx + p[1])
        else:
            y1 = y0 - float(np.mean(y0))
        window = np.hanning(len(y1))
        y1 = y1 * window
        spec = np.fft.rfft(y1)
        power = (spec.real * spec.real + spec.imag * spec.imag) / max(1.0, float(len(y1)))
        freqs = np.fft.rfftfreq(len(y1), d=float(np.median(np.diff(xx))) if len(xx) > 1 else dt)
        mask = freqs >= f_min
        if f_max is not None:
            mask &= freqs <= f_max
        if not np.any(mask):
            return float("nan"), float("nan"), float("nan"), 1, float("nan")
        valid_indices = np.flatnonzero(mask)
        local_i = int(np.argmax(power[mask]))
        i = int(valid_indices[local_i])
        raw_nu = float(freqs[i])
        pk = float(power[i])

        # Quadratic interpolation of log power gives a sub-bin peak estimate
        # without pretending that the finite-record resolution df disappeared.
        nu = raw_nu
        if 0 < i < len(power) - 1 and power[i - 1] > 0 and power[i] > 0 and power[i + 1] > 0:
            a = math.log(float(power[i - 1]))
            b = math.log(float(power[i]))
            c = math.log(float(power[i + 1]))
            denom = a - 2.0 * b + c
            if math.isfinite(denom) and abs(denom) > 1e-15:
                delta = 0.5 * (a - c) / denom
                delta = float(np.clip(delta, -0.5, 0.5))
                nu = raw_nu + delta * float(freqs[1] - freqs[0])

        noise_mask = mask.copy()
        noise_mask[max(0, i - 2):min(len(noise_mask), i + 3)] = False
        noise = power[noise_mask]
        noise = noise[np.isfinite(noise) & (noise > 0)]
        noise_floor = float(np.median(noise)) if len(noise) else float("nan")
        peak_snr = pk / noise_floor if math.isfinite(noise_floor) and noise_floor > 0 else float("nan")
        on_boundary = int(local_i == 0 or local_i == len(valid_indices) - 1)
        return nu, pk, peak_snr, on_boundary, raw_nu

    # displacement peak
    nu_d, pk_d, snr_d, boundary_d, raw_d = _psd_peak(t, disp)

    # velocity peak (centered time stamps)
    if len(t) >= 3:
        dt_local = np.diff(t)
        good = dt_local > 0
        if np.any(good):
            v = np.diff(disp)[good] / dt_local[good]
            t_v = 0.5 * (t[1:][good] + t[:-1][good])
            nu_v, pk_v, snr_v, boundary_v, raw_v = _psd_peak(t_v, v)
        else:
            nu_v, pk_v, snr_v, boundary_v, raw_v = (
                float("nan"), float("nan"), float("nan"), 1, float("nan")
            )
    else:
        nu_v, pk_v, snr_v, boundary_v, raw_v = (
            float("nan"), float("nan"), float("nan"), 1, float("nan")
        )

    # choose signal
    peak_signal = str(peak_signal).strip().lower()
    if peak_signal in ("vel", "velocity"):
        nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_v, pk_v, snr_v, boundary_v, raw_v
    elif peak_signal in ("disp", "displacement"):
        nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_d, pk_d, snr_d, boundary_d, raw_d
    else:
        # "best": default heuristic:
        # - if displacement peak snaps close to f_min (typical for drift/random-walk), prefer velocity
        # - otherwise prefer the larger peak_power
        if math.isfinite(nu_v) and math.isfinite(pk_v):
            if math.isfinite(nu_d) and math.isfinite(pk_d):
                if math.isfinite(df) and df > 0 and nu_d <= float(f_min) + 3.0 * df and nu_v > float(f_min) + 3.0 * df:
                    nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_v, pk_v, snr_v, boundary_v, raw_v
                else:
                    if pk_d >= pk_v:
                        nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_d, pk_d, snr_d, boundary_d, raw_d
                    else:
                        nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_v, pk_v, snr_v, boundary_v, raw_v
            else:
                nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_v, pk_v, snr_v, boundary_v, raw_v
        else:
            nu_peak, peak_power, peak_snr, peak_boundary, raw_bin = nu_d, pk_d, snr_d, boundary_d, raw_d

    return nu_peak, peak_power, dt, span, df, peak_snr, peak_boundary, raw_bin


def _debug_psd_and_trace(
    csv_path: Path,
    out_dir: Path,
    *,
    drop_frac: float,
    f_min: float,
    f_max: Optional[float],
) -> Path:
    """
    Write a diagnostic plot for one run:
      - displacement vs time (full + post-transient)
      - power spectrum (FFT of detrended, Hann-windowed displacement)
    """
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        raise SystemExit(f"matplotlib required for debug plots: {e}")

    t, wall_x, disp, left, right = read_wall_csv_5cols(csv_path)
    if len(t) < 16:
        raise SystemExit(f"Not enough samples in {csv_path}")

    # Ensure increasing time
    order = np.argsort(t)
    t = t[order]
    wall_x = wall_x[order]
    disp = disp[order]

    # Drop transient
    n = len(t)
    k0 = int(max(0, min(n - 8, math.floor(drop_frac * n))))
    t2 = t[k0:] - t[k0]
    d2 = disp[k0:]
    if len(t2) < 16:
        raise SystemExit(f"Not enough post-transient samples in {csv_path} (drop_frac={drop_frac})")

    dt = float(np.median(np.diff(t2)))
    span = float(t2[-1] - t2[0])
    df = float(1.0 / span) if span > 0 else float("nan")

    # Detrend (linear) then window (Hann) to reduce leakage at low frequencies.
    x = np.asarray(t2, dtype=float)
    y = np.asarray(d2, dtype=float)
    p = np.polyfit(x, y, deg=1)
    y_detr = y - (p[0] * x + p[1])
    y_win = y_detr * np.hanning(len(y_detr))

    freqs = np.fft.rfftfreq(len(y_win), d=dt)
    spec = np.fft.rfft(y_win)
    power = (spec.real * spec.real + spec.imag * spec.imag) / max(1.0, float(len(y_win)))

    mask = freqs >= float(f_min)
    if f_max is not None:
        mask &= freqs <= float(f_max)
    if not np.any(mask):
        raise SystemExit(f"No frequencies in band for {csv_path} (f_min={f_min}, f_max={f_max})")
    idx = int(np.argmax(power[mask]))
    nu = float(freqs[mask][idx])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 6.8), constrained_layout=True)

    ax1.plot(t - t[0], disp, lw=1.0, alpha=0.30, label="disp (full)")
    ax1.plot(t2, d2, lw=1.2, label=f"disp (post drop {drop_frac:.0%})")
    ax1.axvline(0.0, color="k", lw=0.8, alpha=0.4)
    ax1.set_xlabel("Time after release (σ-time)")
    ax1.set_ylabel("Divider displacement (σ)")
    ax1.set_title(csv_path.stem)
    ax1.grid(True, alpha=0.25)
    # show count constancy
    left0 = int(left[0]) if len(left) else 0
    right0 = int(right[0]) if len(right) else 0
    counts_ok = bool(len(left)) and bool(len(right)) and bool(np.all(left == left0)) and bool(np.all(right == right0))
    ax1.text(
        0.01,
        0.98,
        f"N_L={left0} N_R={right0} (constant={int(counts_ok)})",
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.85),
    )
    ax1.legend(loc="best", fontsize=9)

    ax2.plot(freqs, power, lw=1.0)
    ax2.axvline(nu, color="r", lw=1.1, label=f"peak ν={nu:.6g}")
    xlo = float(f_min)
    xhi = float(f_max) if f_max is not None else float(freqs[-1])
    ax2.set_xlim(xlo, xhi)
    ax2.set_xlabel("Frequency ν (1/σ-time)")
    ax2.set_ylabel("Power (a.u.)")
    ax2.set_title(f"PSD (dt≈{dt:.4g}, span≈{span:.4g}, df≈{df:.4g})")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="best", fontsize=9)

    out_pdf = out_dir / f"debug_trace_psd_{csv_path.stem}.pdf"
    fig.savefig(out_pdf)
    plt.close(fig)
    return out_pdf


def k_root_bisect(alpha: float, *, iters: int = 120) -> float:
    """
    Solve cot(K) = alpha*K for the fundamental root in (0, pi).
    alpha = M/(2*N_side).
    """
    a = 1e-4
    b = math.pi - 1e-4

    def f(x: float) -> float:
        return 1.0 / math.tan(x) - alpha * x

    fa = f(a)
    fb = f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        # Fallback bracket (works for extremely small alpha)
        a = 0.1
        b = math.pi / 2.0 - 1e-4
        fa = f(a)
        fb = f(b)
        if fa * fb > 0:
            raise ValueError("Could not bracket K root")

    for _ in range(iters):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0:
            return m
        if fa * fm <= 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
    return 0.5 * (a + b)


def weighted_linreg(
    x: np.ndarray,
    y: np.ndarray,
    sy: np.ndarray,
    *,
    force_zero_intercept: bool,
) -> Tuple[float, float, float, float]:
    """
    Weighted least squares for y = a*x + b.
    Returns (a, sa, b, sb).
    """
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan"), float("nan")
    w = 1.0 / np.maximum(sy, 1e-12) ** 2

    if force_zero_intercept:
        sxx = float(np.sum(w * x * x))
        sxy = float(np.sum(w * x * y))
        a = sxy / sxx if sxx > 0 else float("nan")
        b = 0.0
        resid = y - a * x
        dof = max(1, len(x) - 1)
        chi2 = float(np.sum(w * resid * resid))
        s2 = chi2 / dof
        sa = math.sqrt(s2 / sxx) if sxx > 0 else float("nan")
        sb = 0.0
        return a, sa, b, sb

    S = float(np.sum(w))
    Sx = float(np.sum(w * x))
    Sy = float(np.sum(w * y))
    Sxx = float(np.sum(w * x * x))
    Sxy = float(np.sum(w * x * y))
    Delta = S * Sxx - Sx * Sx
    if Delta <= 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    a = (S * Sxy - Sx * Sy) / Delta
    b = (Sxx * Sy - Sx * Sxy) / Delta
    resid = y - (a * x + b)
    dof = max(1, len(x) - 2)
    chi2 = float(np.sum(w * resid * resid))
    s2 = chi2 / dof
    var_a = s2 * (S / Delta)
    var_b = s2 * (Sxx / Delta)
    sa = math.sqrt(max(var_a, 0.0))
    sb = math.sqrt(max(var_b, 0.0))
    return a, sa, b, sb


def _mad_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return float("nan")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    if mad <= 0:
        return 0.0
    return 1.4826 * mad


def robust_filter_points(
    x: np.ndarray,
    y: np.ndarray,
    sy: np.ndarray,
    *,
    mad_z: float,
    min_keep: int,
) -> Tuple[np.ndarray, dict]:
    """
    Robustly filter out outlier (x,y) points before fitting.

    We assume the physical model is approximately linear in this plot:
      nu ≈ c_s * x   (possibly with a tiny intercept)

    Outliers in practice often come from a wrong FFT peak pick (e.g. snapping to f_min),
    which produces an implied slope c_i = y_i/x_i orders of magnitude off.

    Strategy:
      - compute implied slopes c_i = y_i / x_i (x>0)
      - MAD-clip c_i around its median
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sy, dtype=float)
    n = len(x)
    if n == 0:
        return np.zeros((0,), dtype=bool), dict(n_total=0, n_used=0, n_dropped=0)

    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sy) & (x > 0) & (sy > 0)
    if np.sum(keep) < max(2, min_keep):
        # Not enough points to do anything robust; return the finite subset.
        out = dict(n_total=n, n_used=int(np.sum(keep)), n_dropped=int(n - np.sum(keep)))
        return keep, out

    ci = np.zeros_like(y, dtype=float)
    ci[keep] = y[keep] / x[keep]
    ci_sub = ci[keep]
    ci_med = float(np.median(ci_sub))
    sigma = _mad_sigma(ci_sub)
    if not math.isfinite(sigma) or sigma <= 0:
        out = dict(
            n_total=n,
            n_used=int(np.sum(keep)),
            n_dropped=int(n - np.sum(keep)),
            ci_med=ci_med,
            ci_sigma=float("nan"),
        )
        return keep, out

    z = (ci - ci_med) / sigma
    keep2 = keep & (np.abs(z) <= max(0.0, float(mad_z)))

    # Ensure we still have enough points to fit; if not, progressively relax.
    if int(np.sum(keep2)) < max(2, min_keep):
        keep2 = keep.copy()
    out = dict(
        n_total=n,
        n_used=int(np.sum(keep2)),
        n_dropped=int(n - np.sum(keep2)),
        ci_med=ci_med,
        ci_sigma=sigma,
    )
    return keep2, out


def cs_theory_family(eta: np.ndarray, *, a: float, kbt: float, m: float) -> np.ndarray:
    # Matches the Henderson/SPT family used in wall_x_FFT.py for the "FINAL" plot.
    term = (1.0 + eta + 3.0 * a * eta * eta - a * eta * eta * eta) / np.maximum(1e-12, (1.0 - eta) ** 3)
    return np.sqrt(2.0 * kbt / m) * np.sqrt(np.maximum(term, 0.0))


def Z_henderson_family(eta: np.ndarray, *, a: float) -> np.ndarray:
    """Legacy Román Eq. (26) sound factor, retained for old result files.

    Despite the historical function name this is not the compressibility
    factor Z itself; it is [Z + eta*dZ/deta] for
    Z=(1+a*eta^2)/(1-eta)^2.
    """
    return (1.0 + eta + 3.0 * a * eta * eta - a * eta * eta * eta) / np.maximum(1e-12, (1.0 - eta) ** 3)


def dZ_henderson_family(eta: np.ndarray, *, a: float) -> np.ndarray:
    # d/dη of (N/D) with N=1+η+3aη²-aη³, D=(1-η)³
    N = 1.0 + eta + 3.0 * a * eta * eta - a * eta * eta * eta
    dN = 1.0 + 6.0 * a * eta - 3.0 * a * eta * eta
    D = np.maximum(1e-12, (1.0 - eta) ** 3)
    dD = -3.0 * np.maximum(1e-12, (1.0 - eta) ** 2)
    return (dN * D - N * dD) / np.maximum(1e-12, D * D)


def Z_spt(eta: np.ndarray) -> np.ndarray:
    """Legacy Román Eq. (26) sound factor for SPT (a=0).

    The historical name is retained for compatibility.  This is
    ``Z_EOS + eta*dZ_EOS/deta``, not the SPT compressibility factor itself.
    """
    eta = np.asarray(eta, dtype=float)
    return (1.0 + eta) / np.maximum(1e-12, (1.0 - eta) ** 3)


def dZ_spt(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    one_minus = np.maximum(1e-12, 1.0 - eta)
    return 1.0 / (one_minus ** 3) + 3.0 * (1.0 + eta) / (one_minus ** 4)


def cs_simple_from_Z(Z: np.ndarray, *, kbt: float, m: float) -> np.ndarray:
    # Legacy/simple mapping used in older plots: c_s = sqrt(2 kBT/m) * sqrt(Z)
    return np.sqrt(2.0 * kbt / m) * np.sqrt(np.maximum(Z, 0.0))


def cs_adiabatic_2d_monatomic(Z: np.ndarray, dZ: np.ndarray, eta: np.ndarray, *, kbt: float, m: float) -> np.ndarray:
    # More physically consistent adiabatic mapping used in wall_x_FFT.py:
    # c_s^2 ≈ (kBT/m) * [ Z(η) + η Z'(η) + Z(η)^2 ]
    cs2 = (kbt / m) * np.maximum(Z + eta * dZ + Z * Z, 0.0)
    return np.sqrt(cs2)


def Z_henderson_eos(eta: np.ndarray, *, a: float = 0.125) -> np.ndarray:
    """Hard-disk compressibility factor used by Román et al., Eq. (25)."""
    eta = np.asarray(eta, dtype=float)
    return (1.0 + a * eta * eta) / np.maximum(1e-12, (1.0 - eta) ** 2)


def dZ_henderson_eos(eta: np.ndarray, *, a: float = 0.125) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    one_minus = np.maximum(1e-12, 1.0 - eta)
    return (2.0 * a * eta) / (one_minus ** 2) + 2.0 * (1.0 + a * eta * eta) / (one_minus ** 3)


def Z_spt_eos(eta: np.ndarray) -> np.ndarray:
    """Scaled-particle-theory hard-disk compressibility factor."""
    eta = np.asarray(eta, dtype=float)
    return 1.0 / np.maximum(1e-12, (1.0 - eta) ** 2)


def dZ_spt_eos(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    return 2.0 / np.maximum(1e-12, (1.0 - eta) ** 3)


# Kolafa & Rottner (2006), rho_max=0.90 fit.  Their x is eta/(1-eta).
# The source cautions that the final reduced-density interval rho=0.89..0.90,
# equivalent to eta≈0.699..0.707, may retain finite-size effects.
KR2006_COEFFICIENTS = {
    0: 1.0,
    1: 2.0,
    2: 1.12801775,
    3: 0.00181895291,
    4: -0.0526134737,
    5: 0.0504960168,
    6: -0.0325537792,
    7: 0.0134578632,
    8: 0.00140888182,
    9: -0.00834273601,
    10: 0.00694127367,
    11: -0.00262254723,
    12: 0.000355746352,
    22: -5.24672938e-9,
    57: 5.88054639e-23,
}
KR2006_ETA_MAX = math.pi * 0.90 / 4.0
# Stay just below the liquid side of the accepted liquid/hexatic coexistence
# interval.  The EOS paper itself reaches eta≈0.707, but its final interval is
# flagged as potentially affected by finite-size effects.
# Sound speed depends on a derivative of the pressure fit, which is less robust
# than Z itself near the transition.  Use a conservative precritical cutoff
# rho*=4*eta/pi≈0.879 rather than the pressure fit's eta≈0.707 endpoint.
KR2006_PLOT_ETA_MAX = 0.690


def Z_kolafa_rottner_2006(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    x = eta / np.maximum(1e-12, 1.0 - eta)
    out = np.zeros_like(x)
    for power, coefficient in KR2006_COEFFICIENTS.items():
        out += coefficient * np.power(x, power)
    return out


def dZ_kolafa_rottner_2006(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    one_minus = np.maximum(1e-12, 1.0 - eta)
    x = eta / one_minus
    dZdx = np.zeros_like(x)
    for power, coefficient in KR2006_COEFFICIENTS.items():
        if power:
            dZdx += power * coefficient * np.power(x, power - 1)
    return dZdx / (one_minus ** 2)


def cs_roman_gamma2_from_eos(
    Z: np.ndarray, dZ: np.ndarray, eta: np.ndarray, *, kbt: float, m: float
) -> np.ndarray:
    """Román Eq. (1)+(24) with their constant gamma=2 assumption."""
    cs2 = (2.0 * kbt / m) * np.maximum(Z + eta * dZ, 0.0)
    return np.sqrt(cs2)


# ---------------------------------------------------------------------------
# Hard-disk phase boundaries (thermodynamic limit).
#
# Single source of truth: the regime shading, the theory-curve masks and the
# Liu branch switch all read these.  Previously the shading carried one set of
# literals (0.700/0.716/0.720) while the KR curve was clipped at a different,
# unrelated 0.690.
#
# Values follow Bernard & Krauth (2011) / Engel et al. (2013) and are consistent
# with the transition point eta_t=0.720 used by Liu (2021).
# ---------------------------------------------------------------------------
ETA_FLUID_MAX = 0.700        # upper edge of the stable isotropic fluid
ETA_COEX_LO = 0.700          # liquid-hexatic coexistence, lower edge
ETA_COEX_HI = 0.716          # liquid-hexatic coexistence, upper edge
ETA_HEXATIC_SOLID = 0.720    # continuous hexatic->solid transition
ETA_CLOSE_PACKED = math.pi / (2.0 * math.sqrt(3.0))   # 0.9069


# ---------------------------------------------------------------------------
# Liu (2021) global hard-disk equation of state.
#   H. Liu, "Global equation of state and the phase transitions of the hard disc
#   system", Mol. Phys. 119(10) e1905897 (2021); arXiv:2010.10624.
#   Stable-fluid branch from H. Liu, Mol. Phys. 119(9) e1886364 (2021);
#   arXiv:2010.14357, Eq. (14).
#
# Unlike SPT/Henderson/Kolafa-Rottner this EOS is *phase aware*: it covers the
# stable liquid, the liquid-hexatic coexistence loop and the hexatic branch, and
# joins tangentially onto a solid branch at eta_t = 0.720.
#
# The published PDF loses superscripts under text extraction, so the fluid-branch
# coefficients here were reconstructed analytically from the paper's own integral
# Eq. (21) and then validated against four independent anchors -- see
# liu_eos_self_check() below, which reproduces all of them.
# ---------------------------------------------------------------------------
LIU2021_B1 = -1.04191e8
LIU2021_B2 = 2.66813e8
LIU2021_M1 = 53
LIU2021_M2 = 56
LIU2021_C = 1.0 / 0.75        # pole at the random-close-packing fit value eta=0.75
LIU2021_ETA_T = ETA_HEXATIC_SOLID
# Liu quotes the solid branch as valid for eta = 0.715 .. 0.9069.
LIU2021_SOLID_ETA_MAX = ETA_CLOSE_PACKED
# The liquid/hexatic branch has a pole at eta = 0.75; never evaluate it near there.
LIU2021_LH_ETA_MAX = 0.740


def Z_liu_virial(eta: np.ndarray) -> np.ndarray:
    """Liu's Carnahan-Starling-type EOS for the stable hard-disk fluid, Eq. (14).

    Z_v = [1 + eta^2/8 + eta^3/18 - (4/21) eta^4] / (1-eta)^2
    """
    eta = np.asarray(eta, dtype=float)
    num = 1.0 + eta**2 / 8.0 + eta**3 / 18.0 - (4.0 / 21.0) * eta**4
    return num / np.maximum(1e-12, (1.0 - eta) ** 2)


def dZ_liu_virial(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    one_minus = np.maximum(1e-12, 1.0 - eta)
    num = 1.0 + eta**2 / 8.0 + eta**3 / 18.0 - (4.0 / 21.0) * eta**4
    dnum = eta / 4.0 + eta**2 / 6.0 - (16.0 / 21.0) * eta**3
    return (dnum * one_minus + 2.0 * num) / (one_minus ** 3)


def _liu_close_term(eta: np.ndarray) -> np.ndarray:
    """(b1 eta^m1 + b2 eta^m2) / (1 - c eta).

    Factored as eta^m1 (b1 + b2 eta^3) to avoid cancelling two ~2.7e8 numbers.
    """
    eta = np.asarray(eta, dtype=float)
    f = np.power(eta, LIU2021_M1) * (LIU2021_B1 + LIU2021_B2 * eta**3)
    return f / (1.0 - LIU2021_C * eta)


def _dliu_close_term(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    f = np.power(eta, LIU2021_M1) * (LIU2021_B1 + LIU2021_B2 * eta**3)
    fp = (LIU2021_M1 * LIU2021_B1 * np.power(eta, LIU2021_M1 - 1)
          + LIU2021_M2 * LIU2021_B2 * np.power(eta, LIU2021_M2 - 1))
    g = 1.0 - LIU2021_C * eta
    return (fp * g + LIU2021_C * f) / (g * g)


def Z_liu_lh(eta: np.ndarray) -> np.ndarray:
    """Liu Eq. (9): stable-liquid + liquid-hexatic-coexistence + hexatic branch."""
    return Z_liu_virial(eta) + _liu_close_term(eta)


def dZ_liu_lh(eta: np.ndarray) -> np.ndarray:
    return dZ_liu_virial(eta) + _dliu_close_term(eta)


def _liu_alpha(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    return ETA_CLOSE_PACKED / np.maximum(1e-12, eta) - 1.0


def Z_liu_solid(eta: np.ndarray) -> np.ndarray:
    """Liu Eq. (13): solid branch, valid eta = 0.715 .. 0.9069."""
    a = _liu_alpha(eta)
    a = np.maximum(a, 1e-12)
    return 2.0 / a + 1.9 + a - 5.2 * a**2 + 114.48 * a**4


def dZ_liu_solid(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    a = np.maximum(_liu_alpha(eta), 1e-12)
    dZ_da = -2.0 / (a * a) + 1.0 - 10.4 * a + 457.92 * a**3
    da_deta = -ETA_CLOSE_PACKED / np.maximum(1e-12, eta * eta)
    return dZ_da * da_deta


def Z_liu_global(eta: np.ndarray) -> np.ndarray:
    """Liu Eq. (14): Z_lh below eta_t=0.720, Z_solid above."""
    eta = np.asarray(eta, dtype=float)
    return np.where(eta < LIU2021_ETA_T, Z_liu_lh(np.minimum(eta, LIU2021_LH_ETA_MAX)),
                    Z_liu_solid(eta))


def dZ_liu_global(eta: np.ndarray) -> np.ndarray:
    eta = np.asarray(eta, dtype=float)
    return np.where(eta < LIU2021_ETA_T, dZ_liu_lh(np.minimum(eta, LIU2021_LH_ETA_MAX)),
                    dZ_liu_solid(eta))


def liu_rigidity(eta: np.ndarray) -> np.ndarray:
    """Liu Eq. (10): omega = d(Z eta)/d eta = Z + eta dZ/d eta.

    Negative omega marks the van der Waals loop of the first-order liquid-hexatic
    transition, where a single homogeneous phase does not exist and the fluid
    sound-speed expression is meaningless (c_s^2 would be negative).  This is how
    the coexistence gap in the figure is *derived* rather than asserted.
    """
    eta = np.asarray(eta, dtype=float)
    return Z_liu_lh(eta) + eta * dZ_liu_lh(eta)


def liu_eos_self_check(verbose: bool = False) -> list[str]:
    """Validate the reconstructed coefficients against published anchors.

    Returns a list of failure messages (empty when everything checks out).
    Anchors, all from Liu (2021):
      * Z_lh(0.720) = Z_solid(0.720) = 10.0335  (tangential branch join)
      * dZ/deta at 0.720 = 40.9 on both branches
      * virial expansion of Z_v gives B2 = 2, B3 ~ 3.128
      * the close term satisfies f(eta_rcp) = 2 at eta_rcp = 0.75
    """
    problems: list[str] = []

    def check(name, got, want, tol):
        ok = abs(got - want) <= tol
        if verbose:
            print(f"  {'ok ' if ok else 'FAIL'} {name}: {got:.6g} (expected {want:.6g} +- {tol:g})")
        if not ok:
            problems.append(f"{name}: got {got:.6g}, expected {want:.6g} +- {tol:g}")

    et = LIU2021_ETA_T
    check("Z_solid(0.720)", float(Z_liu_solid(np.array([et]))[0]), 10.0335, 2e-3)
    check("Z_lh(0.720)", float(Z_liu_lh(np.array([et]))[0]), 10.0335, 5e-3)
    check("dZ_solid(0.720)", float(dZ_liu_solid(np.array([et]))[0]), 40.9, 0.2)
    check("dZ_lh(0.720)", float(dZ_liu_lh(np.array([et]))[0]), 40.9, 0.2)

    # Virial coefficients of Z_v: expand [1 + a2 eta^2 + a3 eta^3 + a4 eta^4]/(1-eta)^2
    coeffs = [1.0, 0.0, 0.125, 1.0 / 18.0, -4.0 / 21.0]
    b = [sum(coeffs[k] * (n - k + 1) for k in range(min(n, 4) + 1)) for n in range(4)]
    check("B2 of Z_v", b[1], 2.0, 1e-9)
    check("B3 of Z_v", b[2], 3.128, 4e-3)

    # Close-term constraint f(eta_rcp) = D = 2 at eta_rcp = 1/c.
    eta_rcp = 1.0 / LIU2021_C
    f_rcp = eta_rcp ** LIU2021_M1 * (LIU2021_B1 + LIU2021_B2 * eta_rcp ** 3)
    check("f(eta_rcp)", f_rcp, 2.0, 1e-2)

    # The coexistence loop must actually be present.
    loop = liu_rigidity(np.linspace(0.702, 0.7175, 40))
    if not np.any(loop < 0.0):
        problems.append("liquid-hexatic coexistence loop (omega<0) not found in 0.702..0.7175")
    elif verbose:
        print(f"  ok  coexistence loop present, min omega = {float(np.min(loop)):.3f}")

    return problems


def speed_of_sound_theory_curves(
    eta: np.ndarray, *, mapping: str, kbt: float, m: float
) -> dict[str, np.ndarray]:
    """Return theory curves without confusing EOS values with sound factors.

    ``mapping='simple'`` reproduces Román et al.'s constant-gamma=2 mapping
    for the historical SPT and Henderson curves.  The modern Kolafa--Rottner
    EOS is always paired with the thermodynamically consistent hard-disk
    expression c_s^2=(kBT/m)[Z+eta Z'+Z^2].  ``mapping='adiabatic'`` uses that
    expression for every EOS.
    """
    eta = np.asarray(eta, dtype=float)
    if mapping == "simple":
        spt = cs_roman_gamma2_from_eos(
            Z_spt_eos(eta), dZ_spt_eos(eta), eta, kbt=kbt, m=m
        )
        henderson = cs_roman_gamma2_from_eos(
            Z_henderson_eos(eta), dZ_henderson_eos(eta), eta, kbt=kbt, m=m
        )
        # The modern EOS is paired with the thermodynamically consistent 2D
        # adiabatic derivative, not with Román's ideal-gas gamma approximation.
        kolafa = cs_adiabatic_2d_monatomic(
            Z_kolafa_rottner_2006(eta), dZ_kolafa_rottner_2006(eta), eta,
            kbt=kbt, m=m,
        )
    elif mapping == "adiabatic":
        spt = cs_adiabatic_2d_monatomic(
            Z_spt_eos(eta), dZ_spt_eos(eta), eta, kbt=kbt, m=m
        )
        henderson = cs_adiabatic_2d_monatomic(
            Z_henderson_eos(eta), dZ_henderson_eos(eta), eta, kbt=kbt, m=m
        )
        kolafa = cs_adiabatic_2d_monatomic(
            Z_kolafa_rottner_2006(eta), dZ_kolafa_rottner_2006(eta), eta,
            kbt=kbt, m=m,
        )
    else:
        raise ValueError(f"Unknown sound-speed mapping: {mapping}")

    # ------------------------------------------------------------------
    # Liu (2021) phase-aware branches.
    #
    # These are always paired with the thermodynamically consistent adiabatic
    # mapping: the whole point of a phase-aware EOS is that Z and dZ/deta carry
    # the transition, so Román's constant-gamma=2 approximation would defeat it.
    #
    # Each branch is emitted only over the eta range where it describes an actual
    # equilibrium phase, and NaN elsewhere, so a plot can simply draw all three
    # and get the correct gaps for free.
    # ------------------------------------------------------------------
    def _masked(values: np.ndarray, keep: np.ndarray) -> np.ndarray:
        out = np.full_like(eta, np.nan, dtype=float)
        np.copyto(out, values, where=keep)
        return out

    liu_all = cs_adiabatic_2d_monatomic(
        Z_liu_global(eta), dZ_liu_global(eta), eta, kbt=kbt, m=m
    )
    # Where the coexistence region is empty, and why.
    #
    # Across a first-order transition there is no single homogeneous phase, so
    # there is no single c_s.  The criterion is *derived*, not a hard-coded
    # window: the rigidity omega = d(Z eta)/deta goes negative inside the loop,
    # which makes c_s^2 negative.  That negative-omega interval is the spinodal
    # region and is genuinely unphysical -- nothing is drawn there.
    #
    # But the flanks of the loop, between the binodal (0.700 / 0.716) and the
    # spinodal (omega = 0), are METASTABLE states that do exist: superheated
    # liquid on the low side, supercooled hexatic on the high side.  A finite,
    # rapidly driven system does not have time or room to phase-separate, so it
    # follows those metastable branches rather than the equilibrium plateau.
    # They are emitted separately so the figure can draw them dashed.
    with np.errstate(invalid="ignore"):
        omega = liu_rigidity(np.minimum(eta, LIU2021_LH_ETA_MAX))
    unstable = omega < 0.0
    # Split the loop at its most unstable point to tell the two flanks apart.
    loop_eta = eta[unstable & (eta < LIU2021_LH_ETA_MAX)]
    if loop_eta.size:
        eta_lo, eta_hi = float(loop_eta.min()), float(loop_eta.max())
    else:
        eta_lo = eta_hi = 0.5 * (ETA_COEX_LO + ETA_COEX_HI)

    stable_fluid = (eta < ETA_COEX_LO) & ~unstable
    meta_liquid = (eta >= ETA_COEX_LO) & (eta <= eta_lo) & ~unstable
    meta_hexatic = (eta >= eta_hi) & (eta < ETA_COEX_HI) & ~unstable
    hexatic = (eta >= ETA_COEX_HI) & (eta < ETA_HEXATIC_SOLID) & ~unstable
    solid = (eta >= ETA_HEXATIC_SOLID) & (eta <= LIU2021_SOLID_ETA_MAX)

    return {
        "ideal": np.full_like(eta, math.sqrt(2.0 * kbt / m), dtype=float),
        "spt": spt,
        "henderson": henderson,
        "kolafa_rottner": kolafa,
        "liu_fluid": _masked(liu_all, stable_fluid),
        "liu_meta_liquid": _masked(liu_all, meta_liquid),
        "liu_meta_hexatic": _masked(liu_all, meta_hexatic),
        "liu_hexatic": _masked(liu_all, hexatic),
        # NB: in the solid this is the EOS *bulk* mode only. A solid has shear
        # rigidity, so the longitudinal mode obeys c_L^2 = c_bulk^2 + mu/rho and
        # this curve is a rigorous LOWER BOUND on the measured piston mode.
        "liu_solid": _masked(liu_all, solid),
        # Equilibrium coexistence: the Maxwell construction fixes eta*Z = const,
        # so (dP/drho)_T = 0 and the equilibrium (fully phase-separated) sound
        # speed vanishes in the thermodynamic limit.  Reported as a scalar so a
        # caller can annotate it; drawing c_s=0 across the band would be honest
        # but visually useless.
        "coexistence_equilibrium_cs": 0.0,
        "coexistence_unstable_eta": (eta_lo, eta_hi),
        # Liu's Z(eta) is a genuinely GLOBAL EOS: it is defined and continuous over
        # the whole range, and the adiabatic expression c_s^2=(kT/m)[Z+eta Z'+Z^2]
        # stays finite through the loop too, because the Z^2 term dominates the
        # negative eta Z'. What fails there is not the arithmetic but the physics:
        # omega = (dP/drho)_T goes negative, so gamma = 1 + Z^2/omega < 0 and Liu's
        # own C_p is negative -- the HOMOGENEOUS state is thermodynamically unstable
        # and a bulk system phase-separates instead.
        #
        # A small box, however, cannot phase-separate (no room for two domains plus
        # an interface), so it stays homogeneous and follows this continuation. That
        # makes the continuous curve the right reference for a 100-disk system, and
        # it is emitted unmasked for exactly that comparison -- drawn faint, because
        # it is not an equilibrium bulk sound speed.
        "liu_homogeneous": liu_all,
    }


ETA_REGION_XMAX = 0.78


def eta_theory_grid(eta_data: np.ndarray) -> np.ndarray:
    eta_max = float(np.nanmax(eta_data)) if len(eta_data) else 0.0
    x_hi = max(ETA_REGION_XMAX, min(0.88, eta_max * 1.08))
    return np.linspace(0.001, x_hi, 500)


def draw_theory_curves(ax, eta_grid, theory, *, theory_label: str) -> None:
    """Draw the full theory-curve set on `ax`, each over its own valid eta range.

    Single implementation shared by every speed-of-sound figure; previously this
    block was copy-pasted across four scripts and had already drifted out of sync.

    Fluid-only EOS (SPT / Henderson / Kolafa-Rottner) stop at the edge of the
    stable isotropic fluid: they still return finite numbers past it, but those
    numbers no longer describe the phase being simulated.  Liu (2021) is drawn as
    three separate branches with a deliberate gap across the first-order
    liquid-hexatic coexistence region.
    """
    eta_grid = np.asarray(eta_grid, dtype=float)
    fluid_mask = eta_grid <= ETA_FLUID_MAX
    kr_mask = eta_grid <= KR2006_PLOT_ETA_MAX

    # Colours are matplotlib's default cycle (C0..C9), matching the earlier
    # figures.  Below eta~0.6 all four fluid EOS lie almost exactly on top of one
    # another, so every theory curve is drawn semi-transparent and the thickest
    # goes down first: that way an overlapped curve still shows through instead of
    # the last-drawn one hiding the rest.
    a = 0.75

    ax.plot(eta_grid[kr_mask], np.asarray(theory["kolafa_rottner"])[kr_mask],
            "-", color="C3", lw=2.6, alpha=a, zorder=2,
            label="Kolafa–Rottner 2006 (thermodynamic; fluid only)")
    ax.plot(eta_grid[fluid_mask], np.asarray(theory["spt"])[fluid_mask],
            "--", color="C1", lw=1.8, alpha=a, zorder=3,
            label=f"SPT ({theory_label}; fluid only)")
    ax.plot(eta_grid[fluid_mask], np.asarray(theory["henderson"])[fluid_mask],
            "-.", color="C2", lw=1.5, alpha=a, zorder=3,
            label=f"Henderson a=0.125 ({theory_label}; fluid only)")

    if "liu_homogeneous" in theory:
        # Faint and underneath: the homogeneous continuation of the global EOS,
        # continuous across the transition. Not an equilibrium bulk sound speed --
        # but it IS the branch a box too small to phase-separate follows.
        ax.plot(eta_grid, theory["liu_homogeneous"], ":", color="C5", lw=1.3,
                alpha=0.55, zorder=2,
                label="Liu 2021 — homogeneous continuation\n(continuous; unstable in bulk)")

    if "liu_fluid" in theory:
        ax.plot(eta_grid, theory["liu_fluid"], "-", color="C5", lw=2.2, alpha=0.9,
                zorder=3, label="Liu 2021 global EOS — stable fluid")
        # Metastable flanks of the first-order loop: superheated liquid and
        # supercooled hexatic.  Dashed because they are metastable, not stable.
        if "liu_meta_liquid" in theory:
            ax.plot(eta_grid, theory["liu_meta_liquid"], "--", color="C5", lw=1.6,
                    alpha=0.9, zorder=3,
                    label="Liu 2021 — metastable (superheated liquid /\nsupercooled hexatic)")
            ax.plot(eta_grid, theory["liu_meta_hexatic"], "--", color="C6", lw=1.6,
                    alpha=0.9, zorder=3)
        ax.plot(eta_grid, theory["liu_hexatic"], "-", color="C6", lw=2.2, alpha=0.9,
                zorder=3, label="Liu 2021 — hexatic")
        ax.plot(eta_grid, theory["liu_solid"], "--", color="C8", lw=2.2, alpha=0.9,
                zorder=3,
                label=r"Liu 2021 — solid (bulk mode; lower bound on $c_L$)")

    ideal_mask = eta_grid <= 0.08
    ax.plot(eta_grid[ideal_mask], np.asarray(theory["ideal"])[ideal_mask],
            ":", color="black", lw=1.5, alpha=0.9, zorder=3,
            label=r"Ideal-gas limit $c_s=\sqrt{2k_BT/m}$")

    # Mark the genuinely empty band: inside the spinodal (omega < 0) c_s^2 would be
    # negative, and at true equilibrium the Maxwell plateau makes the coexistence
    # sound speed vanish.  Saying so beats leaving an unexplained blank.
    span = theory.get("coexistence_unstable_eta")
    if span and np.isfinite(span[0]) and span[1] > span[0]:
        ax.annotate(
            "no single phase\n" + r"($\omega<0$)",
            xy=(0.5 * (span[0] + span[1]), 0.055),
            xycoords=("data", "axes fraction"),
            ha="center", va="bottom", fontsize=6.5, color="#555555",
        )


def eta_regime_index(eta) -> np.ndarray:
    """Classify each eta into 0=stable fluid, 1=coexistence, 2=hexatic, 3=solid."""
    eta = np.asarray(eta, dtype=float)
    idx = np.zeros_like(eta, dtype=int)
    idx[(eta >= ETA_COEX_LO) & (eta < ETA_COEX_HI)] = 1
    idx[(eta >= ETA_COEX_HI) & (eta < ETA_HEXATIC_SOLID)] = 2
    idx[eta >= ETA_HEXATIC_SOLID] = 3
    return idx


def draw_measured_series(ax, eta, cs, cs_err, *, label=r"Simulation fit ($c_s$)",
                         color="C0") -> None:
    """Plot the measured points, connecting them only *within* a phase regime.

    Every point is drawn with its error bar, but consecutive points that straddle a
    phase boundary are not joined: a line across the first-order liquid-hexatic
    transition would assert a continuity that does not exist, which is exactly the
    claim the theory curves are careful not to make.
    """
    eta = np.asarray(eta, dtype=float)
    cs = np.asarray(cs, dtype=float)
    cs_err = np.asarray(cs_err, dtype=float)
    order = np.argsort(eta)
    eta, cs, cs_err = eta[order], cs[order], cs_err[order]

    ax.errorbar(eta, cs, yerr=cs_err, fmt="o", color=color, ms=5,
                capsize=4, zorder=5, label=label, linestyle="none")

    regime = eta_regime_index(eta)
    start = 0
    for i in range(1, len(eta) + 1):
        if i == len(eta) or regime[i] != regime[start]:
            if i - start >= 2:
                ax.plot(eta[start:i], cs[start:i], "-", color=color, lw=1.2, zorder=4)
            start = i


def add_eta_regime_shading(ax, *, x_max: float = ETA_REGION_XMAX,
                           x_min: float = 0.0) -> None:
    """Annotate approximate hard-disk regimes on a packing-fraction axis.

    `x_min` must match the axis lower limit: regions are clipped to it and any
    region left entirely outside is skipped, so a zoomed axis does not scatter
    labels off-canvas (which bbox_inches="tight" would then expand to include).

    ##CHRIS 2026-10-01 SOURCE for the three PHYSICAL boundaries, which until now were unattributed:
    Engel, Anderson, Glotzer, Isobe, Bernard & Krauth, Phys. Rev. E 87, 042134 (2013).
      - liquid-hexatic COEXISTENCE ends at eta ~= 0.716 ("the coexistence phase ends at
        eta ~= 0.716, the region eta >~ 0.716 is thus hexatic"), and the system is uniformly
        liquid at eta = 0.700, with a hexatic bubble at 0.704 and a stripe at 0.708;
      - HEXATIC-SOLID near eta = 0.720 ("the positional order increases drastically at
        eta = 0.720", where C_q0(r) reaches the KTHNY power-law stability limit r^(-1/3));
      - their Mayer-Wood loop has extrema fixed at eta ~= 0.702 and 0.714 independent of N.
    The other boundaries (0.020 / 0.200 / 0.500) are ANALYSIS LABELS, not phase transitions --
    they mark where Z departs from 1 and where the curve stiffens, and nothing depends on them.
    """
    regions = [
        (0.000, 0.020, "ideal-gas\nlimit", "#f4f4f4"),
        (0.020, 0.200, "dilute gas-like\nfluid", "#dceeff"),
        (0.200, 0.500, "dense isotropic\nfluid", "#dff4e3"),
        (0.500, 0.700, "stiff dense\nisotropic fluid", "#fff0c8"),
        (0.700, 0.716, "fluid-hexatic\ncoexist.", "#ffd9b8"),
        (0.716, 0.720, "hexatic", "#f6c7df"),
        (0.720, x_max, "solid", "#ded8f5"),
    ]
    for lo, hi, label, color in regions:
        hi = min(hi, x_max)
        lo = max(lo, x_min)
        if hi <= lo:
            continue
        ax.axvspan(lo, hi, color=color, alpha=0.32, lw=0, zorder=0)
        width = hi - lo
        if width >= 0.030:
            ax.text(
                (lo + hi) / 2.0,
                0.975,
                label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=8,
                color="#333333",
            )
        else:
            narrow_y = 0.73
            if lo >= 0.716:
                narrow_y = 0.38
            elif lo >= 0.700:
                narrow_y = 0.70
            ax.text(
                (lo + hi) / 2.0,
                narrow_y,
                label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="center",
                rotation=90,
                fontsize=6.5,
                color="#333333",
            )
    for x_val in (ETA_COEX_LO, ETA_COEX_HI, ETA_HEXATIC_SOLID):
        if x_min <= x_val <= x_max:
            ax.axvline(x_val, color="#666666", lw=0.7, ls=":", alpha=0.8, zorder=1)


def roman_2002_reference_points(*, radius: float, height: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert the historical "Roman 2002 by L0" values into (eta, cs, cs_err) points.

    The original mapping is keyed by L0 because older pipelines used fixed N_total=100
    (50 per side), fixed r=0.5, fixed height=10, and varied L0 to get different eta.
    Roman's cs values are interpreted as a function of eta, so we place them on the
    corresponding eta axis.
    """
    N_total_ref = 100
    pts = []
    for L0, (cs, err) in ROMAN_2002_BY_L0.items():
        eta = (N_total_ref * math.pi * radius * radius) / max(1e-12, (2.0 * float(L0) * height))
        pts.append((eta, cs, err))
    pts.sort(key=lambda t: t[0])
    eta = np.array([p[0] for p in pts], dtype=float)
    cs = np.array([p[1] for p in pts], dtype=float)
    cs_err = np.array([p[2] for p in pts], dtype=float)
    return eta, cs, cs_err


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=str, default="", help="Input simulation folder (defaults to latest reduced-units run).")
    ap.add_argument("--prefix", type=str, default="edmd_simulation_", help="Folder prefix to auto-detect (default: edmd_simulation_).")
    ap.add_argument(
        "--base",
        type=str,
        default="hspist3/experiments_speed_of_sound",
        help="Base directory for auto-detect (searches recursively for simulation_* folders).",
    )
    ap.add_argument("--drop-frac", type=float, default=0.20, help="Drop this fraction of initial samples after release.")
    ap.add_argument("--f-min", type=float, default=0.01, help="Min frequency for peak search.")
    ap.add_argument("--f-max", type=float, default=0.60, help="Max frequency for peak search (set <=0 to disable).")
    ap.add_argument(
        "--peak-signal",
        type=str,
        default="vel",
        help="Which signal to use for peak-finding: disp|vel|best (default: vel).",
    )
    ap.add_argument("--radius", type=float, default=0.5, help="Particle radius (sigma units).")
    ap.add_argument("--height", type=float, default=10.0, help="Box height (sigma units).")
    ap.add_argument("--wall-thickness", type=float, default=0.05, help="Divider thickness (sigma units).")
    ap.add_argument("--kbt", type=float, default=1.0, help="kB*T in reduced units.")
    ap.add_argument("--m", type=float, default=1.0, help="Particle mass.")
    ap.add_argument(
        "--theory-cs",
        type=str,
        default="adiabatic",
        choices=["adiabatic", "simple"],
        help="Z->c_s mapping: adiabatic (default, thermodynamic c_s^2=(kBT/m)[Z+eta Z'+Z^2]) "
             "or simple (legacy Roman gamma=2, kept for the historical comparison).",
    )
    ap.add_argument("--force-zero-intercept", action="store_true", help="Force nu = c_s * x (no intercept).")
    ap.add_argument(
        "--no-outlier-filter",
        action="store_true",
        help="Disable robust outlier filtering before per-L0 linear fits (not recommended).",
    )
    ap.add_argument(
        "--outlier-mad-z",
        type=float,
        default=6.0,
        help="MAD z-threshold for outlier filtering on implied slopes nu/x (default: 6).",
    )
    ap.add_argument(
        "--outlier-min-keep",
        type=int,
        default=3,
        help="Minimum number of (M) points to keep for each L0 fit after filtering (default: 3).",
    )
    ap.add_argument("--out", type=str, default="", help="Output directory (default: input dir).")
    ap.add_argument("--out-plots", action="store_true", help="Write outputs into a 'plots/' subfolder inside the input dir.")
    ap.add_argument("--tables-only", action="store_true", help="Only write CSV extraction tables; skip all PDF plots.")
    ap.add_argument("--write-final", action="store_true", help="Also write 'FINAL speed_of_sound_on_packing_fracture.pdf' into output dir.")
    ap.add_argument("--roman-ref", action="store_true", help="Overlay Roman et al. (2002) reference points in cs(eta) plot.")
    ap.add_argument(
        "--combined-loglog",
        action="store_true",
        help="Use log-log axes for the combined nu vs K/(2πL_eff) figure (default: linear axes).",
    )
    ap.add_argument("--strict-counts", action="store_true", help="Fail if any run changes Left/Right_Count over time.")
    ap.add_argument(
        "--min-measured-oscillations", type=float, default=0.0,
        help="Frequency-quality gate: require nu*T_used to reach this many cycles (0 disables).",
    )
    ap.add_argument(
        "--max-relative-bin-width", type=float, default=0.0,
        help="Frequency-quality gate: require (df/nu) <= this value (0 disables).",
    )
    ap.add_argument(
        "--min-peak-snr", type=float, default=0.0,
        help="Frequency-quality gate: minimum peak/median-noise power ratio (0 disables).",
    )
    ap.add_argument(
        "--strict-frequency-quality", action="store_true",
        help="Fail closed if any trajectory violates an enabled frequency-quality gate.",
    )
    ap.add_argument(
        "--min-fit-r2", type=float, default=0.0,
        help="Fit-quality gate: minimum weighted R^2 for nu versus K/(2*pi*L_eff) (0 disables).",
    )
    ap.add_argument(
        "--strict-fit-quality", action="store_true",
        help="Fail closed if a per-L0 mass-frequency fit is nonphysical or violates --min-fit-r2.",
    )
    ap.add_argument(
        "--min-group-acceptance", type=float, default=0.0,
        help="Require this accepted-trajectory fraction for every (L0,M) group (0 disables).",
    )
    ap.add_argument(
        "--min-group-repeats", type=int, default=1,
        help="Require at least this many accepted trajectories per (L0,M) group.",
    )
    ap.add_argument(
        "--strict-group-quality", action="store_true",
        help="Fail closed if any (L0,M) group violates its acceptance gates.",
    )
    ap.add_argument(
        "--debug-run",
        type=str,
        default="",
        help="Write a diagnostic wall-trace+PSD plot for one run. "
        "Value can be either a CSV filename (relative to --dir) or 'L0,M,run' "
        "(e.g. '20,200,0' for L0=20, wall mass factor=200, run=0).",
    )
    ap.add_argument(
        "--debug-only",
        action="store_true",
        help="With --debug-run: write the debug plot and exit without aggregating the full sweep.",
    )
    args = ap.parse_args()

    base_dir = Path(args.base)
    in_dir: Optional[Path] = Path(args.dir) if args.dir else None
    if in_dir and not in_dir.exists():
        raise SystemExit(f"Input dir not found: {in_dir}")
    if in_dir is None:
        # Prefer explicit prefix, but our speed_of_sound batches usually use "simulation_*".
        prefixes = [args.prefix, "simulation_"]
        candidates: list[Path] = []
        for pref in prefixes:
            # 1) direct children (legacy layout)
            p = find_latest_sim_dir(base_dir, pref)
            if p is not None:
                candidates.append(p)
            # 2) recursive (new layout: experiments_speed_of_sound/EDMD/.../simulation_*)
            for d in base_dir.rglob(f"{pref}*"):
                if d.is_dir():
                    candidates.append(d)
        if not candidates:
            raise SystemExit(f"No simulation folders found under {base_dir}")
        in_dir = max(candidates, key=lambda p: p.stat().st_mtime)

    if args.out:
        out_dir = Path(args.out)
    elif bool(args.out_plots):
        out_dir = in_dir / "plots"
    else:
        out_dir = in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    f_max = None if args.f_max <= 0 else float(args.f_max)

    # Optional one-run debug plot (before the full aggregation).
    if args.debug_run:
        token = str(args.debug_run).strip()
        debug_csv: Optional[Path] = None
        if token.lower().endswith(".csv"):
            p = Path(token)
            debug_csv = p if p.is_absolute() else (in_dir / p)
        else:
            parts = [p.strip() for p in token.split(",")]
            if len(parts) != 3:
                raise SystemExit("--debug-run must be a CSV filename or 'L0,M,run'")
            L0 = float(parts[0])
            M = int(float(parts[1]))
            run = int(float(parts[2]))
            fname = f"wall_x_positions_L0_{int(round(L0 * 10))}_wallmassfactor_{M}_run{run}.csv"
            debug_csv = in_dir / fname
        if debug_csv is None or not debug_csv.exists():
            raise SystemExit(f"Debug CSV not found: {debug_csv}")
        out_pdf = _debug_psd_and_trace(
            debug_csv,
            out_dir,
            drop_frac=float(args.drop_frac),
            f_min=float(args.f_min),
            f_max=f_max,
        )
        print(f"Wrote:  {out_pdf}")
        if bool(args.debug_only):
            return 0

    # Collect runs
    run_rows = []
    count_warnings = []
    short_warnings = []
    for p in sorted(in_dir.glob("wall_x_positions_L0_*_wallmassfactor_*_run*.csv")):
        key = parse_run_key(p)
        if not key:
            continue
        t, wall_x, disp, left, right = read_wall_csv_5cols(p)
        nu, peak_power, dt, span, df, peak_snr, peak_boundary, raw_bin_nu = dominant_frequency_fft(
            t,
            disp,
            drop_frac=float(args.drop_frac),
            f_min=float(args.f_min),
            f_max=f_max,
            peak_signal=str(args.peak_signal),
        )
        measured_oscillations = nu * span if math.isfinite(nu) and math.isfinite(span) else float("nan")
        relative_bin_width = df / nu if math.isfinite(df) and math.isfinite(nu) and nu > 0 else float("nan")
        metadata = read_wall_csv_metadata(p)
        try:
            l0_exact = float(metadata.get("L0", key.l0))
        except Exception:
            l0_exact = float(key.l0)
        if not math.isfinite(l0_exact) or l0_exact <= 0:
            l0_exact = float(key.l0)
        quality_reasons: list[str] = []
        if not math.isfinite(nu) or nu <= 0:
            quality_reasons.append("no_finite_peak")
        if peak_boundary:
            quality_reasons.append("peak_on_search_boundary")
        if args.min_measured_oscillations > 0 and (
            not math.isfinite(measured_oscillations) or
            measured_oscillations < float(args.min_measured_oscillations)
        ):
            quality_reasons.append("too_few_measured_oscillations")
        if args.max_relative_bin_width > 0 and (
            not math.isfinite(relative_bin_width) or
            relative_bin_width > float(args.max_relative_bin_width)
        ):
            quality_reasons.append("fft_bins_too_coarse")
        if args.min_peak_snr > 0 and (
            not math.isfinite(peak_snr) or peak_snr < float(args.min_peak_snr)
        ):
            quality_reasons.append("peak_snr_too_low")
        frequency_quality_pass = not quality_reasons
        # Heuristic: if df is too coarse, FFT peak becomes quantized and the cs fit collapses.
        if math.isfinite(df) and df > 0.02:
            short_warnings.append(
                f"{p.name}: span={span:.3f} → df≈{df:.4f} (increase --steps and/or --fixed-dt)"
            )
        counts_constant = bool(len(left)) and bool(len(right)) and bool(np.all(left == left[0])) and bool(np.all(right == right[0]))
        if not counts_constant:
            # give a compact hint: first differing sample (if any)
            idx_bad = None
            if len(left) and len(right):
                bad = np.where((left != left[0]) | (right != right[0]))[0]
                if len(bad):
                    idx_bad = int(bad[0])
            msg = f"{p.name}: Left/Right_Count changed over time"
            if idx_bad is not None and idx_bad < len(t):
                msg += f" (first change at t={float(t[idx_bad]):.3f})"
            count_warnings.append(msg)
        # N/side from first valid counts (divider should conserve counts)
        n_total = int(left[0] + right[0]) if len(left) else 0
        n_side = int(round(n_total / 2)) if n_total else 0
        run_rows.append(
            dict(
                file=str(p.name),
                # The legacy filename stores only L0*10 and truncates values
                # such as 5.61.  The CSV metadata is the canonical exact value.
                L0=l0_exact,
                wall_mass_factor=key.wall_mass_factor,
                run=key.run,
                nu=nu,
                peak_power=peak_power,
                peak_snr=peak_snr,
                peak_on_search_boundary=int(peak_boundary),
                raw_bin_nu=raw_bin_nu,
                dt=dt,
                span=span,
                df=df,
                measured_oscillations=measured_oscillations,
                relative_bin_width=relative_bin_width,
                frequency_quality_pass=int(frequency_quality_pass),
                frequency_quality_reasons=";".join(quality_reasons),
                target_oscillations=metadata.get("Target_Oscillations", ""),
                predicted_frequency=metadata.get("Predicted_Frequency", ""),
                planned_steps=metadata.get("Planned_Steps", ""),
                planned_duration=metadata.get("Planned_Duration", ""),
                left=int(left[0]) if len(left) else 0,
                right=int(right[0]) if len(right) else 0,
                N_total=n_total,
                N_side=n_side,
                counts_constant=int(1 if counts_constant else 0),
            )
        )

    if not run_rows:
        raise SystemExit(f"No run CSVs found in {in_dir}")
    if count_warnings:
        if bool(args.strict_counts):
            raise SystemExit("Counts validation failed:\n  - " + "\n  - ".join(count_warnings[:25]))
        print("⚠️ Counts validation warnings (possible tunneling / invalid runs):")
        for w in count_warnings[:25]:
            print(f"  - {w}")
        if len(count_warnings) > 25:
            print(f"  ... ({len(count_warnings)-25} more)")
    if short_warnings:
        print("⚠️ Run-duration warnings (FFT resolution too coarse; cs fit may be meaningless):")
        for w in short_warnings[:15]:
            print(f"  - {w}")
        if len(short_warnings) > 15:
            print(f"  ... ({len(short_warnings)-15} more)")

    quality_failures = [r for r in run_rows if not int(r["frequency_quality_pass"])]
    quality_csv = out_dir / "speed_of_sound_frequency_quality_failures.csv"
    if quality_failures:
        with quality_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(run_rows[0].keys()))
            w.writeheader()
            w.writerows(quality_failures)
        print(f"⚠️ Frequency-quality failures: {len(quality_failures)}/{len(run_rows)} -> {quality_csv}")

    # Write per-run table
    runs_csv = out_dir / "speed_of_sound_runs.csv"
    with runs_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(run_rows[0].keys()))
        w.writeheader()
        w.writerows(run_rows)
    if quality_failures and bool(args.strict_frequency_quality):
        raise SystemExit(
            f"Frequency-quality gate failed for {len(quality_failures)}/{len(run_rows)} trajectories. "
            f"See {quality_csv}; complete per-run diagnostics are in {runs_csv}."
        )

    # Group by (L0, M)
    groups = {}
    for row in run_rows:
        k = (row["L0"], row["wall_mass_factor"])
        groups.setdefault(k, []).append(row)

    group_rows = []
    for (L0, M), rows in sorted(groups.items()):
        accepted_rows = [r for r in rows if int(r["frequency_quality_pass"])]
        nus = np.array([r["nu"] for r in accepted_rows], dtype=float)
        nus = nus[np.isfinite(nus)]
        n = int(len(nus))
        mean_nu = float(np.mean(nus)) if n else float("nan")
        std_nu = float(np.std(nus, ddof=1)) if n >= 2 else float("nan")
        stderr_nu = float(std_nu / math.sqrt(n)) if n >= 2 else float("nan")
        # FFT frequency resolution (df=1/span). Even if the peak bin is stable across repeats,
        # the true frequency is only resolved to ~±df/2. Use this later as a floor for error bars.
        dfs = np.array([r["df"] for r in accepted_rows], dtype=float)
        dfs = dfs[np.isfinite(dfs) & (dfs > 0)]
        df_mean = float(np.mean(dfs)) if len(dfs) else float("nan")
        # counts sanity
        N_side = int(rows[0]["N_side"])
        N_total = int(rows[0]["N_total"])
        acceptance_fraction = n / len(rows) if rows else 0.0
        group_reasons: list[str] = []
        if n < int(args.min_group_repeats):
            group_reasons.append("too_few_accepted_repeats")
        if args.min_group_acceptance > 0 and acceptance_fraction < float(args.min_group_acceptance):
            group_reasons.append("acceptance_fraction_below_threshold")
        group_rows.append(
            dict(
                L0=L0,
                wall_mass_factor=M,
                n=n,
                n_requested=len(rows),
                n_quality_failed=len(rows) - len(accepted_rows),
                acceptance_fraction=acceptance_fraction,
                group_quality_pass=int(not group_reasons),
                group_quality_reasons=";".join(group_reasons),
                nu_mean=mean_nu,
                nu_std=std_nu,
                nu_stderr=stderr_nu,
                df_mean=df_mean,
                N_side=N_side,
                N_total=N_total,
            )
        )

    groups_csv = out_dir / "speed_of_sound_groups.csv"
    with groups_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(group_rows[0].keys()))
        w.writeheader()
        w.writerows(group_rows)
    group_failures = [r for r in group_rows if not int(r["group_quality_pass"])]
    group_failure_csv = out_dir / "speed_of_sound_group_quality_failures.csv"
    if group_failures:
        with group_failure_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(group_rows[0].keys()))
            w.writeheader()
            w.writerows(group_failures)
        print(f"⚠️ Group-quality failures: {len(group_failures)}/{len(group_rows)} -> {group_failure_csv}")
        if bool(args.strict_group_quality):
            raise SystemExit(
                f"Group-quality gate failed for {len(group_failures)}/{len(group_rows)} (L0,M) groups. "
                f"See {group_failure_csv}."
            )

    # Fit per L0
    by_l0 = {}
    for row in group_rows:
        by_l0.setdefault(row["L0"], []).append(row)

    summary_rows = []
    fit_failure_rows = []
    for L0, rows in sorted(by_l0.items()):
        # keep only finite means and stderrs
        rows2 = [
            r for r in rows
            if math.isfinite(r["nu_mean"]) and int(r.get("group_quality_pass", 1))
        ]
        if len(rows2) < 3:
            continue
        N_side = int(rows2[0]["N_side"])
        N_total = int(rows2[0]["N_total"])
        radius = float(args.radius)
        height = float(args.height)
        diameter = 2.0 * radius
        if L0 <= diameter:
            print(f"Skipping L0={L0:g}: L_eff=L0-2r <= 0 (r={radius:g})")
            continue

        # nominal packing fraction (geometric area)
        eta = (N_total * math.pi * radius * radius) / max(1e-12, (2.0 * L0 * height))

        # Match the historical convention used in `wall_x_FFT.py` and older plots:
        # for diameter σ=1 (radius=0.5), L_eff = L0 - 1.
        # We intentionally do NOT correct for wall thickness here to keep comparisons consistent.
        L_eff = L0 - diameter

        xs = []
        ys = []
        sys = []
        for r in rows2:
            M = float(r["wall_mass_factor"])
            alpha = M / max(1e-12, 2.0 * float(N_side))
            K = k_root_bisect(alpha)
            x = K / (2.0 * math.pi * L_eff)
            xs.append(x)
            ys.append(float(r["nu_mean"]))
            # use stderr; if missing, fall back to std or a tiny epsilon
            s = r["nu_stderr"]
            if not math.isfinite(s) or s <= 0:
                s = r["nu_std"]
            # Impose a resolution floor based on FFT bin width: nu is only known to ~±df/2.
            # This prevents the weighted fit collapsing when nu repeats quantize to the same bin.
            df_mean = float(r.get("df_mean", float("nan")))
            if math.isfinite(df_mean) and df_mean > 0:
                s = max(float(s), 0.5 * df_mean)
            if not math.isfinite(s) or s <= 0:
                s = 1e-6
            sys.append(float(s))

        x = np.asarray(xs, dtype=float)
        y = np.asarray(ys, dtype=float)
        sy = np.asarray(sys, dtype=float)
        if bool(args.no_outlier_filter):
            keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sy)
            filt_meta = dict(n_total=int(len(x)), n_used=int(np.sum(keep)), n_dropped=int(len(x) - np.sum(keep)))
        else:
            keep, filt_meta = robust_filter_points(
                x,
                y,
                sy,
                mad_z=float(args.outlier_mad_z),
                min_keep=int(args.outlier_min_keep),
            )
        x_fit = x[keep]
        y_fit = y[keep]
        sy_fit = sy[keep]
        cs, cs_err, b0, b0_err = weighted_linreg(
            x_fit, y_fit, sy_fit, force_zero_intercept=bool(args.force_zero_intercept)
        )

        fit_weights = 1.0 / np.maximum(sy_fit, 1e-12) ** 2
        fit_pred = cs * x_fit + b0
        fit_ybar = (
            float(np.sum(fit_weights * y_fit) / np.sum(fit_weights))
            if len(y_fit) and float(np.sum(fit_weights)) > 0 else float("nan")
        )
        fit_ss_res = float(np.sum(fit_weights * (y_fit - fit_pred) ** 2))
        fit_ss_tot = float(np.sum(fit_weights * (y_fit - fit_ybar) ** 2))
        fit_r2 = 1.0 - fit_ss_res / fit_ss_tot if fit_ss_tot > 0 else float("nan")
        fit_reasons: list[str] = []
        if len(x_fit) < 3:
            fit_reasons.append("too_few_mass_points")
        if not math.isfinite(cs) or cs <= 0:
            fit_reasons.append("nonpositive_sound_speed_slope")
        if not math.isfinite(fit_r2):
            fit_reasons.append("nonfinite_fit_r2")
        elif args.min_fit_r2 > 0 and fit_r2 < float(args.min_fit_r2):
            fit_reasons.append("fit_r2_below_threshold")
        if fit_reasons:
            fit_failure_rows.append(
                dict(
                    L0=L0,
                    eta=eta,
                    cs=cs,
                    cs_err=cs_err,
                    intercept=b0,
                    intercept_err=b0_err,
                    fit_r2=fit_r2,
                    n_points=int(filt_meta.get("n_total", len(x))),
                    n_points_used=int(filt_meta.get("n_used", len(x_fit))),
                    fit_quality_reasons=";".join(fit_reasons),
                )
            )
            print(
                f"⚠️ Rejecting L0={L0:g} fit: {','.join(fit_reasons)} "
                f"(slope={cs:.6g}, weighted R^2={fit_r2:.6g})"
            )
            continue

        kbt = float(args.kbt)
        mval = float(args.m)
        eta_arr = np.array([eta], dtype=float)
        theory_at_eta = speed_of_sound_theory_curves(
            eta_arr, mapping=str(args.theory_cs), kbt=kbt, m=mval
        )
        cs_spt = float(theory_at_eta["spt"][0])
        cs_hend = float(theory_at_eta["henderson"][0])
        cs_kr = (
            float(theory_at_eta["kolafa_rottner"][0])
            if eta <= KR2006_PLOT_ETA_MAX else float("nan")
        )
        rel_spt = (cs / cs_spt - 1.0) if math.isfinite(cs) and math.isfinite(cs_spt) and cs_spt > 0 else float("nan")
        rel_hend = (cs / cs_hend - 1.0) if math.isfinite(cs) and math.isfinite(cs_hend) and cs_hend > 0 else float("nan")
        rel_kr = (cs / cs_kr - 1.0) if math.isfinite(cs) and math.isfinite(cs_kr) and cs_kr > 0 else float("nan")

        summary_rows.append(
            dict(
                L0=L0,
                eta=eta,
                cs=cs,
                cs_err=cs_err,
                intercept=b0,
                intercept_err=b0_err,
                fit_r2=fit_r2,
                n_points=int(filt_meta.get("n_total", len(x))),
                n_points_used=int(filt_meta.get("n_used", len(x_fit))),
                cs_SPT=cs_spt,
                cs_Henderson=cs_hend,
                cs_Kolafa_Rottner_2006=cs_kr,
                rel_err_SPT=rel_spt,
                rel_err_Henderson=rel_hend,
                rel_err_Kolafa_Rottner_2006=rel_kr,
                N_total=N_total,
                N_side=N_side,
                L_eff=L_eff,
            )
        )

    fit_failure_csv = out_dir / "speed_of_sound_fit_quality_failures.csv"
    if fit_failure_rows:
        with fit_failure_csv.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(fit_failure_rows[0].keys()))
            w.writeheader()
            w.writerows(fit_failure_rows)
        print(f"⚠️ Fit-quality failures: {len(fit_failure_rows)} -> {fit_failure_csv}")
        if bool(args.strict_fit_quality):
            raise SystemExit(
                f"Fit-quality gate failed for {len(fit_failure_rows)} L0 fit(s). "
                f"See {fit_failure_csv}."
            )

    if not summary_rows:
        print("No L0 fit entries produced (need >=3 wall masses per L0 and finite peaks); wrote runs/groups only.")
        return 0

    summary_csv = out_dir / "speed_of_sound_summary.csv"
    with summary_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    # Plot cs vs eta + theory. In tables-only mode, do not import Matplotlib;
    # per-eta extraction should remain cheap and avoid repeated font-cache work.
    if not bool(args.tables_only):
        os.environ.setdefault("MPLCONFIGDIR", str(out_dir / ".mplconfig"))
        os.environ.setdefault("MPLBACKEND", "Agg")
        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            raise SystemExit(f"matplotlib required for plotting: {e}")

    eta = np.array([r["eta"] for r in summary_rows], dtype=float)
    cs = np.array([r["cs"] for r in summary_rows], dtype=float)
    cs_err = np.array([r["cs_err"] for r in summary_rows], dtype=float)

    # Sort by eta
    order = np.argsort(eta)
    eta = eta[order]
    cs = cs[order]
    cs_err = cs_err[order]

    kbt = float(args.kbt)
    mval = float(args.m)
    eta_grid = eta_theory_grid(eta)
    theory = speed_of_sound_theory_curves(
        eta_grid, mapping=str(args.theory_cs), kbt=kbt, m=mval
    )
    cs_spt = theory["spt"]
    cs_hend = theory["henderson"]
    cs_kr = theory["kolafa_rottner"]
    theory_label = "Román γ=2" if args.theory_cs == "simple" else "thermodynamic"

    out_pdf = None
    per_l0_pdf = None
    combined_pdf = None
    if not bool(args.tables_only):
        plt.figure(figsize=(8.2, 5.1))
        add_eta_regime_shading(plt.gca())
        draw_theory_curves(plt.gca(), eta_grid, theory, theory_label=theory_label)
        if bool(args.roman_ref):
            eta_r, cs_r, cs_r_err = roman_2002_reference_points(radius=float(args.radius), height=float(args.height))
            plt.errorbar(eta_r, cs_r, yerr=cs_r_err, fmt="^--", color="C4", capsize=3, alpha=0.85, label="Román et al. (2002)")
        # measured last so it sits on top of every theory curve
        draw_measured_series(plt.gca(), eta, cs, cs_err, label=r"$c_s$ (fit)")
        plt.xlabel(r"Packing fraction $\eta$")
        plt.ylabel(r"Speed of sound $c_s$")
        plt.title(r"Speed of sound $c_s$ vs packing fraction $\eta$ (fit + error bars)")
        plt.xlim(0.0, ETA_REGION_XMAX)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        out_pdf = out_dir / "speed_of_sound_cs_vs_eta.pdf"
        plt.savefig(out_pdf, dpi=300)
        plt.close()

    # Plot: nu vs K/(2πL_eff) per L0 (diagnostic) with error bars + fitted slope
    # This is a compact multi-panel "per-L0" report (useful for debugging fit quality).
    if not bool(args.tables_only):
      try:
        l0_list = sorted(by_l0.keys())
        n_panels = len(l0_list)
        if n_panels > 0:
            ncols = 2
            nrows = int(math.ceil(n_panels / ncols))
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12.0, 3.2 * nrows), squeeze=False)
            for idx, L0 in enumerate(l0_list):
                ax = axes[idx // ncols][idx % ncols]
                rows = [r for r in by_l0[L0] if math.isfinite(r["nu_mean"])]
                if len(rows) < 2:
                    ax.set_axis_off()
                    continue

                N_side = int(rows[0]["N_side"])
                N_total = int(rows[0]["N_total"])
                radius = float(args.radius)
                height = float(args.height)
                eta_here = (N_total * math.pi * radius * radius) / max(1e-12, (2.0 * float(L0) * height))
                if float(L0) <= 2.0 * radius:
                    ax.set_axis_off()
                    continue
                L_eff = float(L0) - 2.0 * radius

                xs = []
                ys = []
                yerrs = []
                for r in sorted(rows, key=lambda d: d["wall_mass_factor"]):
                    M = float(r["wall_mass_factor"])
                    alpha = M / max(1e-12, 2.0 * float(N_side))
                    K = k_root_bisect(alpha)
                    x = K / (2.0 * math.pi * L_eff)
                    y = float(r["nu_mean"])
                    s = r["nu_stderr"]
                    if not math.isfinite(s) or s <= 0:
                        s = r["nu_std"]
                    if not math.isfinite(s) or s <= 0:
                        s = 1e-6
                    xs.append(x)
                    ys.append(y)
                    yerrs.append(float(s))
                    ax.annotate(f"M={int(M)}", (x, y), fontsize=8, alpha=0.75)

                x = np.asarray(xs, dtype=float)
                y = np.asarray(ys, dtype=float)
                sy = np.asarray(yerrs, dtype=float)
                if bool(args.no_outlier_filter):
                    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sy)
                else:
                    keep, _ = robust_filter_points(
                        x,
                        y,
                        sy,
                        mad_z=float(args.outlier_mad_z),
                        min_keep=int(args.outlier_min_keep),
                    )
                x_fit = x[keep]
                y_fit = y[keep]
                sy_fit = sy[keep]
                cs_fit, cs_fit_err, b0, _ = weighted_linreg(
                    x_fit, y_fit, sy_fit, force_zero_intercept=bool(args.force_zero_intercept)
                )
                x_line = np.linspace(np.min(x) * 0.95, np.max(x) * 1.05, 200)
                y_line = cs_fit * x_line + (0.0 if bool(args.force_zero_intercept) else b0)

                ax.errorbar(x_fit, y_fit, yerr=sy_fit, fmt="o", capsize=3)
                if np.any(~keep):
                    ax.scatter(
                        x[~keep],
                        y[~keep],
                        marker="x",
                        s=30,
                        linewidths=1.0,
                        alpha=0.25,
                    )
                ax.plot(x_line, y_line, "--", linewidth=1.5, label=fr"$c_s$={cs_fit:.3f}±{cs_fit_err:.3f}")
                if bool(args.combined_loglog):
                    ax.set_xscale("log")
                    ax.set_yscale("log")
                ax.grid(True, linestyle=":", alpha=0.5)
                ax.set_xlabel(r"$K/(2\pi L_{\mathrm{eff}})$")
                ax.set_ylabel(r"Peak $\nu$")
                ax.set_title(fr"$L_0$={float(L0):.1f}, $\eta$={eta_here:.3f}")
                ax.legend(fontsize=9)

            # Hide unused axes
            for j in range(n_panels, nrows * ncols):
                axes[j // ncols][j % ncols].set_axis_off()

            fig.suptitle(r"Speed-of-sound fits: $\nu$ vs $K/(2\pi L_{\mathrm{eff}})$ (per $L_0$)", y=0.995)
            fig.tight_layout()
            per_l0_pdf = out_dir / "speed_of_sound_per_L0_with_eta.pdf"
            fig.savefig(per_l0_pdf, dpi=300)
            plt.close(fig)
        else:
            per_l0_pdf = None
      except Exception:
        per_l0_pdf = None

    # Plot: combined (all L0) frequency vs K-term, similar to wall_x_FFT.py legacy figure.
    if not bool(args.tables_only):
      try:
        import matplotlib.pyplot as plt

        l0_list = sorted(by_l0.keys())
        if l0_list:
            fig, ax = plt.subplots(figsize=(12.0, 6.0))
            colors = plt.cm.viridis(np.linspace(0.0, 1.0, max(1, len(l0_list))))

            for idx, L0 in enumerate(l0_list):
                rows = [r for r in by_l0[L0] if math.isfinite(r["nu_mean"])]
                if len(rows) < 2:
                    continue

                N_side = int(rows[0]["N_side"])
                N_total = int(rows[0]["N_total"])
                radius = float(args.radius)
                height = float(args.height)
                eta_here = (N_total * math.pi * radius * radius) / max(1e-12, (2.0 * float(L0) * height))

                # For legacy comparison use L0-2r (== L0-1 when r=0.5).
                if float(L0) <= 2.0 * radius:
                    continue
                L_eff = float(L0) - 2.0 * radius

                xs = []
                ys = []
                yerrs = []
                for r in sorted(rows, key=lambda d: d["wall_mass_factor"]):
                    M = float(r["wall_mass_factor"])
                    alpha = M / max(1e-12, 2.0 * float(N_side))
                    K = k_root_bisect(alpha)
                    x = K / (2.0 * math.pi * L_eff)
                    y = float(r["nu_mean"])
                    s = r["nu_stderr"]
                    if not math.isfinite(s) or s <= 0:
                        s = r["nu_std"]
                    if not math.isfinite(s) or s <= 0:
                        s = 1e-6
                    xs.append(x)
                    ys.append(y)
                    yerrs.append(float(s))
                    ax.annotate(f"M={int(M)}", (x, y), fontsize=7, alpha=0.7)

                x = np.asarray(xs, dtype=float)
                y = np.asarray(ys, dtype=float)
                sy = np.asarray(yerrs, dtype=float)
                if bool(args.no_outlier_filter):
                    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sy)
                else:
                    keep, _ = robust_filter_points(
                        x,
                        y,
                        sy,
                        mad_z=float(args.outlier_mad_z),
                        min_keep=int(args.outlier_min_keep),
                    )
                x_fit = x[keep]
                y_fit = y[keep]
                sy_fit = sy[keep]
                cs_fit, cs_fit_err, b0, _ = weighted_linreg(
                    x_fit, y_fit, sy_fit, force_zero_intercept=bool(args.force_zero_intercept)
                )
                x_line = np.linspace(np.min(x) * 0.95, np.max(x) * 1.05, 200)
                y_line = cs_fit * x_line + (0.0 if bool(args.force_zero_intercept) else b0)

                col = colors[idx]
                # Plot inliers (used in fit) prominently; outliers as faint x markers.
                ax.errorbar(x_fit, y_fit, yerr=sy_fit, fmt="o", capsize=3, color=col, ecolor=col, alpha=0.9)
                if np.any(~keep):
                    ax.scatter(
                        x[~keep],
                        y[~keep],
                        marker="x",
                        s=28,
                        linewidths=1.0,
                        color=col,
                        alpha=0.25,
                        label=None,
                    )
                ax.plot(
                    x_line,
                    y_line,
                    "--",
                    color=col,
                    linewidth=1.6,
                    label=fr"$L_0$={float(L0):.1f} ($\eta$={eta_here:.3f}), $c_s$={cs_fit:.2f}±{cs_fit_err:.2f}",
                )

            ax.set_xlabel(r"$K/(2\pi L_{\mathrm{eff}})$  (≈ $K/(2\pi(L_0-1))$ for $r=0.5$)")
            ax.set_ylabel(r"Peak frequency $\nu$")
            ax.set_title(r"Fundamental frequency vs theoretical $K$ term (fits per $L_0$)")
            if bool(args.combined_loglog):
                ax.set_xscale("log")
                ax.set_yscale("log")
            ax.grid(True, linestyle=":", alpha=0.5)
            ax.legend(fontsize=8, title="Per-$L_0$ fit", title_fontsize=9)
            fig.tight_layout()
            combined_pdf = out_dir / "speed_of_sound_freq_vs_kterm_combined.pdf"
            fig.savefig(combined_pdf, dpi=300)
            plt.close(fig)
      except Exception:
        combined_pdf = None

    # Optionally write a "FINAL ..." copy for presentations (matching legacy filename).
    final_pdf = None
    if bool(args.write_final) and not bool(args.tables_only):
        final_pdf = out_dir / "FINAL speed_of_sound_on_packing_fracture.pdf"
        # Re-render with a slightly more "final-figure" title.
        # This used to render a reduced curve set (no Kolafa-Rottner, no ideal-gas
        # line) and swallow any failure into final_pdf=None. It now shares the same
        # curve set as every other figure, and a failure is reported rather than
        # silently producing no output.
        import matplotlib.pyplot as plt  # noqa: F401
        plt.figure(figsize=(9.0, 5.6))
        add_eta_regime_shading(plt.gca())
        draw_theory_curves(plt.gca(), eta_grid, theory, theory_label=theory_label)
        if bool(args.roman_ref):
            eta_r, cs_r, cs_r_err = roman_2002_reference_points(radius=float(args.radius), height=float(args.height))
            plt.errorbar(eta_r, cs_r, yerr=cs_r_err, fmt="^--", color="C4", capsize=3, alpha=0.85, label="Román et al. (2002)")
        draw_measured_series(plt.gca(), eta, cs, cs_err)
        plt.xlabel(r"Packing fraction $\eta$")
        plt.ylabel(r"Speed of sound $c_s$")
        plt.title(r"Speed of sound $c_s$ vs packing fraction $\eta$")
        plt.xlim(0.0, ETA_REGION_XMAX)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend(fontsize=7.5, loc="upper left", framealpha=0.92)
        plt.tight_layout()
        plt.savefig(final_pdf, dpi=300)
        plt.close()

    # Print quick accuracy stats (relative to theory lines)
    rel_spt = np.array([r["rel_err_SPT"] for r in summary_rows], dtype=float)[order]
    rel_hend = np.array([r["rel_err_Henderson"] for r in summary_rows], dtype=float)[order]
    if np.any(np.isfinite(rel_spt)):
        mean_abs = float(np.nanmean(np.abs(rel_spt)))
        print(f"Accuracy vs SPT: mean(|Δ|) = {100.0*mean_abs:.2f}%")
    if np.any(np.isfinite(rel_hend)):
        mean_abs = float(np.nanmean(np.abs(rel_hend)))
        print(f"Accuracy vs Henderson(a=0.125): mean(|Δ|) = {100.0*mean_abs:.2f}%")

    print(f"Input:  {in_dir}")
    print(f"Output: {out_dir}")
    print(f"Wrote:  {runs_csv}")
    print(f"Wrote:  {groups_csv}")
    print(f"Wrote:  {summary_csv}")
    if out_pdf is not None:
        print(f"Wrote:  {out_pdf}")
    if per_l0_pdf is not None:
        print(f"Wrote:  {per_l0_pdf}")
    if combined_pdf is not None:
        print(f"Wrote:  {combined_pdf}")
    if final_pdf is not None:
        print(f"Wrote:  {final_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
