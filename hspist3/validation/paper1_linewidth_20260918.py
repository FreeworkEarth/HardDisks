#!/usr/bin/env python3
"""##CHRIS 2026-09-18: item 2 -- Paper 1's linewidth against Mansour's piston Q. Analysis only.

Malek Mansour, Garcia & Baras 2006 linearise the adiabatic piston between two compartments
(Eqs. 33-39): a damped oscillator with, in their scaled time, omega_0 = 4 and damping beta = 2 mu,

    mu = Gamma L_y sqrt( 2 / (Mhat N) ),    Mhat = M + m N / 3,    Gamma = zeta_bulk + eta_shear,

so Q = omega_0 / (2 beta) = 1 / mu and the fractional linewidth of the resonance is Delta f / f = 1/Q.

The A1 v2 divider IS that piston, so Paper 1's spectra measure Gamma for free -- and Gamma is what
sets Paper 2's hydrodynamic friction zeta_hyd = L_y Gamma / X_p. Two things are tested here:

  (a) the SCALING: 1/Q depends on mass only through Mhat^(-1/2), so Delta f/f * sqrt(Mhat) must be
      one constant across the nine masses. That test needs no value of Gamma at all.
  (b) the VALUE: the constant gives Gamma = (Delta f/f) sqrt(Mhat N / 2) / L_y, to be compared with
      the audit's Enskog Gamma = 0.331 at eta = 0.10.

The linewidth is the FWHM of the seed-averaged position spectrum, interpolated between bins; at 200
predicted periods one bin is 0.5 % of f, so a 1.5 % line spans three bins and the interpolation
matters. Trajectories are the A1 v2 records, same health contract as everywhere else.
"""
import json, math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tests_20260913 as T

TD, X_EDGE, LY, NPART = 200, 2.5, 10.0, 100
GAMMA_AUDIT = 0.331          # Enskog zeta_bulk + eta_shear at eta = 0.100 (audit Eq. 5)


def fwhm(P, df, k_lo):
    """FWHM of the tallest peak at or above bin k_lo, linearly interpolated on the half-max level."""
    k = k_lo + int(np.argmax(P[k_lo:]))
    pk = P[k]
    half = 0.5 * pk
    i = k
    while i > k_lo and P[i] > half:
        i -= 1
    if P[i] > half:
        return np.nan, k * df
    lo = (i + (half - P[i]) / (P[i + 1] - P[i])) * df
    j = k
    while j < len(P) - 1 and P[j] > half:
        j += 1
    if P[j] > half:
        return np.nan, k * df
    hi = (j - (half - P[j]) / (P[j - 1] - P[j])) * df
    return hi - lo, k * df


def cell(task):
    eta, L0, M, runs = task
    spec, df0, n = None, None, 0
    for r, p, disc in runs:
        if disc:
            continue
        t, x, nup = T._load(p)
        dt = (t[-1] - t[0]) / (len(t) - 1)
        nn = T._prefix(t, nup, TD)
        P, df = T._spectrum(x[:nn], dt)
        if spec is None:
            spec, df0 = np.zeros_like(P), df
        m = min(len(spec), len(P))
        spec[:m] += P[:m]; n += 1
    if n == 0:
        return None
    spec /= n
    k_lo = max(1, int(round(TD / X_EDGE)))
    w, f0 = fwhm(spec, df0, k_lo)
    # ##CHRIS 2026-09-18: take the bin out in quadrature before comparing with 1/Q. Binning can only
    # BROADEN a line, so the raw FWHM is an upper bound; one bin is 1/TD of f by construction.
    rel_raw = float(w / f0) if np.isfinite(w) else float("nan")
    rel_bin = 1.0 / TD
    rel = float(math.sqrt(max(0.0, rel_raw ** 2 - rel_bin ** 2))) if np.isfinite(rel_raw) else float("nan")
    return dict(eta=eta, M=M, n=n, f0=float(f0), fwhm=float(w), rel_raw=rel_raw, rel=rel)


def main():
    table = T.a1_leaf_table()
    want = [min((l for l in table), key=lambda l: abs(l["eta"] - x)) for x in (0.10, 0.30, 0.50)]
    from multiprocessing import Pool
    tasks = []
    for leaf in want:
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            if runs:
                tasks.append((leaf["eta"], float(leaf["L0"]), M, runs))
    with Pool(9) as pool:
        cells = [c for c in pool.map(cell, tasks, chunksize=1) if c]
    out = []
    for leaf in want:
        cs = sorted([c for c in cells if c["eta"] == leaf["eta"]], key=lambda c: c["M"])
        if not cs:
            continue
        print(f"\n### η = {leaf['eta']:.6f}  (bin width = {100.0 / TD:.2f} % of f, 200 periods, 25 seeds)\n")
        print("| M | α = M/N | M̂ | raw Δf/f [%] | bin-deconvolved [%] | 1/Q with Γ = 0.331 [%] | ratio | Γ implied |")
        print("|---|---|---|---|---|---|---|---|")
        rows = []
        for c in cs:
            Mh = c["M"] + NPART / 3.0
            mu = GAMMA_AUDIT * LY * math.sqrt(2.0 / (Mh * NPART))
            gimp = c["rel"] * math.sqrt(Mh * NPART / 2.0) / LY
            rows.append((c["M"], c["M"] / NPART, Mh, c["rel"], mu, gimp))
            print(f"| {c['M']} | {c['M'] / NPART:g} | {Mh:.1f} | {100 * c['rel_raw']:.2f} | {100 * c['rel']:.2f} | "
                  f"{100 * mu:.2f} | {c['rel'] / mu:.2f} | {gimp:.3f} |")
        g = np.array([r[5] for r in rows]); a = np.array([r[1] for r in rows])
        k = a > 2
        print(f"\nΓ implied, all nine masses: {g.mean():.3f} ± {g.std(ddof=1) / math.sqrt(len(g)):.3f} "
              f"(spread {100 * g.std(ddof=1) / g.mean():.1f} %)")
        if k.sum() >= 3:
            print(f"Γ implied, α > 2 only (where K ≪ 1 holds): {g[k].mean():.3f} ± "
                  f"{g[k].std(ddof=1) / math.sqrt(k.sum()):.3f} (spread {100 * g[k].std(ddof=1) / g[k].mean():.1f} %)")
        out.extend(rows)
    json.dump([{"M": r[0], "alpha": r[1], "Mhat": r[2], "rel": r[3], "mu": r[4], "gamma": r[5]} for r in out],
              open(os.path.join(T.PLOTS, "260918_A1v2_linewidth.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
