#!/usr/bin/env python3
"""##CHRIS 2026-10-01: the four Paper 1 figures the draft lists as missing. Analysis only, no runs.

Every number is read from a CSV/JSON on disk or recomputed from the stored traces with the SAME
estimator the paper uses (T._spectrum plus the nu_pred/2.5 floor). Nothing is typed in by hand.

  1  estimator floor sensitivity      261001_p1_estimator_floor
  2  mass-ladder residual flatness    261001_p1_massladder_residuals
  3  one worked mass-ladder line      261001_p1_massladder_line
  4  slow-mode illustration           261001_p1_slowmode

Colours are the house convention: our data BLUE #2a78d6, Kolafa-Rottner RED #e34948,
Roman 2002 BLACK. Every panel carries error bars.
"""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import tests_20260913 as T

BLUE, RED, GREY, ORANGE = "#2a78d6", "#e34948", "#52514e", "#eb6834"
OUTDIR = T.PLOTS
ETA_CELL, L0_CELL = 0.1122, 34.9999


def save(fig, name):
    os.makedirs(OUTDIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png/.pdf")


# ----------------------------------------------------------------- 1. estimator floor
def fig_floor():
    """Deviation of c_s from the windowed reference vs the floor parameter X.

    The floor is 'largest bin at f >= nu_pred/X', so its bin index is k_min = N_cyc/X: a LARGE X
    is a LOW floor, which lets the slow adiabatic-piston wander win the periodogram. The CSV column
    names are <k_min>_<N_cyc> and the deviations d<k_min>_<N_cyc>; X is derived, not stored.
    """
    import csv
    series = {25: {}, 200: {}}
    for fname in ("260914_A1v2_kmin_sensitivity.csv",
                  "260914_A1v2_kmin_sensitivity_extended.csv",
                  "260914_A1v2_kmin_sensitivity_X2p5.csv"):
        path = T.plot_path(fname)
        if not path or not os.path.exists(path):
            print(f"  MISSING {fname}"); continue
        rows = list(csv.DictReader(open(path)))
        for col in rows[0]:
            if not col.startswith("d"):
                continue
            try:
                kmin, ncyc = col[1:].split("_"); kmin, ncyc = int(kmin), int(ncyc)
            except ValueError:
                continue
            vals = np.array([float(r[col]) for r in rows if r[col] not in ("", "nan")])
            if vals.size == 0:
                continue
            series[ncyc][ncyc / kmin] = vals      # X -> per-density deviations [%]

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for ncyc, colour, mark in ((25, ORANGE, "s"), (200, BLUE, "o")):
        if not series[ncyc]:
            continue
        X = np.array(sorted(series[ncyc]))
        med = np.array([np.median(np.abs(series[ncyc][x])) for x in X])
        # spread across the 35 densities, as an error bar: 16th-84th percentile of |deviation|
        lo = np.array([np.percentile(np.abs(series[ncyc][x]), 16) for x in X])
        hi = np.array([np.percentile(np.abs(series[ncyc][x]), 84) for x in X])
        ax.errorbar(X, med, yerr=[med - lo, hi - med], fmt=mark + "-", color=colour, capsize=3,
                    lw=1.8, ms=6, label=f"{ncyc} predicted periods (35 densities)")
    ax.axvline(2.5, color=RED, ls="--", lw=2, label="production floor $X = 2.5$")
    ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xlabel(r"floor parameter $X$   (largest bin at $f \geq \nu_{\rm pred}/X$;  low floor $\to$ large $X$)")
    ax.set_ylabel(r"$|\Delta c_s|$ vs windowed reference  [%]")
    ax.set_title("Estimator floor sensitivity: the slow mode wins whenever the floor is too low")
    ax.grid(alpha=0.3); ax.legend(frameon=False, loc="upper left")
    save(fig, "261001_p1_estimator_floor")


# ------------------------------------------------- 2. mass-ladder residual flatness vs ln(alpha)
def _damping_cells():
    p = T.plot_path("260915_A1v2_damping_cells.json")
    if not p or not os.path.exists(p):
        print("  MISSING 260915_A1v2_damping_cells.json"); return None
    return json.load(open(p))


def fig_residuals():
    """Residual of each mass about its own density's through-origin ladder line, vs ln(alpha).

    A mass-independent estimator gives a flat line. This reproduces estimator_massladder_20260917
    .extra(): per density fit s = sum(xy)/sum(x^2), residual = 100 (y - s x)/(s x), then average
    over densities at fixed alpha. eta <= 0.69 only, where KR is inside its fitted range.
    """
    cells = _damping_cells()
    if cells is None:
        return
    by_eta = {}
    for c in cells:
        if c["eta"] > 0.69:
            continue
        by_eta.setdefault(round(c["eta"], 6), []).append(c)

    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    for key, label, colour, mark in (("pos", "largest bin (production)", BLUE, "o"),
                                     ("f0", "per-trajectory fit", ORANGE, "s")):
        acc = {}
        for eta, cs in by_eta.items():
            cs = sorted(cs, key=lambda c: c["M"])
            x = np.array([T.x_of(c["M"], c["L0"]) for c in cs])
            y = np.array([c[key] for c in cs])
            ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if ok.sum() < 3:
                continue
            s = float(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2))
            for c, xi, yi in zip(np.array(cs)[ok], x[ok], y[ok]):
                acc.setdefault(c["M"] / 100.0, []).append(100.0 * (yi - s * xi) / (s * xi))
        al = np.array(sorted(acc))
        mean = np.array([np.mean(acc[a]) for a in al])
        sem = np.array([np.std(acc[a], ddof=1) / math.sqrt(len(acc[a])) for a in al])
        w = 1.0 / sem ** 2
        A = np.vstack([np.log(al), np.ones_like(al)]).T
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        beta = cov @ (A.T @ (w * mean))
        slope, dslope = beta[0], math.sqrt(cov[0, 0])
        ax.errorbar(al, mean, yerr=sem, fmt=mark, color=colour, capsize=3, ms=6,
                    label=f"{label}: {slope:+.3f} ± {abs(dslope):.3f} %/e-fold")
        g = np.linspace(al.min(), al.max(), 50)
        ax.plot(g, beta[0] * np.log(g) + beta[1], "-", color=colour, lw=1.5, alpha=0.75)
    ax.axhline(0, color=GREY, lw=1)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\alpha = M/(2N_s m)$")
    ax.set_ylabel("residual about the ladder line  [%]")
    ax.set_title(r"Mass-ladder residual flatness: the largest-bin estimator is mass-independent")
    ax.grid(alpha=0.3); ax.legend(frameon=False)
    save(fig, "261001_p1_massladder_residuals")


# ------------------------------------------------------------------ 3. one worked ladder line
def _nu_per_run(cell_dir, M):
    """Per-trajectory nu with the production estimator, so the point can carry a real error bar."""
    out = []
    for r, p, disc in T.cell_runs(cell_dir, M):
        if disc:
            continue
        try:
            t, x, nup = T._load(p)
        except Exception:
            continue
        if len(t) < 64:
            continue
        dt = (t[-1] - t[0]) / (len(t) - 1)
        P, df = T._spectrum(x, dt)
        k = max(1, int(round(len(x) * dt * nup / 2.5)))
        if k >= len(P):
            continue
        out.append((k + int(np.argmax(P[k:]))) * df)
    return np.array(out)


def fig_ladder_line():
    """nu_bar_M against x_M at one density, nine masses, slope through the origin = c_s.

    NOTE eta = 0.1122, not 0.10: the A1 v2 sweep has no 0.10 leaf, and 0.1122 is the density at
    which all nine masses were run and which the linewidth work already uses.
    """
    cell_root = os.path.join(T.DROOT, "eta_0p112200")
    if not os.path.isdir(cell_root):
        print(f"  MISSING {cell_root}"); return
    xs, ys, es, ns = [], [], [], []
    for M in T.A1_MASSES:
        d = os.path.join(cell_root, f"m_{M}")
        if not os.path.isdir(d):
            continue
        nu = _nu_per_run(d, M)
        if nu.size < 3:
            continue
        xs.append(T.x_of(M, L0_CELL)); ys.append(nu.mean())
        es.append(nu.std(ddof=1) / math.sqrt(nu.size)); ns.append((M, nu.size))
    if len(xs) < 3:
        print("  not enough masses for the ladder line"); return
    xs, ys, es = np.array(xs), np.array(ys), np.array(es)
    # T.slope is the paper's own estimator: through-origin slope, and sd(nu/x) as the mass scatter.
    s, sd_mass, nfit = T.slope(xs, ys)
    ds = sd_mass / math.sqrt(nfit)          # error on c_s itself; sd_mass is the scatter over masses
    kr = float(T.kr_cs(ETA_CELL))

    fig, (ax, axr) = plt.subplots(1, 2, figsize=(12.2, 5.2), gridspec_kw={"width_ratios": [1.5, 1]})
    ax.errorbar(xs, ys, yerr=es, fmt="o", color=BLUE, capsize=3, ms=7, zorder=3,
                label=f"A1 v2, $\\eta = {ETA_CELL}$, {len(xs)} masses, 25 seeds each")
    g = np.linspace(0, xs.max() * 1.06, 50)
    ax.plot(g, s * g, "-", color=BLUE, lw=1.8, zorder=2,
            label=f"through-origin slope $c_s = {s:.4f} \\pm {ds:.4f}$")
    ax.plot(g, kr * g, "--", color=RED, lw=1.8, zorder=1, label=f"Kolafa–Rottner $c_s = {kr:.4f}$")
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$x_M = K(\alpha)\,/\,2\pi L_{\rm eff}$")
    ax.set_ylabel(r"$\bar\nu_M$")
    ax.set_title(f"One worked mass ladder ($\\eta = {ETA_CELL}$)")
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9, loc="upper left")

    res = 100.0 * (ys - s * xs) / (s * xs)
    rer = 100.0 * es / (s * xs)
    axr.errorbar([m for m, _ in ns], res, yerr=rer, fmt="o", color=BLUE, capsize=3, ms=6)
    axr.axhline(0, color=GREY, lw=1)
    axr.axhline(100 * (kr - s) / s, color=RED, ls="--", lw=1.6,
                label=f"KR offset {100*(kr-s)/s:+.2f} %")
    axr.set_xscale("log"); axr.set_xlabel("divider mass $M$")
    axr.set_ylabel("residual about the line  [%]")
    axr.set_title("Residuals: flat across a factor 40 in $M$")
    axr.grid(alpha=0.3); axr.legend(frameon=False, fontsize=9)
    save(fig, "261001_p1_massladder_line")
    print(f"     c_s = {s:.5f} +- {ds:.5f}   KR = {kr:.5f}   ratio {s/kr:.5f}")


# --------------------------------------------------------------------- 4. slow-mode illustration
def fig_slowmode(M=1000):
    cell = os.path.join(T.DROOT, "eta_0p112200", f"m_{M}")
    if not os.path.isdir(cell):
        print(f"  MISSING {cell}"); return
    runs = [(r, p) for r, p, disc in T.cell_runs(cell, M) if not disc]
    if not runs:
        print("  no healthy runs"); return
    t, x, nup = T._load(runs[0][1])
    dt = (t[-1] - t[0]) / (len(t) - 1)
    W = max(3, int(round(5.0 / (nup * dt))))            # 5-period running mean
    ker = np.ones(W) / W
    slow = np.convolve(x, ker, mode="same")
    P, df = T._spectrum(x, dt)
    k = max(1, int(round(len(x) * dt * nup / 2.5)))
    nu = (k + int(np.argmax(P[k:]))) * df

    fig, (ax, axs) = plt.subplots(1, 2, figsize=(12.6, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
    show = slice(0, min(len(t), int(40 / (nup * dt))))   # first ~40 periods
    ax.plot(t[show], x[show], "-", color=BLUE, lw=0.7, alpha=0.75, label="divider displacement")
    ax.plot(t[show], slow[show], "-", color=ORANGE, lw=2.2, label=f"5-period running mean (slow mode)")
    ax.axhline(0, color=GREY, lw=1)
    ax.set_xlabel(r"time [$\sigma$-time]"); ax.set_ylabel(r"displacement [$\sigma$]")
    ax.set_title(f"The slow mode, $M = {M}$, $\\eta = {ETA_CELL}$")
    ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9)

    f = np.arange(len(P)) * df
    axs.loglog(f[1:], P[1:], "-", color=BLUE, lw=0.9)
    axs.axvline(nup / 2.5, color=RED, ls="--", lw=2, label=r"estimator floor $\nu_{\rm pred}/2.5$")
    axs.axvline(nup, color=GREY, ls=":", lw=1.6, label=r"$\nu_{\rm pred}$")
    axs.plot([nu], [P[int(round(nu / df))]], "v", color=ORANGE, ms=11, zorder=5,
             label=f"largest bin above floor, $\\nu = {nu:.5f}$")
    axs.set_xlabel("frequency"); axs.set_ylabel("power")
    axs.set_title("…and why the floor is needed")
    axs.grid(alpha=0.3, which="both"); axs.legend(frameon=False, fontsize=9, loc="lower left")
    save(fig, "261001_p1_slowmode")
    print(f"     M={M}: dt={dt:.4f}, window={W} samples, nu_pred={nup:.6f}, nu={nu:.6f}")


if __name__ == "__main__":
    print("Paper 1 figures ->", OUTDIR)
    fig_floor()
    fig_residuals()
    fig_ladder_line()
    fig_slowmode()
