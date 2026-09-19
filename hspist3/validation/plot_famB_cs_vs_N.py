#!/usr/bin/env python3
"""##CHRIS 2026-09-12: famB finite-size figure -- c_s vs N at fixed eta and fixed
aspect ratio, one panel per eta, plus the slope d c_s / d(1/sqrt N).

famB geometry: N particles split N/2 per side, so N_side = N/2 and the Roman piston
parameter is alpha = M/(2 N_side) = M/N. L0 and H both scale as sqrt(N), so the
aspect ratio is fixed at each eta and N is the only thing changing.

Per (eta, N):
  accept   strict health contract (forced_advance = clamp_repair = overlap_repair
           = wall_overdue = 0) AND fit quality sigma_nu/nu < 3e-4
  per mass mean nu over surviving repeats, corrected by 1/sqrt(T_i) with T_i read
           from the HD_KE_TRACE audit (drift-first makes this exactly 1)
  c_s      weighted fit nu = c_s x through the origin, x = K/(2 pi L_eff)
  error    1 sigma scatter over the per-mass implied c_s -- this is the plotted bar
Runs partial: any cell with fewer than 3 accepted repeats is skipped, any (eta, N)
with fewer than 3 masses is skipped, and both are reported.
"""
import csv, glob, math, os, re, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fit_nu_damped import fit_trace
import plot_speed_of_sound_edmd as sos

CAMP = sys.argv[1]
OUTBASE = sys.argv[2]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 4
CUT = float(sys.argv[4]) if len(sys.argv) > 4 else 3e-4
R = 0.5
# ##CHRIS: ACCURACY guard, distinct from the precision cut above. The damped-cosine fit
# can lock onto a high-frequency component of the strided trace instead of the piston
# mode; such a fit is very well determined and completely wrong (two traces in famB came
# out 258x the predicted frequency with sigma_nu/nu of 6e-7 and 3e-5, i.e. the tightest
# formal errors in their cells). sigma_nu/nu cannot see this, and neither can the fit rms.
# Measured over the 444 famB traces that pass the 0.03 % cut, nu/nu_predicted spans
# [0.977, 1.237] for every genuine fit; the two aliases sit at 258.4 and 258.7. A window
# of [1/3, 3] is therefore ~2.4x clear of the real data and ~86x clear of the aliases --
# two orders of magnitude wider than the ~1 % effect being measured, so it cannot bias
# c_s. Every rejection is printed, never silent.
NU_RATIO_LO, NU_RATIO_HI = 1.0 / 3.0, 3.0


def scan_cell(cell):
    """Return (eta, N, M, [(nu, T_i), ...], n_traces, n_health_excluded, L0)."""
    log_p = os.path.join(cell, "run.log")
    health, ti_of, pend = {}, {}, None
    if os.path.exists(log_p):
        for ln in open(log_p, errors="replace"):
            k = re.search(r"2b after per-segment equalize\s+N=(\d+)\s+KE_tot=\S+\s+"
                          r"KE_left=(\S+)\s+KE_right=(\S+)", ln)
            if k:
                ns = int(k.group(1)) / 2.0
                pend = (float(k.group(2)) / ns, float(k.group(3)) / ns)
            m = re.search(r"Running: L0 = [\d.]+, M = (\d+)\*m, run = (\d+), seed = (\d+)", ln)
            if m and pend:
                ti_of[(int(m.group(1)), int(m.group(2)))] = pend; pend = None
            h = re.search(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): "
                          r"forced_advance=(\d+) wall_clamp_repairs=(\d+) "
                          r"overlap_repairs=(\d+) wall_overdue=(\d+)", ln)
            if h:
                health[(int(h.group(1)), int(h.group(2)))] = tuple(int(h.group(i)) for i in (4, 5, 6, 7))
    out, n_tr, n_hl, n_cut, L0 = [], 0, 0, 0, None
    aliased = []
    for f in sorted(glob.glob(os.path.join(cell, "wall_x_positions_*.csv"))):
        mm = re.search(r"wallmassfactor_(\d+)_run(\d+)\.csv$", f)
        if not mm: continue
        M, run = int(mm.group(1)), int(mm.group(2))
        n_tr += 1
        if any(health.get((M, run), (0, 0, 0, 0))):
            n_hl += 1; continue
        try:
            with open(f) as fh:
                row0 = next(csv.DictReader(fh), None)
        except Exception:
            continue
        if row0 is None: continue
        if L0 is None:
            try: L0 = float(row0["L0"])
            except Exception: pass
        res = fit_trace(f)
        if not res or not (res["nu"] > 0) or res["sigma_nu"] / res["nu"] >= CUT:
            n_cut += 1; continue
        try: pred = float(row0["Predicted_Frequency"])
        except (KeyError, ValueError, TypeError): pred = 0.0
        if pred > 0:
            ratio = res["nu"] / pred
            if not (NU_RATIO_LO < ratio < NU_RATIO_HI):
                aliased.append((os.path.relpath(f), ratio, res["sigma_nu"] / res["nu"]))
                continue
        ti = ti_of.get((M, run))
        out.append((M, res["nu"], 0.5 * (ti[0] + ti[1]) if ti else float("nan")))
    return (cell, out, n_tr, n_hl, n_cut, L0, aliased)


def main():
    cells = sorted(d for d in glob.glob(f"{CAMP}/eta_*/N*/m_*") if os.path.isdir(d))
    print(f"cells found: {len(cells)}")
    with Pool(NW) as p:
        results = p.map(scan_cell, cells, chunksize=1)

    data, n_tr, n_hl, n_cut, ti_all, alias = {}, 0, 0, 0, [], []
    for cell, out, a, b, c, L0, al in results:
        n_tr += a; n_hl += b; n_cut += c; alias.extend(al)
        m = re.search(r"eta_([\dp]+)/N(\d+)/m_(\d+)$", cell)
        if not m or not out or L0 is None: continue
        eta = float(m.group(1).replace("p", ".")); N = int(m.group(2))
        d = data.setdefault((eta, N), {"L0": L0, "M": {}})
        for M, nu, ti in out:
            d["M"].setdefault(M, []).append((nu, ti))
            if ti == ti: ti_all.append(ti)
    print(f"traces {n_tr}  excluded by health {n_hl}  excluded by the {CUT*100:.3f}% cut {n_cut}"
          f"  excluded as frequency aliases {len(alias)}")
    for f_, r_, q_ in sorted(alias, key=lambda z: -abs(z[1])):
        print(f"  ALIAS nu/nu_pred={r_:.4g}  sigma_nu/nu={q_:.2e}  {f_}")
    if ti_all:
        print(f"T_i over {len(ti_all)} trajectories: min={min(ti_all):.9g} "
              f"max={max(ti_all):.9g} mean={sum(ti_all)/len(ti_all):.9g}")

    pts, skipped = {}, []
    for (eta, N), d in sorted(data.items()):
        L_eff = d["L0"] - 2 * R
        xs, ys, ws, cs_mass, nrun = [], [], [], [], 0
        for M, pairs in sorted(d["M"].items()):
            if len(pairs) < 3:
                skipped.append(f"eta={eta} N={N} M={M}: {len(pairs)} accepted repeats"); continue
            nu = np.array([q[0] for q in pairs]); ti = np.array([q[1] for q in pairs])
            corr = nu / np.sqrt(np.where(np.isfinite(ti) & (ti > 0), ti, 1.0))
            K = sos.k_root_bisect(M / float(N))      # alpha = M/(2 N_side), N_side = N/2
            x = K / (2 * math.pi * L_eff); y = float(np.mean(corr))
            xs.append(x); ys.append(y)
            ws.append(max(float(np.std(corr, ddof=1)) / math.sqrt(len(corr)), 1e-12))
            cs_mass.append(y / x); nrun += len(pairs)
        if len(xs) < 3:
            skipped.append(f"eta={eta} N={N}: only {len(xs)} masses -- no point"); continue
        cs, cs_err = sos.weighted_linreg(np.array(xs), np.array(ys), np.array(ws),
                                         force_zero_intercept=True)[:2]
        pts.setdefault(eta, []).append(dict(N=N, c_s=cs, err_fit=cs_err,
                                            scatter=float(np.std(np.array(cs_mass), ddof=1)),
                                            n_masses=len(xs), n_runs=nrun, L0=d["L0"]))
    for s_ in skipped: print("  skip:", s_)
    if not pts:
        print("nothing complete enough to plot"); return

    etas = sorted(pts)
    ncol = len(etas)
    fig, axes = plt.subplots(1, ncol, figsize=(3.5 * ncol, 4.0), squeeze=False)
    table = []
    for j, eta in enumerate(etas):
        ax = axes[0][j]
        rows = sorted(pts[eta], key=lambda r: r["N"])
        Ns = np.array([r["N"] for r in rows], float)
        cs = np.array([r["c_s"] for r in rows])
        sc = np.array([r["scatter"] for r in rows])
        a_ = np.array([eta])
        Zk = sos.Z_kolafa_rottner_2006(a_); h = 1e-5
        dZk = (sos.Z_kolafa_rottner_2006(a_ + h) - sos.Z_kolafa_rottner_2006(a_ - h)) / (2 * h)
        kr = float(sos.cs_adiabatic_2d_monatomic(Zk, dZk, a_, kbt=1, m=1)[0]) if eta <= 0.69 else float("nan")
        ax.errorbar(Ns, cs, yerr=sc, fmt="o", ms=5, capsize=3, color="#1f4e79",
                    label="famB (bar = 1σ over masses)")
        if kr == kr:
            ax.axhline(kr, color="red", lw=1.4, label=f"KR 2006 = {kr:.4f}")
        ax.set_xscale("log")
        ax.set_xticks(Ns); ax.set_xticks([], minor=True)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlim(Ns.min() / 1.6, Ns.max() * 1.6)
        ax.set_xlabel("N"); ax.set_title(f"η = {eta:.2f}", fontsize=10)
        if j == 0: ax.set_ylabel(r"$c_s$  ($\sigma/t$, $k_BT=m=\sigma=1$)")
        ax.grid(alpha=0.3, which="major")
        ax.legend(fontsize=6.5, loc="best")
        nm = ", ".join(str(r["n_masses"]) for r in rows)
        ax.text(0.02, 0.03, f"masses in fit: {nm}", transform=ax.transAxes, fontsize=6,
                color="0.35", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        # slope of c_s against 1/sqrt(N), weighted by the mass scatter
        if len(rows) >= 3:
            u = 1.0 / np.sqrt(Ns)
            sl, sl_err, ic, ic_err = sos.weighted_linreg(u, cs, np.maximum(sc, 1e-12),
                                                         force_zero_intercept=False)[:4]
            table.append((eta, sl, sl_err, ic, ic_err, len(rows)))
        else:
            table.append((eta, float("nan"), float("nan"), float("nan"), float("nan"), len(rows)))
    fig.suptitle("famB: speed of sound vs system size at fixed η and fixed aspect ratio\n"
                 "drift-first seeding · strict health contract · 0.03 % fit-quality cut · "
                 "bar = 1σ scatter over piston masses",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0.02, 1, 0.90])
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUTBASE}.{ext}", dpi=170)
    print(f"\nwrote {OUTBASE}.pdf / .png")

    print(f"\n{'eta':>6} {'d c_s / d(1/sqrt N)':>21} {'c_s(N->inf)':>20} {'nN':>3}")
    for eta, sl, se, ic, ie, nN in table:
        print(f"{eta:>6.2f} {sl:>12.4f} ± {se:<7.4f} {ic:>11.4f} ± {ie:<7.4f} {nN:>3}")
    with open(f"{OUTBASE}_slopes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["eta", "slope_dcs_dinvsqrtN", "slope_err", "c_s_infinite_N", "c_s_inf_err", "n_N"])
        for r in table: w.writerow([f"{r[0]:.2f}"] + [f"{v:.6g}" for v in r[1:5]] + [r[5]])
    with open(f"{OUTBASE}_points.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["eta", "N", "L0", "c_s", "c_s_err_fit", "c_s_scatter_mass", "n_masses", "n_runs"])
        for eta in etas:
            for r in sorted(pts[eta], key=lambda z: z["N"]):
                w.writerow([f"{eta:.2f}", r["N"], f"{r['L0']:.6f}", f"{r['c_s']:.5f}",
                            f"{r['err_fit']:.5f}", f"{r['scatter']:.5f}", r["n_masses"], r["n_runs"]])
    print(f"wrote {OUTBASE}_slopes.csv and {OUTBASE}_points.csv")


if __name__ == "__main__":
    main()
