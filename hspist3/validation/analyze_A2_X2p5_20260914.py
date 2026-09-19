#!/usr/bin/env python3
"""##CHRIS 2026-09-14: A2 (fixed-eta N ladder) re-analysed with the Paper 1 primary estimator.

Estimator, per trajectory: mean-subtracted FFT of the whole post-release divider trace (no transient
drop), largest bin at f >= nu_pred / 2.5, i.e. bins below k_min = round(N_cyc / 2.5) excluded, where
N_cyc = T * nu_pred is the record length in predicted periods. One-sided lower edge, no upper edge.
Mean and SD over health-clean seeds per (eta, N, M). c_s per (eta, N) = through-origin slope of the
per-mass mean nu against x_M = K(M/N) / (2 pi (L0 - 2r)), cot K = (M/N) K; error = 1 sigma scatter of
the per-mass nu_M / x_M. Finite size: c_s = a + b / sqrt(N) per eta, weighted by that scatter.
Also carries the w = 3 window value per trajectory as a cross-check.

usage: analyze_A2_X2p5_20260914.py ROOT1,ROOT2,... OUT_PREFIX
"""
import csv, glob, math, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

X_EDGE = 2.5
R = 0.5


def x_of(M, N, L0):
    return T.k_root(M / float(N)) / (2 * math.pi * (L0 - 2 * R))


def cell(task):
    root, cdir = task
    parts = cdir.rstrip("/").split("/")
    eta = float(parts[-3].split("_")[1].replace("p", ".")); N = int(parts[-2][1:]); M = int(parts[-1][2:])
    out = dict(root=os.path.basename(root), eta=eta, N=N, M=M, nu=[], nu_w3=[], ncyc=[], disc=0, L0=None)
    runs = list(T.cell_runs(cdir, M))
    # ##CHRIS 2026-09-16: campaigns run one exact-seed invocation per folder (m_M/r<k>/, trace named run0,
    # health line in that folder's stdout.log) -- e.g. A2_dilute_20260916. Same trace check and the same
    # rule for the health contract; a folder layout and a flat layout never coexist in one cell.
    for k, rd in enumerate(sorted(glob.glob(os.path.join(cdir, "r*")))):
        g = glob.glob(os.path.join(rd, f"wall_x_positions_L0_*_wallmassfactor_{M}_run0.csv"))
        if not (g and T.trace_check(g[0])[0]):
            continue
        log = os.path.join(rd, "stdout.log")
        disc = os.path.exists(log) and "EDMD-HEALTH" in open(log, errors="replace").read()
        runs.append((k, g[0], disc))
    for r, p, disc in runs:
        if disc:
            out["disc"] += 1
            continue
        import pandas as pd
        d = pd.read_csv(p, usecols=["Time", "Displacement(σ)", "Predicted_Frequency", "L0"])
        t = d["Time"].to_numpy(float); x = d["Displacement(σ)"].to_numpy(float); nup = float(d["Predicted_Frequency"].iloc[0])
        out["L0"] = float(d["L0"].iloc[0])
        dt = (t[-1] - t[0]) / (len(t) - 1)
        P, df = T._spectrum(x, dt)
        ncyc = len(x) * dt * nup; kmin = max(1, int(round(ncyc / X_EDGE)))
        out["nu"].append((kmin + int(np.argmax(P[kmin:]))) * df); out["ncyc"].append(ncyc)
        f = np.arange(len(P)) * df
        sel = np.where((f >= nup / 3) & (f <= 3 * nup) & (np.arange(len(P)) >= 1))[0]
        out["nu_w3"].append(float(f[sel[int(np.argmax(P[sel]))]]))
    return out


def fit_point(cells, key):
    xs, ys, rows = [], [], []
    for c in sorted(cells, key=lambda c: c["M"]):
        if len(c[key]) < 3:
            continue
        a = np.array(c[key]); x = x_of(c["M"], c["N"], c["L0"])
        xs.append(x); ys.append(a.mean())
        rows.append(dict(M=c["M"], n=len(a), nu_mean=a.mean(), nu_sd=a.std(ddof=1), c_s_mass=a.mean() / x))
    if len(xs) < 3:
        return None
    xs, ys = np.array(xs), np.array(ys)
    cs = float((xs * ys).sum() / (xs * xs).sum()); sc = float(np.std(ys / xs, ddof=1))
    return dict(c_s=cs, scatter=sc, n_masses=len(xs), n_runs=sum(r["n"] for r in rows), masses=rows)


if __name__ == "__main__":
    roots = [r for r in sys.argv[1].split(",") if r]; prefix = sys.argv[2]
    tasks = [(root, c) for root in roots for c in sorted(glob.glob(f"{root}/eta_*/N*/m_*"))]
    with Pool(10) as pool:
        cells = pool.map(cell, tasks, chunksize=1)
    # merge leaves of the same (eta, N, M) across roots (top-up traces are distinct seeds)
    merged = {}
    for c in cells:
        k = (c["eta"], c["N"], c["M"])
        m = merged.setdefault(k, dict(eta=c["eta"], N=c["N"], M=c["M"], nu=[], nu_w3=[], ncyc=[], disc=0, L0=c["L0"], roots=[]))
        m["nu"] += c["nu"]; m["nu_w3"] += c["nu_w3"]; m["ncyc"] += c["ncyc"]; m["disc"] += c["disc"]; m["roots"].append(c["root"])
        if m["L0"] is None: m["L0"] = c["L0"]
    pts = {}
    for (eta, N) in sorted({(k[0], k[1]) for k in merged}):
        cs = [v for k, v in merged.items() if k[0] == eta and k[1] == N]
        a = fit_point(cs, "nu"); b = fit_point(cs, "nu_w3")
        if a:
            pts[(eta, N)] = dict(a, w3=b["c_s"] if b else float("nan"), L0=cs[0]["L0"],
                                 disc=sum(c["disc"] for c in cs), ncyc_min=min(min(c["ncyc"]) for c in cs if c["ncyc"]),
                                 ncyc_max=max(max(c["ncyc"]) for c in cs if c["ncyc"]))
    with open(prefix + "_cs_per_mass.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["eta", "N", "L0", "M", "n_runs", "nu_mean", "nu_sd", "c_s_mass"])
        for (eta, N), p in sorted(pts.items()):
            for r in p["masses"]:
                w.writerow([f"{eta:.2f}", N, f"{p['L0']:.6f}", r["M"], r["n"], f"{r['nu_mean']:.10g}", f"{r['nu_sd']:.6g}", f"{r['c_s_mass']:.6f}"])
    print("### Per (η, N): c_s = through-origin slope over masses, error = 1σ scatter of per-mass c_s")
    print("| η | N | L0 | masses | trajectories used / discarded | N_cyc | c_s (X = 2.5) | ± scatter | c_s (w = 3) | Δ X=2.5 vs w=3 [%] | KR | dev from KR [%] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    ext_rows = []
    for (eta, N), p in sorted(pts.items()):
        kr = T.kr_cs(eta)
        print(f"| {eta:.2f} | {N} | {p['L0']:.4f} | {p['n_masses']} | {p['n_runs']} / {p['disc']} | {p['ncyc_min']:.2f}–{p['ncyc_max']:.2f} | {p['c_s']:.4f} | {p['scatter']:.4f} | {p['w3']:.4f} | {100*(p['c_s']-p['w3'])/p['w3']:+.3f} | {kr:.4f} | {100*(p['c_s']-kr)/kr:+.2f} |")
    print("\n### Finite size: c_s = c_∞ + b/√N per η, weighted by the mass scatter")
    print("| η | N values | c_∞ | ± | b | ± | χ² | dof | KR | (c_∞ − KR)/KR [%] | (c_∞ − KR)/σ |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    fits = {}
    for eta in sorted({e for e, _ in pts}):
        rs = sorted((N, v) for (e, N), v in pts.items() if e == eta)
        if len(rs) < 3:
            continue
        Ns = np.array([N for N, _ in rs], float); c = np.array([v["c_s"] for _, v in rs]); s = np.array([max(v["scatter"], 1e-12) for _, v in rs])
        sl, sl_err, ic, ic_err = sos.weighted_linreg(1 / np.sqrt(Ns), c, s, force_zero_intercept=False)[:4]
        chi2 = float((((c - (ic + sl / np.sqrt(Ns))) / s) ** 2).sum()); dof = len(Ns) - 2
        kr = T.kr_cs(eta); fits[eta] = dict(ic=ic, ic_err=ic_err, sl=sl, sl_err=sl_err, Ns=Ns, c=c, s=s, kr=kr, chi2=chi2)
        ext_rows.append(dict(eta=f"{eta:.2f}", N_values=" ".join(str(int(n)) for n in Ns), c_inf=f"{ic:.5f}", c_inf_err=f"{ic_err:.5f}",
                             slope=f"{sl:.5f}", slope_err=f"{sl_err:.5f}", chi2=f"{chi2:.3f}", dof=dof, KR=f"{kr:.5f}",
                             dev_KR_pct=f"{100*(ic-kr)/kr:+.3f}", dev_KR_sigma=f"{(ic-kr)/ic_err:+.2f}"))
        print(f"| {eta:.2f} | {' '.join(str(int(n)) for n in Ns)} | {ic:.4f} | {ic_err:.4f} | {sl:.4f} | {sl_err:.4f} | {chi2:.3f} | {dof} | {kr:.4f} | {100*(ic-kr)/kr:+.2f} | {(ic-kr)/ic_err:+.2f} |")
    T.write_csv(prefix + "_cs_vs_N_extrapolation.csv", ext_rows)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    etas = sorted(fits); nc = min(3, len(etas)); nr = math.ceil(len(etas) / nc)
    fig, axs = plt.subplots(nr, nc, figsize=(4.6 * nc, 3.9 * nr), squeeze=False)
    for ax, eta in zip(axs.flat, etas):
        q = fits[eta]; u = 1 / np.sqrt(q["Ns"])
        ax.errorbar(u, q["c"], yerr=q["s"], fmt="o", color="#2a78d6", capsize=3, ms=6, label="EDMD A2, per-N c_s ± mass scatter")
        uu = np.linspace(0, u.max() * 1.08, 50)
        ax.plot(uu, q["ic"] + q["sl"] * uu, "-", color="#52514e", lw=1.2, label=f"c_∞ + b/√N, χ² = {q['chi2']:.2f}")
        ax.errorbar([0], [q["ic"]], yerr=[q["ic_err"]], fmt="s", color="#eb6834", capsize=3, ms=7, label=f"c_∞ = {q['ic']:.4f} ± {q['ic_err']:.4f}")
        ax.axhline(q["kr"], color="#e34948", lw=1.4, label=f"KR 2006 = {q['kr']:.4f}")
        for n, ui, ci in zip(q["Ns"], u, q["c"]):
            ax.annotate(f"N={int(n)}", (ui, ci), textcoords="offset points", xytext=(4, 6), fontsize=7, color="#52514e")
        ax.set_title(f"η = {eta:.2f}  ({100*(q['ic']-q['kr'])/q['kr']:+.2f} % vs KR, {(q['ic']-q['kr'])/q['ic_err']:+.1f}σ)", fontsize=10)
        ax.set_xlabel("1/√N"); ax.set_ylabel("c_s"); ax.grid(True, ls=":", alpha=0.6); ax.legend(fontsize=6.8, loc="best")
    for ax in list(axs.flat)[len(etas):]:
        ax.axis("off")
    fig.suptitle("A2, fixed η, N ladder: largest FFT bin at f ≥ ν_pred/2.5, mean over seeds, through-origin fit over masses", fontsize=10.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}_cs_vs_N.{ext}", dpi=170)
    print("\nhealth: trajectories used", sum(p["n_runs"] for p in pts.values()), "discarded", sum(p["disc"] for p in pts.values()))
    print("wrote", prefix + "_cs_vs_N.png/.pdf, _cs_per_mass.csv, _cs_vs_N_extrapolation.csv")
