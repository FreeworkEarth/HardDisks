#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57, addendum 1, task 2; notes 261012 sec. 4.7.137): THE P2 REPORT ON THE MAC -- its tables and
figures, from the fetched tree (hspist3/experiments_gen3_stage2_koa_261016/P2; fetch_p2_koa.sh, runsheet block E). No interpretation: the
numbers, and the rules stated here. Writes into a NEW folder --out:
  tables/      p2_tables_mac.py on the same data, unchanged: E1 per (cell, mass) exactly as s1_verdict_261012.per_mass (identical to it on
               stage 1's folder N 900, eta 0.69, M 4500, 25 seeds), E1 per cell without stage 1's flag, the per-trajectory table p2_seeds.csv;
  p2_report.md the report's tables, in the order of task 2:
    0. per cell: runs, gen3 clean, reduction errors, fit failures (fit_ok 0), stalls, the arrangement and vacancies as run (the run header),
       the M4C line against the pin, the largest deviation of the record's orientation from the pin, seconds per run;
    1. E1 per (cell, G, pin, mass): used / present, nu_d-bar, tau_r, Q, c with its jackknife SE; stage 1's flag (more than 2 of the 10
       trajectories lost: the mass leaves the cell's E1; s1_verdict_261012.per_mass);
    2. E1 per cell over the masses: s1_verdict_261012.e1 (the inverse-variance mean of 5 unflagged masses with an SE, the SE scaled by
       sqrt(chi2/dof) when that exceeds 1), chi2/dof, stage 1's low-Q mark (a mass with Q < 5);
    3. the tests:
       (i)   eta 0.72: c(pin 0) - c(turned) per G; z = the difference / sqrt(sig_0^2 + sig_turned^2), for the cells' E1 and per mass;
       (ii)  the same at eta 0.78;
       (iii) G1 vs G2 per eta and pin, the same z;
       (iv)  the end strain e (the last snapshot, per compartment) over the 10 seeds of each (cell, mass), and over the cell's 50 runs: the
             mean, sd and range; the thermal sd: the median over those runs of each run's own sd of e over its record's snapshots
             (p2_frames.csv, phase 1); R = the sd over the seeds / the thermal sd. The rule: R <= 2, the seeds end in one state within its
             thermal fluctuation (reproducible); R > 2, not;
       (v)   Q per mass (the seed-averaged fit: E1's) and per trajectory (its own fit, p2_reduce.py): median and range; the stalls; the fit
             failures; the spectrum share (the periodogram's power within 0.8-1.2 nu_d): median and smallest;
       (vi)  nu_d / nu_pred (nu_d the registered item-4a estimator, nud_runs.csv; nu_pred the binary's own): per cell the median and range,
             per mass the median; and 200 nu_d / nu_pred, the damped periods in a record of 200 predicted periods;
       (vii) sd / sd_E1 and the drift per 1000 sigma-time: per cell the median and range.
  p2_cells_e1.csv, p2_tests.csv, p2_strain.csv   the numbers of 2, 3 (i)-(iii) and 3 (iv);
  fig_p2_traces_fits_eta0.72.png, fig_p2_traces_fits_eta0.78.png   one seed per cell (s = 0 at alpha 5), as Report D56-1's trace figures:
               the whole record (200 predicted periods) with its running mean over one predicted period, the first 10 predicted periods,
               the ACF and its item-4a fit, refitted here from acf_runs.npz exactly as p2_reduce.py fits it (the refit's nu_d and tau_r
               against p2_runs.csv's are printed: equal when scipy agrees);
  fig_p2_c_vs_alpha.png   c per mass with its SE against alpha, one panel per eta (colour: the pin; marker and line: G); the cells' E1 as
               horizontal lines with their SE band;
  fig_p2_end_strain.png   the end strain per seed, one panel per cell (colour: alpha; filled: the left compartment, open: the right); the
               cell's mean end strain +- the thermal sd grey; the start's lattice strain dashed (in the title when off the scale).
usage (from hspist3/): python3 cluster/gen3_stage2_261016/p2_report_mac.py --data <fetched P2 folder> --out <new folder>
"""
import glob, math, os, subprocess, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
for p in (HERE, os.path.join(CL, "gen3_stage1_261012"), os.path.join(CL, "gen3_testR_261011"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import tests_20260913 as T
import resched_testT_followup_261007 as FU
import s1_verdict_261012 as V
from paper1_populate_cs_err_20261002 import TD
from p2_reduce import load_trace, one
warnings.filterwarnings("ignore")

ETAS = (0.72, 0.78); GEOMS = ("G1", "G2"); PINS = ("0", "turned"); ALPHAS = (1, 2, 5, 10, 20); NSEED = 10; N = 900
FIG_ALPHA, FIG_SEED, R_RULE = 5, 0, 2.0
C_TR, C_FIT = "#1f3b5a", "#c0392b"                                  # Report D56-1's trace and fit colours
C_PIN = {"0": "#2a6fbb", "turned": "#d9822b"}; MK_G = {"G1": ("o", "-"), "G2": ("s", "--")}
GP = [(G, pin) for G in GEOMS for pin in PINS]; CELLS = [(eta, G, pin) for eta in ETAS for G, pin in GP]
M4C = {"0": "left lines along x, right lines along x", "turned": "left lines along y, right lines along y"}


def cid(eta, G, pin): return f"P2_N{N:04d}_eta{eta:.6f}_{G}_pin{pin}"
def fx(v, d=4): return "-" if v is None or not np.isfinite(v) else f"{v:.{d}f}"
def pmx(c, s, d=3): return "-" if not np.isfinite(c) else (f"{c:.{d}f} ± {s:.{d}f}" if np.isfinite(s) else f"{c:.{d}f} (no SE)")


def mr(v, d=3, scale=1.0):
    v = scale * np.asarray(v, float); v = v[np.isfinite(v)]
    return "-" if not len(v) else f"{np.median(v):.{d}f} ({v.min():.{d}f} to {v.max():.{d}f})"


def vac(v): return "? vacancies" if not np.isfinite(v) else (f"{int(v)} vacancy" if int(v) == 1 else f"{int(v)} vacancies")


def zdiff(a, b):
    d = a[0] - b[0]; s = math.hypot(a[1], b[1]); return d, s, (d / s if s > 0 else math.nan)


def mass_row(Mt, c, a):
    g = Mt[(Mt["cell"] == c) & np.isclose(Mt["alpha"], a)]
    return (float(g["c"].iloc[0]), float(g["sig_c"].iloc[0])) if len(g) else (math.nan, math.nan)


# ------------------------------------------------------------------------------------------------------------------- the tables
def e1_cells(Mt):
    Mt = Mt.copy(); Mt["lost"] = NSEED - Mt["used"]; Mt["flagged"] = Mt["lost"] > 2; E = {}
    for eta, G, pin in CELLS:
        rows = [dict(alpha=r["alpha"], flagged=bool(r["flagged"]), c=r.get("c", math.nan), sig_c=r.get("sig_c", math.nan), Q=r.get("Q", math.nan),
                     M=int(r["M"])) for r in Mt[Mt["cell"] == cid(eta, G, pin)].to_dict("records")]
        e = V.e1(rows) if rows else None
        usable = sum(1 for r in rows if r["alpha"] >= 1 and not r["flagged"] and np.isfinite(r["c"]) and np.isfinite(r["sig_c"]) and r["sig_c"] > 0)
        E[(eta, G, pin)] = dict(e or dict(cs=math.nan, sig=math.nan, chi2=math.nan, dof=max(0, len(rows) - 1), n=0, lowQ=False, masses=[]), usable=usable)
    return Mt, E


def tests_pz(Mt, E):
    out = []
    for eta in ETAS:
        for G in GEOMS:
            out.append(dict(test="i" if eta == 0.72 else "ii", eta=eta, by=G, a=cid(eta, G, "0"), b=cid(eta, G, "turned"),
                            ea=E[(eta, G, "0")], eb=E[(eta, G, "turned")]))
    for eta in ETAS:
        for pin in PINS:
            out.append(dict(test="iii", eta=eta, by=f"pin {pin}", a=cid(eta, "G1", pin), b=cid(eta, "G2", pin), ea=E[(eta, "G1", pin)], eb=E[(eta, "G2", pin)]))
    rows = []
    for t in out:
        d, s, z = zdiff((t["ea"]["cs"], t["ea"]["sig"]), (t["eb"]["cs"], t["eb"]["sig"]))
        row = dict(test=t["test"], eta=t["eta"], by=t["by"], a=t["a"], b=t["b"], c_a=t["ea"]["cs"], sig_a=t["ea"]["sig"], c_b=t["eb"]["cs"], sig_b=t["eb"]["sig"],
                   diff=d, sig_diff=s, z=z)
        for a in ALPHAS:
            row[f"z_alpha{a}"] = zdiff(mass_row(Mt, t["a"], a), mass_row(Mt, t["b"], a))[2]
        rows.append(row)
    return pd.DataFrame(rows)


def thermal(data):
    th = []
    for d in sorted(glob.glob(os.path.join(data, "merged", "stage2", "gen3", "*", "m_*"))):
        fp = os.path.join(d, "p2_frames.csv")
        if not os.path.exists(fp): continue
        F = pd.read_csv(fp); F = F[F["phase"] == 1]
        for r, g in F.groupby("run"):
            th.append(dict(cell=d.split(os.sep)[-2], M=int(d.split("_")[-1]), run=int(r), rec_frames=len(g),
                           sd_t_0=float(g["e_0"].std(ddof=1)), sd_t_1=float(g["e_1"].std(ddof=1))))
    return pd.DataFrame(th, columns=["cell", "M", "run", "rec_frames", "sd_t_0", "sd_t_1"])


def strain_table(S, TH):
    X = S.merge(TH, on=["cell", "M", "run"], how="left"); rows = []
    for eta, G, pin in CELLS:
        g = X[X["cell"] == cid(eta, G, pin)]
        for a in list(ALPHAS) + ["all"]:
            h = g if a == "all" else g[np.isclose(g["alpha"], a)]
            row = dict(cell=cid(eta, G, pin), eta=eta, G=G, pin=pin, alpha=str(a), n=len(h), rec_frames=float(h["rec_frames"].median()) if len(h) else math.nan)
            for k in (0, 1):
                v = h[f"end_e_{k}"].to_numpy(float); v = v[np.isfinite(v)]; t = h[f"sd_t_{k}"].to_numpy(float); t = t[np.isfinite(t)]
                sd = float(np.std(v, ddof=1)) if len(v) > 1 else math.nan; tm = float(np.median(t)) if len(t) else math.nan
                row.update({f"mean_{k}": float(v.mean()) if len(v) else math.nan, f"sd_{k}": sd, f"min_{k}": float(v.min()) if len(v) else math.nan,
                            f"max_{k}": float(v.max()) if len(v) else math.nan, f"thermal_{k}": tm, f"R_{k}": sd / tm if tm > 0 else math.nan})
            rows.append(row)
    return pd.DataFrame(rows), X


# ------------------------------------------------------------------------------------------------------------------- the figures
def refit(d, r):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")); z = np.load(os.path.join(d, "acf_runs.npz"))
    nl = min(len(z[f"run{k}"]) for k in red["run"]); a = np.asarray(z[f"run{r}"][:nl], float)
    p = FU.fit_acf(a, float(red["dt"].iloc[0]), float(red["nu"].mean())); cn = a / a[0]; lag = np.arange(len(cn)) * float(red["dt"].iloc[0])
    return p, cn, lag, FU._g(lag, *p)


def fig_traces(data, eta, S, path):
    fig = plt.figure(figsize=(20, 3.7 * len(GP) + 1.2))
    gs = fig.add_gridspec(len(GP), 3, width_ratios=[1.45, 1.0, 1.2], wspace=0.24, hspace=0.62, top=0.94, bottom=0.04)
    checks = []; M = int(round(FIG_ALPHA * N))
    for i, (G, pin) in enumerate(GP):
        c = cid(eta, G, pin); d = os.path.join(data, "merged", "stage2", "gen3", c, f"m_{M}")
        a0, a1, a2 = (fig.add_subplot(gs[i, j]) for j in range(3))
        pr = pd.read_csv(os.path.join(d, "p2_runs.csv")) if os.path.exists(os.path.join(d, "p2_runs.csv")) else pd.DataFrame()
        rw = pr[pr["run"] == FIG_SEED] if len(pr) else pr
        if not len(rw):
            a0.text(0.5, 0.5, f"{c}: no run s = {FIG_SEED} at alpha {FIG_ALPHA}", ha="center", va="center", transform=a0.transAxes); continue
        rw = rw.iloc[0]; sr = S[(S["cell"] == c) & (S["M"] == M) & (S["run"] == FIG_SEED)]
        sde = float(sr["sd_over_E1"].iloc[0]) if len(sr) else math.nan
        tr = one(d, f"wall_x_positions_L0_*_run{FIG_SEED}.csv*")
        if tr:
            t, x, nup = load_trace(tr); n = T._prefix(t, nup, TD); t, x = t[:n], x[:n]; ddt = (t[-1] - t[0]) / (len(t) - 1)
            w = max(1, int(1.0 / (nup * ddt))); sl = np.convolve(x, np.ones(w) / w, mode="same")
            a0.plot(t, x, color=C_TR, lw=0.35, alpha=0.55, label="divider displacement")
            a0.plot(t[w:-w], sl[w:-w], color=C_FIT, lw=1.4, label="running mean over one predicted period")
            k = t <= t[0] + 10 / nup; a1.plot(t[k], x[k], color=C_TR, lw=0.9)
        else:
            for ax in (a0, a1): ax.text(0.5, 0.5, "the trace is not in the fetched tree", ha="center", va="center", transform=ax.transAxes, fontsize=9)
        a0.set_title(f"{c.replace('P2_N0900_', '')}, alpha {FIG_ALPHA}, s = {FIG_SEED} (seed {int(rw['seed'])}): {rw.get('arrangement', '?')}, "
                     f"{vac(rw.get('n_vac', math.nan))}; H {rw['H']:.4f}, L0 {rw['L0']:.4f}\n"
                     f"the whole record ({TD} predicted periods): mean {rw['x_mean']:+.3f}, sd {rw['x_sd']:.3f}, drift {rw['drift_per_1000']:+.4f} per 1000; "
                     f"sd / sd_E1 {sde:.2f}", fontsize=8)
        a1.set_title("the first 10 predicted periods", fontsize=8)
        try:
            p, cn, lag, g = refit(d, FIG_SEED); nd = p[4] / (2 * math.pi); nup_r = float(rw["nu_pred"])
            a2.plot(lag * nup_r, cn, color=C_TR, lw=1.0, label="ACF of the displacement (normalised)")
            a2.plot(lag * nup_r, g, color=C_FIT, lw=1.2, ls="--", label="item-4a fit: A e^(-t/tau_T) + B e^(-t/tau_r) cos(wt + phi)")
            dn, dtr = nd / float(rw["fit_nu_d"]) - 1, p[3] / float(rw["fit_tau_r"]) - 1
            checks.append(dict(cell=c, M=M, run=FIG_SEED, nu_d_refit=nd, nu_d_koa=float(rw["fit_nu_d"]), tau_r_refit=p[3], tau_r_koa=float(rw["fit_tau_r"]),
                               rel_nu_d=dn, rel_tau_r=dtr))
            a2.set_title(f"nu_d {nd:.5f}, tau_r {p[3]:.1f}, Q {math.pi * nd * p[3]:.1f};  A {p[0]:.2f}, tau_T {p[1]:.0f}, B {p[2]:.2f}; spectrum share "
                         f"{100 * rw['spectrum_share']:.0f} %\nthe refit against KOA's (p2_runs.csv): nu_d {100 * dn:+.2g} %, tau_r {100 * dtr:+.2g} %", fontsize=8)
        except Exception as e:
            a2.text(0.5, 0.5, f"the refit raised {type(e).__name__}", ha="center", va="center", transform=a2.transAxes)
        a2.axhline(0, color="0.5", lw=0.6, ls=":"); a2.set_xlabel("lag  [predicted periods]", fontsize=8)
        a0.set_xlabel("record time  [sigma sqrt(m/kT)]", fontsize=8); a0.set_ylabel("displacement  [sigma]", fontsize=8); a1.set_xlabel("record time", fontsize=8)
        for ax in (a0, a1): ax.axhline(0, color="0.5", lw=0.6, ls=":")
        if i == 0: a0.legend(fontsize=7, loc="upper right", framealpha=0.9); a2.legend(fontsize=7, loc="upper right", framealpha=0.9)
        for ax in (a0, a1, a2): ax.tick_params(labelsize=7); ax.grid(alpha=0.25)
    fig.suptitle(f"P2, eta {eta}: the trace and the fit of one seed per cell (s = {FIG_SEED}, alpha {FIG_ALPHA})", fontsize=11, y=0.995)
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return checks


def fig_c(Mt, E, path):
    fig, axs = plt.subplots(1, len(ETAS), figsize=(15, 5.4))
    for ax, eta in zip(axs, ETAS):
        for j, (G, pin) in enumerate(GP):
            g = Mt[Mt["cell"] == cid(eta, G, pin)].sort_values("alpha"); mk, ls = MK_G[G]; col = C_PIN[pin]; off = 1 + 0.035 * (j - 1.5)
            ax.errorbar(g["alpha"] * off, g["c"], yerr=g["sig_c"].fillna(0), marker=mk, ls=ls, color=col, mfc=col if G == "G1" else "white", mec=col,
                        ms=6, capsize=3, lw=1.1, label=f"{G}, pin {pin}")
            e = E[(eta, G, pin)]
            if np.isfinite(e["cs"]):
                ax.axhline(e["cs"], color=col, ls=ls, lw=0.9, alpha=0.85); ax.axhspan(e["cs"] - e["sig"], e["cs"] + e["sig"], color=col, alpha=0.08, lw=0)
        ax.set_xscale("log"); ax.set_xticks(ALPHAS); ax.set_xticklabels([str(a) for a in ALPHAS]); ax.minorticks_off()
        ax.set_xlabel("alpha = M / N", fontsize=9); ax.set_ylabel("c  (E1 per mass, with its jackknife SE)", fontsize=9)
        ax.set_title(f"eta {eta}: the cells' E1 over the masses as horizontal lines (band: their SE)", fontsize=9); ax.grid(alpha=0.25); ax.legend(fontsize=8)
    fig.suptitle("P2: the E1 speed per mass and cell (colour: the pin; filled circles and solid lines: G1; open squares and dashed lines: G2)", fontsize=10)
    fig.savefig(path, dpi=120, bbox_inches="tight"); plt.close(fig)


def fig_strain(X, ST, path):
    fig, axs = plt.subplots(len(ETAS), len(GP), figsize=(20, 8.6)); cols = plt.cm.Blues(np.linspace(0.42, 0.95, len(ALPHAS)))
    for i, eta in enumerate(ETAS):
        for j, (G, pin) in enumerate(GP):
            ax = axs[i, j]; c = cid(eta, G, pin); g = X[X["cell"] == c]; al = ST[(ST["cell"] == c) & (ST["alpha"] == "all")]
            for ai, a in enumerate(ALPHAS):
                h = g[np.isclose(g["alpha"], a)]; xo = h["run"] + (ai - 2) * 0.13
                ax.plot(xo, 100 * h["end_e_0"], "o", color=cols[ai], ms=5, label=f"alpha {a}" if (i, j) == (0, 0) else None)
                ax.plot(xo, 100 * h["end_e_1"], "o", mfc="white", mec=cols[ai], ms=5)
            ys = 100 * g[["end_e_0", "end_e_1"]].to_numpy(float).ravel(); ys = ys[np.isfinite(ys)]
            if len(al):
                for k, ls_ in ((0, "-"), (1, ":")):
                    m, s = 100 * al[f"mean_{k}"].iloc[0], 100 * al[f"thermal_{k}"].iloc[0]
                    if np.isfinite(m) and np.isfinite(s):
                        ax.axhspan(m - s, m + s, color="0.5", alpha=0.12, lw=0); ax.axhline(m, color="0.45", lw=0.8, ls=ls_)
            el = 100 * float(g["e_lattice"].mean()) if len(g) else math.nan
            lo, hi = (ys.min(), ys.max()) if len(ys) else (-1.0, 1.0); pad = 0.15 * (hi - lo) + 0.05; ax.set_ylim(lo - pad, hi + pad)
            if np.isfinite(el) and lo - pad <= el <= hi + pad: ax.axhline(el, color="0.2", ls="--", lw=1.0, label="the start's lattice strain" if (i, j) == (0, 0) else None)
            note = "" if not np.isfinite(el) or lo - pad <= el <= hi + pad else " (off the scale)"
            ax.set_title(f"eta {eta}, {G}, pin {pin}: the start's lattice e {el:+.2f} %{note}", fontsize=8.5)
            ax.set_xticks(range(NSEED)); ax.set_xlabel("seed index s", fontsize=8); ax.set_ylabel("end strain e  [%]", fontsize=8)
            ax.axhline(0, color="0.7", lw=0.6, ls=":"); ax.tick_params(labelsize=7); ax.grid(alpha=0.2)
    axs[0, 0].legend(fontsize=7, loc="best", framealpha=0.9)
    fig.suptitle("P2: the end strain per seed (the last snapshot; filled: left compartment, open: right; colour: alpha). Grey: the cell's mean end "
                 "strain +- the thermal sd (solid: left, dotted: right); dashed: the start's lattice strain", fontsize=10)
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------- the report
def main():
    av = sys.argv; data = os.path.abspath(av[av.index("--data") + 1]); out = os.path.abspath(av[av.index("--out") + 1])
    if os.path.exists(out): sys.exit(f"REFUSED: {out} exists (a new folder only)")
    os.makedirs(out); tdir = os.path.join(out, "tables")
    pt = subprocess.run([sys.executable, os.path.join(HERE, "p2_tables_mac.py"), "--data", data, "--out", tdir], capture_output=True, text=True)
    open(os.path.join(out, "tables_stdout.txt"), "w").write(pt.stdout + pt.stderr)
    if pt.returncode != 0: sys.exit(f"STOP: p2_tables_mac.py exit {pt.returncode}: {(pt.stdout + pt.stderr)[-600:]}")
    Mt = pd.read_csv(os.path.join(tdir, "p2_masses.csv")); S = pd.read_csv(os.path.join(tdir, "p2_seeds.csv"))
    Mt, E = e1_cells(Mt); TZ = tests_pz(Mt, E); TH = thermal(data); ST, X = strain_table(S, TH)
    pd.DataFrame([dict(cell=cid(*k), eta=k[0], G=k[1], pin=k[2], E1_c=v["cs"], E1_sig=v["sig"], chi2=v["chi2"], dof=v["dof"], masses=v["n"],
                       lowQ=int(v["lowQ"])) for k, v in E.items()]).to_csv(os.path.join(out, "p2_cells_e1.csv"), index=False, float_format="%.8g")
    TZ.to_csv(os.path.join(out, "p2_tests.csv"), index=False, float_format="%.6g"); ST.to_csv(os.path.join(out, "p2_strain.csv"), index=False, float_format="%.6g")
    mc = os.path.join(data, "merged", ".merge_complete"); rs = os.path.join(data, "merged", "p2_reduce_summary.txt")
    L = ["# The P2 report tables (p2_report_mac.py; decision 57 addendum 1, task 2)", "",
         f"- data: {os.path.relpath(data, HS) if data.startswith(HS) else data}; merge record: {open(mc).readline().strip() if os.path.exists(mc) else 'MISSING'}",
         f"- p2_reduce on KOA: {open(rs).read().strip().splitlines()[-1] if os.path.exists(rs) else 'MISSING'}",
         f"- runs in the fetched tables: {len(S)} of {len(CELLS) * len(ALPHAS) * NSEED} in {Mt['cell'].nunique()} cells and {len(Mt)} (cell, mass) folders", ""]
    S["err"] = S["error"].fillna("").astype(str).str.len() > 0 if "error" in S else False
    L += ["## 0. Completeness and health per cell", "",
          "| eta | G | pin | runs | gen3 clean | reduction errors | fit failures | stalls | the arrangement as run (vacancies) | M4C as pinned | largest abs(orientation - pin) [deg] | seconds per run: median (max) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = S[S["cell"] == cid(eta, G, pin)]
        arr = "; ".join(f"{a} ({int(v)})" for a, v in sorted({(str(r["arrangement"]), r["n_vac"]) for r in g.to_dict("records") if np.isfinite(r.get("n_vac", math.nan))}))
        dev = np.abs(g[["orient_rec_0_dev", "orient_rec_1_dev"]].to_numpy(float)); dev = dev[np.isfinite(dev)]
        m4 = int(g["m4c"].fillna("").astype(str).map(lambda s: s.startswith(pin + " ") and M4C[pin] in s).sum())
        L.append(f"| {eta} | {G} | {pin} | {len(g)} | {int((g['clean'] == 1).sum())} | {int(g['err'].sum())} | {int((g['fit_ok'] == 0).sum())} | "
                 f"{int(g['stall'].fillna(0).sum())} | {arr or '-'} | {m4} of {len(g)} | {fx(dev.max() if len(dev) else math.nan, 2)} | "
                 f"{fx(g['seconds'].median(), 0)} ({fx(g['seconds'].max(), 0)}) |")
    L += ["", "The state per cell: means over the runs and both compartments (strain e in %; the start: the initial state; the hold's end and the end: "
          "the snapshots); psi6 and the orientation over the record (psi6_t); Z the hold's last quarter.", "",
          "| eta | G | pin | e of the start's lattice | e start | e at the hold's end | e at the end | psi6 (record) | orientation - pin [deg]: median | Z |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = S[S["cell"] == cid(eta, G, pin)]; mn = lambda cols: 100 * float(np.nanmean(g[cols].to_numpy(float))) if len(g) else math.nan
        od = g[["orient_rec_0_dev", "orient_rec_1_dev"]].to_numpy(float).ravel(); od = od[np.isfinite(od)]
        L.append(f"| {eta} | {G} | {pin} | {fx(mn(['e_lattice']), 3)} | {fx(mn(['start_e_0', 'start_e_1']), 3)} | {fx(mn(['hold_e_0', 'hold_e_1']), 3)} | "
                 f"{fx(mn(['end_e_0', 'end_e_1']), 3)} | {fx(mn(['psi6_rec_0', 'psi6_rec_1']) / 100, 3)} | {fx(np.median(od) if len(od) else math.nan, 2)} | "
                 f"{fx(g['Z'].mean() if len(g) else math.nan, 2)} |")
    L += ["", "## 1. E1 per (cell, G, pin, mass)", "",
          "Exactly as s1_verdict_261012.per_mass (the used trajectories: gen3 clean and fit_ok; nu_d-bar their mean; tau_r of the item-4a fit of their "
          "seed-averaged ACF; nu_0 = sqrt(nu_d-bar^2 + (1 / (2 pi tau_r))^2); Q = pi nu_d-bar tau_r; c = nu_0 / x, x = K(alpha) / (2 pi (L0 - 1.025)); "
          "the delete-one jackknife SE). Flag: more than 2 of the 10 trajectories lost (stage 1's rule): the mass leaves the cell's E1.", "",
          "| eta | G | pin | alpha | used / present | flag | nu_d-bar | tau_r | Q | c ± SE |", "|---|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        for r in Mt[Mt["cell"] == cid(eta, G, pin)].sort_values("alpha").to_dict("records"):
            L.append(f"| {eta} | {G} | {pin} | {r['alpha']:g} | {int(r['used'])} / {int(r['present'])} | {'FLAG' if r['flagged'] else '-'} | "
                     f"{fx(r.get('nud', math.nan), 5)} | {fx(r.get('tau_r', math.nan), 1)} | {fx(r.get('Q', math.nan), 1)} | {pmx(r.get('c', math.nan), r.get('sig_c', math.nan))} |")
    L += ["", "## 2. E1 per cell over the masses", "", "s1_verdict_261012.e1: the inverse-variance mean of the 5 unflagged masses with an SE; the SE "
          "scaled by sqrt(chi2/dof) when that exceeds 1; no value with fewer than 5 such masses. Low Q: a mass with Q < 5.", "",
          "| eta | G | pin | H [sigma] | L0 [sigma] | masses in E1 | E1: c ± SE | chi2 / dof | low Q |", "|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        e = E[(eta, G, pin)]; g = Mt[Mt["cell"] == cid(eta, G, pin)]
        nm = str(e["n"]) if np.isfinite(e["cs"]) else f"none: {e['usable']} unflagged with an SE, 5 needed"
        L.append(f"| {eta} | {G} | {pin} | {fx(g['H'].iloc[0] if len(g) else math.nan, 4)} | {fx(g['L0'].iloc[0] if len(g) else math.nan, 4)} | {nm} | "
                 f"{pmx(e['cs'], e['sig'])} | {fx(e['chi2'], 2)} / {e['dof']} | {'yes' if e['lowQ'] else 'no'} |")
    L += ["", "## 3. The tests", "", "z = the difference / sqrt(SE_a^2 + SE_b^2), for the cells' E1 (section 2) and per mass (section 1).", ""]
    for tt, head in (("i", "(i) eta 0.72: pin 0 against turned, per G"), ("ii", "(ii) eta 0.78: pin 0 against turned, per G"),
                     ("iii", "(iii) G1 against G2, per eta and pin")):
        g = TZ[TZ["test"] == tt]; a_lab, b_lab = ("pin 0", "turned") if tt != "iii" else ("G1", "G2")
        L += [f"### {head}", "", f"| eta | {'G' if tt != 'iii' else 'pin'} | c({a_lab}) | c({b_lab}) | difference ± SE | z | z per mass (alpha 1, 2, 5, 10, 20) |",
              "|---|---|---|---|---|---|---|"]
        for r in g.to_dict("records"):
            L.append(f"| {r['eta']} | {r['by']} | {pmx(r['c_a'], r['sig_a'])} | {pmx(r['c_b'], r['sig_b'])} | {pmx(r['diff'], r['sig_diff'])} | {fx(r['z'], 2)} | "
                     + ", ".join(fx(r[f'z_alpha{a}'], 2) for a in ALPHAS) + " |")
        L.append("")
    L += ["### (iv) the end strain over the seeds", "",
          f"The end strain: the last snapshot, per compartment (p2_lattice.config_strain). sd: over the 10 seeds of one mass (or over the cell's 50 runs); "
          f"the thermal sd: the median over the same runs of each run's own sd of e over its record's snapshots; R = sd / the thermal sd. The rule: "
          f"R <= {R_RULE:g}: the seeds end in one state within its thermal fluctuation (reproducible); R > {R_RULE:g}: not. Strains in %.", "",
          "| eta | G | pin | compartment | end e: mean over the cell's runs (range) | sd over the 10 seeds at alpha 1, 2, 5, 10, 20 | thermal sd (median) | "
          "R at alpha 1, 2, 5, 10, 20 | R over the cell's runs | R <= 2 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = ST[ST["cell"] == cid(eta, G, pin)]; al = g[g["alpha"] == "all"]; pm_ = g[g["alpha"] != "all"]
        for k, side in ((0, "left"), (1, "right")):
            if not len(al): continue
            a = al.iloc[0]; ok = int((pm_[f"R_{k}"] <= R_RULE).sum()); nf = int(np.isfinite(pm_[f"R_{k}"].to_numpy(float)).sum())
            L.append(f"| {eta} | {G} | {pin} | {side} | {fx(100 * a[f'mean_{k}'], 3)} ({fx(100 * a[f'min_{k}'], 2)} to {fx(100 * a[f'max_{k}'], 2)}) | "
                     + ", ".join(fx(100 * v, 3) for v in pm_[f"sd_{k}"]) + f" | {fx(100 * a[f'thermal_{k}'], 3)} | " + ", ".join(fx(v, 2) for v in pm_[f"R_{k}"])
                     + f" | {fx(a[f'R_{k}'], 2)} | {f'{ok} of {nf} masses' if nf else '-'} |")
    L += ["", f"- record snapshots per run (median over each cell's 50): " + ", ".join(
        f"{cid(*k).replace('P2_N0900_', '')} {fx(ST[(ST['cell'] == cid(*k)) & (ST['alpha'] == 'all')]['rec_frames'].iloc[0], 0)}" for k in CELLS), ""]
    L += ["### (v) Q, the stalls and the spectrum share per cell", "",
          "| eta | G | pin | Q per mass, the seed-averaged fit (alpha 1, 2, 5, 10, 20) | Q per trajectory, its own fit: median (range) | stalls / runs | "
          "fit failures | spectrum share: median (smallest) |", "|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = S[S["cell"] == cid(eta, G, pin)]; m = Mt[Mt["cell"] == cid(eta, G, pin)].sort_values("alpha"); sh = 100 * g["spectrum_share"].to_numpy(float)
        sh = sh[np.isfinite(sh)]
        L.append(f"| {eta} | {G} | {pin} | " + ", ".join(fx(v, 1) for v in m.get("Q", pd.Series(dtype=float))) + f" | {mr(g['fit_Q'], 1)} | "
                 f"{int(g['stall'].fillna(0).sum())} / {len(g)} | {int((g['fit_ok'] == 0).sum())} | "
                 f"{fx(np.median(sh) if len(sh) else math.nan, 1)} % ({fx(sh.min() if len(sh) else math.nan, 1)} %) |")
    L += ["", "### (vi) nu_d / nu_pred per cell (the run-length rule)", "",
          "| eta | G | pin | nu_d / nu_pred: median (range) over the runs | median per mass (alpha 1, 2, 5, 10, 20) | damped periods in the record, 200 nu_d / nu_pred: median (range) |",
          "|---|---|---|---|---|---|"]
    S["nud_over_nup"] = S["nu_d"] / S["nu_pred"]
    for eta, G, pin in CELLS:
        g = S[S["cell"] == cid(eta, G, pin)]
        L.append(f"| {eta} | {G} | {pin} | {mr(g['nud_over_nup'])} | " + ", ".join(fx(g[np.isclose(g['alpha'], a)]['nud_over_nup'].median(), 3) for a in ALPHAS)
                 + f" | {mr(g['nud_over_nup'], 0, 200.0)} |")
    L += ["", "### (vii) sd / sd_E1 and the drift per cell", "",
          "sd: the trace's sd over the record; sd_E1 = sqrt(L_eff / (2 H rho c^2)) with c the (cell, mass) E1 (p2_tables_mac.py); the drift: the least-squares "
          "slope of the trace, per 1000 sigma-time [sigma].", "",
          "| eta | G | pin | sd / sd_E1: median (range) | drift per 1000: median (range) |", "|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = S[S["cell"] == cid(eta, G, pin)]
        L.append(f"| {eta} | {G} | {pin} | {mr(g['sd_over_E1'], 2)} | {mr(g['drift_per_1000'], 4)} |")
    ck = []
    for eta in ETAS: ck += fig_traces(data, eta, S, os.path.join(out, f"fig_p2_traces_fits_eta{eta}.png"))
    fig_c(Mt, E, os.path.join(out, "fig_p2_c_vs_alpha.png")); fig_strain(X, ST, os.path.join(out, "fig_p2_end_strain.png"))
    C = pd.DataFrame(ck)
    L += ["", "## Figures", "",
          f"- fig_p2_traces_fits_eta0.72.png, fig_p2_traces_fits_eta0.78.png: s = {FIG_SEED} at alpha {FIG_ALPHA} in each cell. The ACF refitted here from "
          f"acf_runs.npz as p2_reduce.py fits it, against KOA's p2_runs.csv over the {len(C)} runs shown: largest abs relative difference nu_d "
          f"{fx(C['rel_nu_d'].abs().max() if len(C) else math.nan, 6)}, tau_r {fx(C['rel_tau_r'].abs().max() if len(C) else math.nan, 6)}.",
          "- fig_p2_c_vs_alpha.png: c per mass and the cells' E1.", "- fig_p2_end_strain.png: the end strain per seed."]
    open(os.path.join(out, "p2_report.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
