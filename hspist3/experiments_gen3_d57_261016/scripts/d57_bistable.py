#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 57, addendum 3, task 2; notes 261012 sec. 4.7.140): THE BISTABLE RUN -- P2, eta 0.78, G2, pin 0,
alpha 2 (M 1800), run 9 (seed 3661301716), whose left compartment relaxed to e of about +2 %, against its nine siblings (runs 0-8, the same
cell and mass). No run.
  1. Per run: the per-trajectory c (sqrt(nu_d^2 + (1 / (2 pi tau_r))^2) / x: nu_d the registered estimator, tau_r its own fit;
     p2_tables_mac.py's c_seed), nu_d, tau_r, Q, the spectrum share, sd / sd_E1, the divider's mean displacement and drift, Z; psi6 and the
     orientation per compartment over the record (psi6_t); the strain per compartment at the start, the hold's end, over the record (mean) and
     at the end; the line direction and ordered fraction at the end (p2_runs.csv, p2_frames.csv, report_261016/tables/p2_seeds.csv).
  2. The (cell, mass) E1 exactly as s1_verdict_261012.per_mass (p2_tables_mac.per_mass's steps) with all 10 runs and without run 9, and
     the cell's E1 over the 5 masses (s1_verdict_261012.e1) with the alpha-2 value from each.
  3. A scan of all 800 compartments of P2: the record-mean strain of each against the median of its (cell, compartment) over the 50 runs, in %
     and in units of the thermal sd (the median of the runs' own sd over their record's snapshots); every compartment more than 1 % away.
  4. A figure: e(t) and psi6(t) per compartment for the ten runs of (0.78, G2, pin 0, alpha 2), run 9 marked, the hold's end dashed.
Writes bistable/: bistable_runs.csv, bistable_scan.csv, bistable.md, fig_d57_bistable.png.
usage (from hspist3/): python3 experiments_gen3_d57_261016/scripts/d57_bistable.py [--out <a new folder>]
"""
import glob, math, os, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d57_261016")
for p in (os.path.join(HS, "cluster", "gen3_stage2_261016"), os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "cluster", "gen3_testR_261011"),
          os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import tests_20260913 as T
import s1_verdict_261012 as V
from p2_tables_mac import seedavg
warnings.filterwarnings("ignore")
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016", "P2", "merged", "stage2", "gen3")
CELL, M, RUN = "P2_N0900_eta0.780000_G2_pin0", 1800, 9
SEEDS = os.path.join(HS, "experiments_gen3_p2_261016", "report_261016", "tables", "p2_seeds.csv")
MASSES = os.path.join(HS, "experiments_gen3_p2_261016", "report_261016", "tables", "p2_masses.csv")
OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else os.path.join(E, "bistable")


def per_mass(d, exclude=()):
    """p2_tables_mac.per_mass's steps (= s1_verdict_261012.per_mass) on the folder's runs less `exclude`"""
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run"); red = red[~red["run"].isin(exclude)]
    nud = pd.read_csv(os.path.join(d, "nud_runs.csv")); hl = pd.read_csv(os.path.join(d, "runs_health.csv")); z = np.load(os.path.join(d, "acf_runs.npz"))
    clean = set(hl[(hl["gen3_clean"] == 1) & (hl["gen3_records"] == 1)]["run"]); ok = set(nud[nud["fit_ok"] == 1]["run"])
    used = sorted(set(red["run"]) & clean & ok); dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    nd = nud.set_index("run").loc[used, "nu_d"].to_numpy(float); acfs = [z[f"run{r}"] for r in used]; p = seedavg(acfs, dt, nu_ref)
    nu0 = math.sqrt(nd.mean() ** 2 + (1 / (2 * math.pi * p[3])) ** 2); jk = []
    for i in range(len(used)):
        keep = [k for k in range(len(used)) if k != i]; fi = seedavg([acfs[k] for k in keep], dt, nu_ref)
        if fi is not None: jk.append(math.sqrt(nd[keep].mean() ** 2 + (1 / (2 * math.pi * fi[3])) ** 2))
    n = len(jk); se = math.sqrt((n - 1) / n * sum((v - np.mean(jk)) ** 2 for v in jk)) if n >= 3 else math.nan
    return dict(used=len(used), nud=float(nd.mean()), tau_r=float(p[3]), Q=math.pi * nd.mean() * p[3], nu0=nu0, nu0_se=se)


def main():
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    d = os.path.join(P2, CELL, f"m_{M}"); R = pd.read_csv(os.path.join(d, "p2_runs.csv")); F = pd.read_csv(os.path.join(d, "p2_frames.csv"))
    S = pd.read_csv(SEEDS); S = S[(S["cell"] == CELL) & (S["M"] == M)]
    rows = []
    for r in R.sort_values("run").to_dict("records"):
        s = S[S["run"] == r["run"]].iloc[0]; f = F[(F["run"] == r["run"]) & (F["phase"] == 1)]
        rows.append(dict(run=int(r["run"]), seed=int(r["seed"]), c_run=s["c_seed"], nu_d=r["nu_d"], tau_r=r["fit_tau_r"], Q=r["fit_Q"], spectrum_share=r["spectrum_share"],
                         sd_over_E1=s["sd_over_E1"], x_mean=r["x_mean"], drift_per_1000=r["drift_per_1000"], Z=r["Z"],
                         psi6_rec_0=r["psi6_rec_0"], psi6_rec_1=r["psi6_rec_1"], orient_dev_0=s["orient_rec_0_dev"], orient_dev_1=s["orient_rec_1_dev"],
                         e_start_0=r["start_e_0"], e_start_1=r["start_e_1"], e_hold_0=r["hold_e_0"], e_hold_1=r["hold_e_1"],
                         e_rec_0=f["e_0"].mean(), e_rec_1=f["e_1"].mean(), e_rec_sd_0=f["e_0"].std(ddof=1), e_rec_sd_1=f["e_1"].std(ddof=1),
                         e_end_0=r["end_e_0"], e_end_1=r["end_e_1"], line_end_0=r["end_line_0"], line_end_1=r["end_line_1"],
                         ordered_end_0=r["end_ordered_0"], ordered_end_1=r["end_ordered_1"], div_end=r["end_div_minus_L0"]))
    B = pd.DataFrame(rows); os.makedirs(OUT); B.to_csv(os.path.join(OUT, "bistable_runs.csv"), index=False, float_format="%.8g")
    # 2. E1 with and without run 9
    a = M / 900.0; L0 = float(R["L0"].iloc[0]); x = T.k_root(a) / (2 * math.pi * (L0 - 1.025))
    full, less = per_mass(d), per_mass(d, exclude=(RUN,))
    Mt = pd.read_csv(MASSES); Mt = Mt[Mt["cell"] == CELL]
    def cell_e1(c2, s2):
        rr = [dict(alpha=q["alpha"], flagged=False, c=(c2 if int(q["M"]) == M else q["c"]), sig_c=(s2 if int(q["M"]) == M else q["sig_c"]), Q=q["Q"], M=int(q["M"]))
              for q in Mt.to_dict("records")]
        return V.e1(rr)
    eA, eB = cell_e1(full["nu0"] / x, full["nu0_se"] / x), cell_e1(less["nu0"] / x, less["nu0_se"] / x)
    # 3. the scan of all 800 compartments
    sc = []
    for dd in sorted(glob.glob(os.path.join(P2, "*", "m_*"))):
        FF = pd.read_csv(os.path.join(dd, "p2_frames.csv")); FF = FF[FF["phase"] == 1]; cid = dd.split(os.sep)[-2]; MM = int(dd.split("_")[-1])
        for r_, g in FF.groupby("run"):
            for k in (0, 1):
                sc.append(dict(cell=cid, M=MM, run=int(r_), side=k, e_rec=g[f"e_{k}"].mean(), sd_t=g[f"e_{k}"].std(ddof=1)))
    SC = pd.DataFrame(sc)
    SC["median_cell_side"] = SC.groupby(["cell", "side"])["e_rec"].transform("median"); SC["thermal_cell_side"] = SC.groupby(["cell", "side"])["sd_t"].transform("median")
    SC["dev"] = SC["e_rec"] - SC["median_cell_side"]; SC["dev_in_thermal_sd"] = SC["dev"] / SC["thermal_cell_side"]
    SC.to_csv(os.path.join(OUT, "bistable_scan.csv"), index=False, float_format="%.8g"); far = SC[SC["dev"].abs() > 0.01]
    # the report
    sib = B[B["run"] != RUN]; me = B[B["run"] == RUN].iloc[0]
    def rng(col, k=1.0, dg=3):
        v = k * sib[col].to_numpy(float); return f"{np.median(v):.{dg}f} ({v.min():.{dg}f} to {v.max():.{dg}f})"
    items = [("c (this run's own)", "c_run", 1, 3), ("nu_d", "nu_d", 1, 5), ("tau_r (own fit)", "tau_r", 1, 1), ("Q (own fit)", "Q", 1, 1),
             ("spectrum share [%]", "spectrum_share", 100, 1), ("sd / sd_E1", "sd_over_E1", 1, 2), ("divider's mean displacement [sigma]", "x_mean", 1, 4),
             ("drift per 1000 [sigma]", "drift_per_1000", 1, 4), ("Z", "Z", 1, 2),
             ("psi6, record: left", "psi6_rec_0", 1, 3), ("psi6, record: right", "psi6_rec_1", 1, 3),
             ("orientation - pin [deg], record: left", "orient_dev_0", 1, 2), ("orientation - pin [deg], record: right", "orient_dev_1", 1, 2),
             ("e at the start [%]: left", "e_start_0", 100, 3), ("e at the start [%]: right", "e_start_1", 100, 3),
             ("e at the hold's end [%]: left", "e_hold_0", 100, 3), ("e at the hold's end [%]: right", "e_hold_1", 100, 3),
             ("e, the record's mean [%]: left", "e_rec_0", 100, 3), ("e, the record's mean [%]: right", "e_rec_1", 100, 3),
             ("e's sd over the record [%]: left", "e_rec_sd_0", 100, 3), ("e's sd over the record [%]: right", "e_rec_sd_1", 100, 3),
             ("e at the end [%]: left", "e_end_0", 100, 3), ("e at the end [%]: right", "e_end_1", 100, 3),
             ("ordered fraction at the end: left", "ordered_end_0", 1, 3), ("ordered fraction at the end: right", "ordered_end_1", 1, 3),
             ("line direction at the end (0 = x): left / right", None, 1, 0), ("divider - L0 at the end [sigma]", "div_end", 1, 4)]
    L = [f"# The bistable run (decision 57 addendum 3, task 2; d57_bistable.py): {CELL}, alpha {a:g} (M {M}), run {RUN} (seed {int(me['seed'])})", "",
         "| quantity | run 9 | its nine siblings: median (range) |", "|---|---|---|"]
    for lab, col, k, dg in items:
        if col is None:
            L.append(f"| {lab} | {int(me['line_end_0'])} / {int(me['line_end_1'])} | {int((sib['line_end_0'] == 0).sum())} of 9 x / {int((sib['line_end_1'] == 0).sum())} of 9 x |"); continue
        L.append(f"| {lab} | {k * me[col]:.{dg}f} | {rng(col, k, dg)} |")
    f9 = F[F["run"] == RUN].sort_values("t"); up = f9[f9["e_0"] > 0.01]
    L += ["", f"- Run {RUN}'s left compartment: the first snapshot with e > +1 % at t = {up['t'].iloc[0]:.0f} sigma-time (snapshots every "
              f"{f9['t'].diff().median():.0f}); e > 0 in {int((f9['e_0'] > 0).sum())} of its {len(f9)} snapshots (the hold and the record); its nine siblings' "
              f"left compartments: e > 0 in {int(((F['run'] != RUN) & (F['e_0'] > 0) & (F['frame'] > 0)).sum())} of their "
              f"{int(((F['run'] != RUN) & (F['frame'] > 0)).sum())} snapshots after the first."]
    L += ["", "**The (cell, mass) E1 with and without run 9** (s1_verdict_261012.per_mass's steps; the cell's E1 over the 5 masses by s1_verdict_261012.e1, "
              "the other four masses unchanged):", "",
          "| the runs | used | nu_d-bar | tau_r | Q | c (alpha 2) ± SE | the cell's E1 over the 5 masses ± SE (chi2 / dof) |", "|---|---|---|---|---|---|---|"]
    for lab, pm, ee in (("all 10", full, eA), ("without run 9", less, eB)):
        L.append(f"| {lab} | {pm['used']} | {pm['nud']:.5f} | {pm['tau_r']:.1f} | {pm['Q']:.1f} | {pm['nu0'] / x:.3f} ± {pm['nu0_se'] / x:.3f} | "
                 f"{ee['cs']:.3f} ± {ee['sig']:.3f} ({ee['chi2']:.2f} / {ee['dof']}) |")
    L += ["", f"**The scan of all {len(SC)} compartments** (the record-mean strain against the median of its cell and compartment over the cell's 50 runs): "
              f"{len(far)} lie more than 1 % away; the largest |deviation| in thermal sd: "
              f"{SC['dev_in_thermal_sd'].abs().max():.1f} (the next {SC['dev_in_thermal_sd'].abs().sort_values().iloc[-2]:.1f}).", ""]
    if len(far):
        L += ["| cell | M | run | compartment | e, the record's mean [%] | the median of its cell and compartment [%] | deviation [thermal sd] |", "|---|---|---|---|---|---|---|"]
        for r in far.to_dict("records"):
            L.append(f"| {r['cell'].replace('P2_N0900_', '')} | {r['M']} | {r['run']} | {'left' if r['side'] == 0 else 'right'} | {100 * r['e_rec']:+.3f} | "
                     f"{100 * r['median_cell_side']:+.3f} | {r['dev_in_thermal_sd']:+.1f} |")
    open(os.path.join(OUT, "bistable.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    # 4. the figure
    fig, axs = plt.subplots(2, 2, figsize=(15, 8.5), sharex=True)
    for k, side in ((0, "left"), (1, "right")):
        for r in sorted(F["run"].unique()):
            f = F[F["run"] == r].sort_values("t"); hot = r == RUN
            axs[0, k].plot(f["t"], 100 * f[f"e_{k}"], color="#c0392b" if hot else "0.6", lw=1.6 if hot else 0.8, alpha=1 if hot else 0.8,
                           label=f"run {RUN} (seed {int(me['seed'])})" if hot else ("the nine siblings" if r == 0 else None), zorder=3 if hot else 2)
            ps = pd.read_csv(glob.glob(os.path.join(d, f"psi6_t_L0_*_run{r}.csv*"))[0])
            axs[1, k].plot(ps["t_sigma"], ps[f"psi6_{k}"], color="#c0392b" if hot else "0.6", lw=1.0 if hot else 0.5, alpha=1 if hot else 0.7, zorder=3 if hot else 2)
        for ax in axs[:, k]: ax.axvline(5000, color="0.3", ls="--", lw=0.9); ax.grid(alpha=0.25); ax.tick_params(labelsize=8)
        axs[0, k].set_title(f"the {side} compartment: strain e per snapshot", fontsize=9); axs[0, k].set_ylabel("e  [%]", fontsize=9)
        axs[1, k].set_title(f"the {side} compartment: |psi6| (psi6_t)", fontsize=9); axs[1, k].set_ylabel("|psi6|", fontsize=9)
        axs[1, k].set_xlabel("time from the start  [sigma-time]  (the hold ends at the dashed line)", fontsize=9)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle(f"P2, eta 0.78, G2, pin 0, alpha {a:g}: the ten runs, run {RUN} in red", fontsize=10.5)
    fig.savefig(os.path.join(OUT, "fig_d57_bistable.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
