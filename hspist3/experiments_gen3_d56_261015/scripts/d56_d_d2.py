#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56, task 2 d; notes 261012 sec. 4.7.133): THE D2 SIZING (eta 0.02, 0.035, 0.075 at N 900 and
1600, block D's registered protocol: T_eq 1000 sigma-time (60000 held steps), 200 predicted periods, nine masses, 25 seeds) AND BLOCK D'S
OWN FIT RECORD. Planning numbers and information; nothing registered.
  1. The boxes: H = 10 sqrt(N/100), L0 = N pi / (8 H eta) per side; the mean free path of the dilute 2D gas, lambda = pi / (4 sqrt 2 eta)
     sigma (lambda = 1 / (sqrt 2 n sigma), n = 4 eta / pi), against H.
  2. The period from c_s = sqrt 2 and L0: P = 2 pi L_eff / (sqrt 2 K(alpha)), L_eff = L0 - 1.025 (E1's x_k); 200 periods; the run, T_eq +
     200 P. The binary's own nu_pred (Henderson, as the protocol runs it) for comparison.
  3. Events and minutes per run, extrapolated from block D (N 100 and 400; 1350 runs): the event rate per disk and sigma-time fitted as
     r = a (Z - 1) + c_w (1/H + 1/L0) (the disk-disk collisions and the wall and divider hits of a dilute gas, whose rate per disk is
     sqrt(kT / (2 pi m)) x the perimeter over the area), Z block D's measured value per eta; the time per event fitted as
     u = u0 + u1 n_x, n_x = (the mean speed sqrt(pi/2) / r) / (the cell, 32 px = 4/3 sigma), the cells a disk crosses between two events,
     on block D's E5-2670 v3 runs. Both checked against B900 and B1600 at eta 0.05 and 0.10 (outside the fit).
  4. Block D's fit record, per (cell, mass), from the fetched tables: the trajectories, the fit failures of the registered fit (KOA's
     nud_runs.csv), the group rule's flag (more than 2 lost), Q_k = pi nu_d-bar tau_r (tau_r from the seed-averaged ACF over the used
     runs, as s1_verdict_261012.per_mass), tau_r against the fitted window (20 predicted periods); the deviation of the damped cosine
     from the data: the item-4a fit (FU.fit_acf, unchanged) of every trajectory's own ACF where KOA's fit converged, the rms and the
     largest |normalised ACF - fit| over the window, and the same for the seed-averaged ACF; the fits that "converge" without leaving
     the initial guess (A 0.3 and B 0.7 within 1e-3: a stall) and the fits that miss the ACF (rms > 0.15); KOA's nu_d of the same runs,
     to show that the registered values carry the same fits; the damping correction 1 / (8 Q^2) that E1 applies.
Writes into experiments_gen3_d56_261015/d_d2/ (a new folder): d2_table.csv, d2_table.md, d_fit_record.csv, d_runs_fit.csv, d_summary.txt,
fig_d1_block_D_acf_fits.png, fig_d2_block_D_Q_and_deviation.png. Six processes for the refits.
usage (from hspist3/): python3 experiments_gen3_d56_261015/scripts/d56_d_d2.py [--out <a new folder instead>]
"""
import concurrent.futures as cf, glob, math, os, re, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d56_261015")
for p in (os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import tests_20260913 as T
import resched_testT_followup_261007 as FU
warnings.filterwarnings("ignore")
KOA = os.path.join(HS, "experiments_gen3_stage1_koa_261012"); TIMING = os.path.join(HS, "experiments_gen3_s1_261012", "interim_261013", "timing")
OUT = os.path.join(E, "d_d2")
ETAS, NS, CELL = (0.02, 0.035, 0.075), (900, 1600), 32.0 / 24.0
V3 = "Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz"


def t_rec(nup): return math.ceil((200.0 / nup) / GEN.DTF) * GEN.DTF


def runs(groups, eta_max):
    tm = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(TIMING, "timing_*.tsv")))]); cpu = dict(zip(tm["node"], tm["cpu"]))
    rows = []
    for g in groups:
        for m in sorted(glob.glob(os.path.join(KOA, g, "merged", "stage1", "gen3", "*", "m_*"))):
            cell = m.split(os.sep)[-2]; N = int(cell[3:7]); eta = float(cell.split("eta")[1])
            if eta > eta_max + 1e-12: continue
            M = int(m.split("_")[-1]); h = pd.read_csv(os.path.join(m, "runs_health.csv")); s1 = pd.read_csv(os.path.join(m, "s1_runs.csv")); red = pd.read_csv(os.path.join(m, "red_nu.csv"))
            log = open(os.path.join(m, "run.log"), errors="replace").read()
            ev = {int(a): int(b) for a, b in re.findall(r"^##RUN \S+ \S+ run (\d+) seed \d+ .*?executed events (\d+)", log, re.M | re.S)}
            hold = {int(a): int(b) for a, b in re.findall(r"^##RUN \S+ \S+ run (\d+) seed \d+ .*?^Wall hold steps: (\d+)$", log, re.M | re.S)}
            H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta)
            for r in h.itertuples():
                rr = red[red["run"] == r.run]; ss = s1[s1["run"] == r.run]
                if not len(rr) or r.run not in ev: continue
                Tt = hold[r.run] * 0.4 / 24.0 + t_rec(float(rr["nu_pred"].iloc[0]))
                rows.append(dict(group=g, N=N, eta=eta, M=M, run=int(r.run), seconds=int(r.seconds), cpu=cpu.get(r.node, "?"), events=ev[r.run], T=Tt,
                                 Z=float(ss["Z"].iloc[0]), H=H, L0=L0))
    D = pd.DataFrame(rows); D["rate"] = D["events"] / (D["N"] * D["T"]); D["us"] = 1e6 * D["seconds"] / D["events"]
    D["wall"] = 1.0 / D["H"] + 1.0 / D["L0"]; D["n_x"] = math.sqrt(math.pi / 2) / D["rate"] / CELL
    return D


def fit_lin(X, y):
    c, *_ = np.linalg.lstsq(X, y, rcond=None); return c


def refit(args):
    """the registered item-4a fit of one trajectory's ACF; the deviation of the fitted damped cosine from the normalised ACF"""
    a, dt, nu_ref = args
    try:
        p = FU.fit_acf(a, dt, nu_ref)
    except Exception:
        return None
    c = a / a[0]; g = FU._g(np.arange(len(c)) * dt, *p); d = c - g
    return dict(nu_d=p[4] / (2 * math.pi), tau_r=p[3], A=p[0], tau_T=p[1], B=p[2], rms=float(np.sqrt(np.mean(d ** 2))), maxdev=float(np.max(np.abs(d))),
                rms_late=float(np.sqrt(np.mean(d[len(d) // 2:] ** 2))))


def fit_record():
    out, per_run, curves = [], [], {}
    folders = sorted(glob.glob(os.path.join(KOA, "D", "merged", "stage1", "gen3", "*", "m_*")), key=lambda p: (p.split(os.sep)[-2], int(p.split("_")[-1])))
    jobs, meta = [], []
    for m in folders:
        cell = m.split(os.sep)[-2]; M = int(m.split("_")[-1]); N = int(cell[3:7])
        nud = pd.read_csv(os.path.join(m, "nud_runs.csv")); red = pd.read_csv(os.path.join(m, "red_nu.csv")); z = np.load(os.path.join(m, "acf_runs.npz"))
        dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean()); nl = min(len(z[f"run{r}"]) for r in red["run"])
        for r, nk in zip(nud[nud["fit_ok"] == 1]["run"], nud[nud["fit_ok"] == 1]["nu_d"]):
            jobs.append((np.asarray(z[f"run{r}"][:nl], float), dt, nu_ref)); meta.append((cell, N, M, int(r), float(nk), nu_ref))
    with cf.ProcessPoolExecutor(6) as ex: res = list(ex.map(refit, jobs, chunksize=4))
    for (cell, N, M, r, nk, nref), x in zip(meta, res):
        per_run.append(dict(cell=cell, N=N, M=M, alpha=M / N, run=r, nu_d_koa=nk, nu_ref=nref, refit_ok=int(x is not None), **(x or {})))
    PR = pd.DataFrame(per_run)
    PR["stall"] = ((PR["A"] - 0.3).abs() < 1e-3) & ((PR["B"] - 0.7).abs() < 1e-3); PR["miss"] = PR["rms"] > 0.15
    PR["nu_d_vs_koa"] = PR["nu_d"] / PR["nu_d_koa"] - 1
    for m in folders:
        cell = m.split(os.sep)[-2]; M = int(m.split("_")[-1]); N = int(cell[3:7]); eta = float(cell.split("eta")[1])
        nud = pd.read_csv(os.path.join(m, "nud_runs.csv")); red = pd.read_csv(os.path.join(m, "red_nu.csv")); z = np.load(os.path.join(m, "acf_runs.npz"))
        h = pd.read_csv(os.path.join(m, "runs_health.csv")); clean = set(h[(h["gen3_clean"] == 1) & (h["gen3_records"] == 1)]["run"])
        used = sorted(set(nud[nud["fit_ok"] == 1]["run"]) & clean); dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
        nl = min(len(z[f"run{r}"]) for r in red["run"]); mean = np.mean([np.asarray(z[f"run{r}"][:nl], float) for r in used], axis=0)
        p = FU.fit_acf(mean, dt, nu_ref); c = mean / mean[0]; g = FU._g(np.arange(nl) * dt, *p); d = c - g
        nd = float(nud.set_index("run").loc[used, "nu_d"].mean()); win = nl * dt
        pr = PR[(PR.cell == cell) & (PR.M == M) & (PR.refit_ok == 1)]
        out.append(dict(cell=cell, N=N, eta=eta, M=M, alpha=M / N, present=len(red), fit_failures=int((nud["fit_ok"] == 0).sum()), used=len(used),
                        flagged=int(25 - len(used) > 2), Q_k=math.pi * nd * p[3], tau_r=p[3], window=win, tau_r_over_window=p[3] / win,
                        decay_over_window=math.exp(-win / p[3]), meanfit_rms=float(np.sqrt(np.mean(d ** 2))), meanfit_maxdev=float(np.max(np.abs(d))),
                        run_rms_median=float(pr["rms"].median()), run_rms_max=float(pr["rms"].max()), run_maxdev_median=float(pr["maxdev"].median()),
                        run_Q_median=float((math.pi * pr["nu_d"] * pr["tau_r"]).median()), run_Q_p10=float((math.pi * pr["nu_d"] * pr["tau_r"]).quantile(0.1)),
                        run_Q_p90=float((math.pi * pr["nu_d"] * pr["tau_r"]).quantile(0.9)), refit_failures=int((PR[(PR.cell == cell) & (PR.M == M)]["refit_ok"] == 0).sum()),
                        stalls=int(pr["stall"].sum()), misses=int(pr["miss"].sum()), damping_corr=1.0 / (8.0 * (math.pi * nd * p[3]) ** 2)))
        curves[(cell, M)] = (np.arange(nl) * dt * nu_ref, c, g)
    return pd.DataFrame(out), PR, curves


def main():
    global OUT
    if "--out" in sys.argv: OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    os.makedirs(OUT); L = ["# d56_d_d2.py: the D2 sizing and block D's fit record (decision 56 task 2 d)", ""]
    Dd = runs(("D",), 0.08); Bv = runs(("B900", "B1600"), 0.10)
    a, cw = fit_lin(np.vstack([Dd["Z"] - 1, Dd["wall"]]).T, Dd["rate"].to_numpy())
    v3 = Dd[Dd["cpu"] == V3]; u0, u1 = fit_lin(np.vstack([np.ones(len(v3)), v3["n_x"]]).T, v3["us"].to_numpy())
    L += ["## The rate and time models, fitted on block D", "",
          f"- r = {a:.4f} (Z - 1) + {cw:.4f} (1/H + 1/L0)   (block D's {len(Dd)} runs; a dilute gas's wall rate would be sqrt(1 / 2 pi) = 0.399 per wall: "
          f"0.798 for the two walls of each direction)",
          f"- microseconds per event on E5-2670 v3: u = {u0:.3f} + {u1:.4f} n_x   ({len(v3)} runs; n_x the cells crossed between two events)", "",
          "per cell, measured against the models (B900 and B1600 at eta 0.05 and 0.10 are outside the fit):", ""]
    allc = pd.concat([Dd, Bv]); pc = allc.groupby(["group", "N", "eta"]).agg(Z=("Z", "mean"), rate=("rate", "mean"), wall=("wall", "first"), n_x=("n_x", "mean")).reset_index()
    us_c = allc[allc["cpu"] == V3].groupby(["group", "N", "eta"])["us"].median()
    for r in pc.itertuples():
        rm = a * (r.Z - 1) + cw * r.wall; um = u0 + u1 * r.n_x; uo = us_c.get((r.group, r.N, r.eta), math.nan)
        L.append(f"  {r.group:5s} N {r.N:5d} eta {r.eta:.3f}: r {r.rate:.4f} (model {rm:.4f}, {100 * (rm / r.rate - 1):+.1f} %); microseconds per event (v3) "
                 f"{uo:.2f} (model {um:.2f}, {100 * (um / uo - 1):+.1f} %)")
    Zd = Dd.groupby("eta")["Z"].mean().to_dict(); rows = []
    for N in NS:
        H = 10.0 * math.sqrt(N / 100.0)
        for eta in ETAS:
            L0 = N * math.pi / (8.0 * H * eta); Leff = L0 - 1.025; lam = math.pi / (4 * math.sqrt(2) * eta)
            r = a * (Zd[eta] - 1) + cw * (1 / H + 1 / L0); nx = math.sqrt(math.pi / 2) / r / CELL; us = u0 + u1 * nx
            per = []
            for al in GEN.ALPHA:
                P = 2 * math.pi * Leff / (math.sqrt(2) * T.k_root(al)); nup = GEN.nu_pred_binary(N, H, L0, al * N)
                per.append(dict(alpha=al, P=P, T=1000.0 + 200 * P, P_pred=1 / nup, T_pred=1000.0 + t_rec(nup)))
            P5 = [x for x in per if x["alpha"] == 5][0]; Ts = np.array([x["T"] for x in per])
            ev5 = r * N * P5["T"]; ev9 = r * N * Ts
            rows.append(dict(N=N, eta=eta, H=H, L0=L0, box_w=2 * L0, box_h=H, lam=lam, Kn_H=lam / H, Z_D=Zd[eta], rate=r, us_per_event=us,
                             P_alpha5=P5["P"], P_min9=min(x["P"] for x in per), P_max9=max(x["P"] for x in per), T200_alpha5=200 * P5["P"],
                             T_run_alpha5=P5["T"], T_run_mean9=float(Ts.mean()), P_pred_alpha5=P5["P_pred"], T_run_pred_alpha5=P5["T_pred"],
                             events_alpha5=ev5, events_mean9=float(ev9.mean()), min_alpha5=ev5 * us / 6e7, min_mean9=float(ev9.mean()) * us / 6e7,
                             min_max9=float(ev9.max()) * us / 6e7, coreh_9x25=float(ev9.sum()) * 25 * us / 3.6e9))
    S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "d2_table.csv"), index=False, float_format="%.6g")
    md = ["| eta | N | box 2 L0 x H [sigma] | mean free path [sigma] (/ H) | period alpha 5 (9 masses) [sigma-time] | 200 periods, alpha 5 | run T_eq + 200 P, "
          "alpha 5 (mean of 9) | events per run, alpha 5 (mean of 9) | microseconds per event (v3) | minutes per run, alpha 5 (mean; largest) | core-hours 9 x 25 |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in S.itertuples():
        md.append(f"| {r.eta:.3f} | {r.N} | {r.box_w:.1f} x {r.box_h:.0f} | {r.lam:.1f} ({r.Kn_H:.2f}) | {r.P_alpha5:.0f} ({r.P_min9:.0f}-{r.P_max9:.0f}) | {r.T200_alpha5:.3g} | "
                  f"{r.T_run_alpha5:.3g} ({r.T_run_mean9:.3g}) | {r.events_alpha5 / 1e6:.0f} M ({r.events_mean9 / 1e6:.0f} M) | {r.us_per_event:.2f} | "
                  f"{r.min_alpha5:.0f} ({r.min_mean9:.0f}; {r.min_max9:.0f}) | {r.coreh_9x25:.0f} |")
    md.append(f"| all 6 | | | | | | | | | | {S.coreh_9x25.sum():.0f} |")
    L += ["", "## D2 (planning; block D's protocol: T_eq 1000, 200 periods; the period from c_s = sqrt 2; minutes on E5-2670 v3; "
          "the binary's Henderson nu_pred gives periods about 4 % shorter)", ""] + md
    R, PR, curves = fit_record(); R.to_csv(os.path.join(OUT, "d_fit_record.csv"), index=False, float_format="%.6g"); PR.to_csv(os.path.join(OUT, "d_runs_fit.csv"), index=False, float_format="%.6g")
    L += ["", "## Block D's fit record (per cell and mass; the registered item-4a fit; the window = 20 predicted periods of the ACF)", "",
          "| cell | alpha | fit failures (KOA) | used | flag | Q_k | tau_r / window | decay over the window | seed-avg fit: rms, largest dev | per run: rms median (max), Q median (10-90 %) | stalls, misses |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in R.itertuples():
        L.append(f"| {r.cell} | {r.alpha:g} | {r.fit_failures} | {r.used} | {'FLAG' if r.flagged else '-'} | {r.Q_k:.1f} | {r.tau_r_over_window:.2f} | {r.decay_over_window:.2f} | "
                 f"{r.meanfit_rms:.4f}, {r.meanfit_maxdev:.3f} | {r.run_rms_median:.4f} ({r.run_rms_max:.3f}), {r.run_Q_median:.0f} ({r.run_Q_p10:.0f}-{r.run_Q_p90:.0f}) | {r.stalls}, {r.misses} |")
    pcq = R.groupby("cell").agg(Q_median=("Q_k", "median"), Q_min=("Q_k", "min"), Q_max=("Q_k", "max"), fails=("fit_failures", "sum"), flags=("flagged", "sum"),
                                dev=("meanfit_rms", "median"), run_dev=("run_rms_median", "median")).reset_index()
    L += ["", "per cell: " + "; ".join(f"{r.cell}: Q_k median {r.Q_median:.0f} ({r.Q_min:.0f}-{r.Q_max:.0f}), {r.fails} fit failures, {r.flags} flagged, "
                                     f"seed-avg rms {r.dev:.4f}, per-run rms {r.run_dev:.4f}" for r in pcq.itertuples()),
          f"the refits here (FU.fit_acf on every trajectory KOA fitted): {int(PR['refit_ok'].sum())} of {len(PR)} converge; their nu_d against KOA's: "
          f"largest |ratio - 1| {PR['nu_d_vs_koa'].abs().max():.1e}",
          f"stalls (A 0.3 and B 0.7 within 1e-3: the fit ends at its initial guess): {int(PR['stall'].sum())}; misses (rms > 0.15): {int(PR['miss'].sum())}, "
          f"{int((PR['stall'] & PR['miss']).sum())} of them stalls; nu_d / nu_ref of the misses: mean {PR[PR.miss]['nu_d'].div(PR[PR.miss]['nu_ref']).mean():.4f}, "
          f"sd {PR[PR.miss]['nu_d'].div(PR[PR.miss]['nu_ref']).std():.4f}; of the rest: sd {PR[~PR.miss & (PR.refit_ok == 1)]['nu_d'].div(PR[~PR.miss & (PR.refit_ok == 1)]['nu_ref']).std():.4f}",
          f"the damping correction E1 applies, nu_0 / nu_d - 1 = 1 / (8 Q^2): {R['damping_corr'].min():.1e} .. {R['damping_corr'].max():.1e} over block D's 54 (cell, mass)"]
    open(os.path.join(OUT, "d2_table.md"), "w").write("\n".join(md) + "\n"); open(os.path.join(OUT, "d_summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    cells = sorted(R["cell"].unique()); fig, axs = plt.subplots(len(cells), 3, figsize=(17, 2.9 * len(cells)), sharex=True)
    for i, cell in enumerate(cells):
        N = int(cell[3:7])
        for j, al in enumerate((1, 5, 20)):
            ax = axs[i, j]; M = int(round(al * N)); x, c, g = curves[(cell, M)]; rr = R[(R.cell == cell) & (R.M == M)].iloc[0]
            ax.plot(x, c, color="#1f3b5a", lw=0.9, label="seed-averaged ACF (used runs)"); ax.plot(x, g, color="#c0392b", lw=1.0, ls="--", label="item-4a fit")
            ax.plot(x, c - g - 1.25, color="#2a9d8f", lw=0.8, label="ACF - fit, shifted by -1.25"); ax.axhline(-1.25, color="0.6", lw=0.5)
            ax.set_ylim(-1.6, 1.05); ax.grid(alpha=0.25); ax.tick_params(labelsize=7)
            ax.set_title(f"{cell}, alpha {al}: Q_k {rr.Q_k:.0f}, tau_r / window {rr.tau_r_over_window:.2f}, rms {rr.meanfit_rms:.4f}; fails {rr.fit_failures}", fontsize=8)
            if i == 0 and j == 0: ax.legend(fontsize=7, loc="lower right")
            if i == len(cells) - 1: ax.set_xlabel("lag [predicted periods]", fontsize=8)
    fig.suptitle("block D: the seed-averaged position ACF and the registered damped-cosine fit, three masses per cell", fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_d1_block_D_acf_fits.png"), dpi=110, bbox_inches="tight"); plt.close(fig)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15, 5.2))
    for cell in cells:
        q = R[R.cell == cell]; ax.plot(q["alpha"], q["Q_k"], "o-", ms=4, lw=1, label=cell)
        f = q[q["flagged"] == 1]; ax.plot(f["alpha"], f["Q_k"], "kx", ms=10, mew=2)
        bx.plot(q["alpha"], q["run_rms_median"], "o-", ms=4, lw=1, label=cell)
    ax.set_xscale("log"); ax.set_xlabel("alpha = M / N"); ax.set_ylabel("Q_k (seed-averaged fit)"); ax.grid(alpha=0.25); ax.legend(fontsize=7)
    ax.set_title("block D: Q per mass (x: flagged by the group rule)", fontsize=10)
    bx.set_xscale("log"); bx.set_xlabel("alpha = M / N"); bx.set_ylabel("per-run |ACF - fit| rms (median over the runs)"); bx.grid(alpha=0.25)
    bx.set_title("block D: the deviation of the damped cosine from each run's ACF", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig_d2_block_D_Q_and_deviation.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
