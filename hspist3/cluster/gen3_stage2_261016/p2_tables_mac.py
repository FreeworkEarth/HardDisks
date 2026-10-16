#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE P2 TABLES ON THE MAC, from the fetched tree
(hspist3/experiments_gen3_stage2_koa_261016/P2; fetch_p2_koa.sh). Reads only the merged tables; writes into a NEW folder --out.
  per (cell, mass)  E1 exactly as s1_verdict_261012.per_mass: the used trajectories (gen3 clean and fit_ok), nu_d-bar (nud_runs.csv),
                    tau_r of the item-4a fit of their seed-averaged ACF (the ACFs cut to the shortest; nu_ref the folder's argmax mean),
                    nu_0 = sqrt(nu_d-bar^2 + (1 / (2 pi tau_r))^2), Q = pi nu_d-bar tau_r, the delete-one jackknife SE (3 or more
                    trajectories), c = nu_0 / x, x = K(alpha) / (2 pi (L0 - 1.025)); the low-Q mark (Q < 5); with fewer than 3 used
                    trajectories nu_0 and c are printed without an SE and marked "n < 3" (stage 1 would give no value);
  per cell          E1 = the inverse-variance mean of the masses' c (stage 1's e1: 5 or more masses with an SE), chi2 and dof, the low-Q flag;
  per trajectory    (p2_runs.csv, p2_reduce.py on KOA) nu_d, its own c = sqrt(nu_d^2 + (1 / (2 pi tau_r,run))^2) / x, Q and tau_r of its
                    own fit, the stall flag, the spectrum share, sd / sd_E1 (sd_E1 = sqrt(L_eff / (2 H rho c^2)), c the (cell, mass) E1,
                    H the cell's height, rho = 4 eta / pi), the drift, Z, psi6 per compartment and the orientation's deviation from the
                    pin (deg, in (-30, 30]: 0 for pin 0, 30 for turned), the strain e of the
                    start's lattice, of the initial state, of the hold's last snapshot and of the last snapshot, per compartment.
Writes p2_masses.csv, p2_cells.csv, p2_seeds.csv, p2_tables.md, fig_p2_c_per_mass.png, fig_p2_strain.png.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/p2_tables_mac.py --data <fetched P2 folder> --out <new folder>
"""
import glob, math, os, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import tests_20260913 as T
import resched_testT_followup_261007 as FU
warnings.filterwarnings("ignore")


def seedavg(acfs, dt, nu_ref):
    nl = min(len(a) for a in acfs); m = np.mean([np.asarray(a[:nl], float) for a in acfs], axis=0)
    try: return FU.fit_acf(m, dt, nu_ref)
    except Exception: return None


def per_mass(d):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run"); nud = pd.read_csv(os.path.join(d, "nud_runs.csv"))
    hl = pd.read_csv(os.path.join(d, "runs_health.csv")); z = np.load(os.path.join(d, "acf_runs.npz"))
    clean = set(hl[(hl["gen3_clean"] == 1) & (hl["gen3_records"] == 1)]["run"]); ok = set(nud[nud["fit_ok"] == 1]["run"])
    used = sorted(set(red["run"]) & clean & ok); dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    out = dict(present=len(red), used=len(used), fit_failures=len((set(red["run"]) & clean) - ok))
    if not used: return out
    nd = nud.set_index("run").loc[used, "nu_d"].to_numpy(float); acfs = [z[f"run{r}"] for r in used]; p = seedavg(acfs, dt, nu_ref)
    out.update(nud=float(nd.mean()), nud_sd=float(nd.std(ddof=1)) if len(nd) > 1 else math.nan)
    if p is None: return out
    nu0 = math.sqrt(out["nud"] ** 2 + (1 / (2 * math.pi * p[3])) ** 2); jk = []
    if len(used) >= 3:
        for i in range(len(used)):
            keep = [k for k in range(len(used)) if k != i]; fi = seedavg([acfs[k] for k in keep], dt, nu_ref)
            if fi is not None: jk.append(math.sqrt(nd[keep].mean() ** 2 + (1 / (2 * math.pi * fi[3])) ** 2))
    n = len(jk); se = math.sqrt((n - 1) / n * sum((v - np.mean(jk)) ** 2 for v in jk)) if n >= 3 else math.nan
    out.update(tau_r=p[3], nu0=nu0, nu0_se=se, Q=math.pi * out["nud"] * p[3], A=p[0], B=p[2])
    return out


def main():
    av = sys.argv; data = os.path.abspath(av[av.index("--data") + 1]); out = os.path.abspath(av[av.index("--out") + 1])
    if os.path.exists(out): sys.exit(f"REFUSED: {out} exists (a new folder only)")
    folders = sorted(glob.glob(os.path.join(data, "merged", "stage2", "gen3", "*", "m_*")))
    if not folders: sys.exit(f"STOP: no merged/stage2/gen3/<cell>/m_<M> under {data}")
    os.makedirs(out); mrows, srows = [], []
    for d in folders:
        cid = d.split(os.sep)[-2]; M = int(d.split("_")[-1]); f = cid.split("_"); eta = float(f[2][3:]); G = f[3]; pin = f[4][3:]
        pr = pd.read_csv(os.path.join(d, "p2_runs.csv")); H = float(pr["H"].iloc[0]); L0 = float(pr["L0"].iloc[0]); a = M / 900.0
        x = T.k_root(a) / (2 * math.pi * (L0 - 1.025)); pm = per_mass(d)
        c = pm.get("nu0", math.nan) / x; sc = pm.get("nu0_se", math.nan) / x
        mrows.append(dict(cell=cid, eta=eta, G=G, pin=pin, H=H, L0=L0, alpha=a, M=M, **pm, x=x, c=c, sig_c=sc,
                          lowQ=int(pm.get("Q", math.inf) < 5), note="" if pm.get("used", 0) >= 3 else "n < 3: no SE, not an E1 of record"))
        rho = 4 * eta / math.pi; leff = L0 - 1.025
        for _, r in pr.iterrows():
            cs = math.sqrt(r.get("nu_d", math.nan) ** 2 + (1 / (2 * math.pi * r.get("fit_tau_r", math.nan))) ** 2) / x if r.get("fit_tau_r", 0) else math.nan
            sde1 = math.sqrt(leff / (2 * H * rho * c ** 2)) if math.isfinite(c) and c > 0 else math.nan
            row = {k: r.get(k, math.nan) for k in ("run", "seed", "node", "seconds", "events", "clean", "nu_pred", "nu_arg", "nu_d", "fit_ok", "fit_tau_r",
                                                  "fit_Q", "stall", "fit_rms", "spectrum_share", "welch_peak_over_nup", "x_mean", "x_sd", "drift_per_1000",
                                                  "slow_sd", "fast_sd", "Z", "psi6_rec_0", "psi6_rec_1", "orient_rec_0", "orient_rec_1", "arrangement",
                                                  "n_vac", "e_lattice", "start_e_0", "start_e_1", "start_div_minus_L0", "hold_e_0", "hold_e_1",
                                                  "end_e_0", "end_e_1", "end_theta0_0", "end_theta0_1", "end_ordered_0", "end_ordered_1",
                                                  "end_div_minus_L0", "m4c", "error")}
            row.update(cell=cid, eta=eta, G=G, pin=pin, H=H, L0=L0, alpha=a, M=M, c_seed=cs, sd_E1=sde1, sd_over_E1=r.get("x_sd", math.nan) / sde1)
            pa = 0.0 if pin == "0" else 30.0                       # the pinned orientation; deviations in (-30, 30]
            for k in (0, 1):
                for col in (f"orient_rec_{k}", f"end_theta0_{k}"):
                    v = row.get(col, math.nan); row[col + "_dev"] = ((v - pa + 30.0) % 60.0) - 30.0 if isinstance(v, (int, float)) and math.isfinite(v) else math.nan
            srows.append(row)
    Mt = pd.DataFrame(mrows); S = pd.DataFrame(srows); crows = []
    for cid, g in Mt.groupby("cell", sort=False):
        u = g[(g["alpha"] >= 1) & np.isfinite(g["c"]) & np.isfinite(g["sig_c"]) & (g["sig_c"] > 0)]
        if len(u) >= 5:
            w = 1 / u["sig_c"] ** 2; cs = float((w * u["c"]).sum() / w.sum()); ch = float((w * (u["c"] - cs) ** 2).sum()); dof = len(u) - 1
            sig = float(w.sum() ** -0.5 * max(1.0, math.sqrt(ch / dof)))
        else: cs = sig = ch = math.nan; dof = len(u) - 1
        sg = S[S["cell"] == cid]
        crows.append(dict(cell=cid, eta=g["eta"].iloc[0], G=g["G"].iloc[0], pin=g["pin"].iloc[0], H=g["H"].iloc[0], L0=g["L0"].iloc[0],
                          masses=len(u), E1_c=cs, E1_sig=sig, chi2=ch, dof=dof, lowQ=int((g["Q"] < 5).any()) if "Q" in g else 0,
                          e_lattice=float(sg["e_lattice"].mean()), e_start=float(sg[["start_e_0", "start_e_1"]].stack().mean()),
                          e_hold=float(sg[["hold_e_0", "hold_e_1"]].stack().mean()), e_end=float(sg[["end_e_0", "end_e_1"]].stack().mean()),
                          psi6_rec=float(sg[["psi6_rec_0", "psi6_rec_1"]].stack().mean()), Z=float(sg["Z"].mean()), runs=len(sg),
                          stalls=int(sg["stall"].sum()), errors=int((sg["error"].fillna("") != "").sum())))
    C = pd.DataFrame(crows)
    for df, n in ((Mt, "p2_masses.csv"), (C, "p2_cells.csv"), (S, "p2_seeds.csv")): df.to_csv(os.path.join(out, n), index=False, float_format="%.8g")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    md = ["# P2 tables (p2_tables_mac.py)", "", f"data: {data}", "", "## per cell", "",
          C.to_string(index=False, float_format=lambda v: f"{v:.5g}"), "", "## per (cell, mass)", "",
          Mt[["cell", "alpha", "present", "used", "fit_failures", "nud", "tau_r", "Q", "nu0", "c", "sig_c", "lowQ", "note"]].to_string(index=False, float_format=lambda v: f"{v:.5g}"),
          "", "## per trajectory", "",
          S[["cell", "alpha", "run", "seed", "clean", "nu_d", "c_seed", "fit_Q", "stall", "spectrum_share", "sd_over_E1", "drift_per_1000", "Z",
             "psi6_rec_0", "psi6_rec_1", "orient_rec_0_dev", "orient_rec_1_dev", "e_lattice", "start_e_0", "start_e_1", "hold_e_0", "hold_e_1", "end_e_0",
             "end_e_1", "start_div_minus_L0", "end_div_minus_L0", "error"]].to_string(index=False, float_format=lambda v: f"{v:.4g}")]
    open(os.path.join(out, "p2_tables.md"), "w").write("\n".join(md) + "\n"); print("\n".join(md))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for cid, g in Mt.groupby("cell", sort=False):
        ax.errorbar(g["alpha"], g["c"], yerr=g["sig_c"].fillna(0), marker="o", ms=5, capsize=3, lw=1, label=cid.replace("P2_N0900_", ""))
    ax.set_xscale("log"); ax.set_xlabel("alpha = M / N"); ax.set_ylabel("c (E1 per mass)"); ax.grid(alpha=0.25); ax.legend(fontsize=7, ncol=2)
    ax.set_title("P2: the E1 speed per mass and cell", fontsize=10); fig.savefig(os.path.join(out, "fig_p2_c_per_mass.png"), dpi=120, bbox_inches="tight"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 5.5)); lab = {"start": "o", "hold": "s", "end": "D"}
    for i, (cid, g) in enumerate(S.groupby("cell", sort=False)):
        for k, mk in lab.items():
            v = 100 * g[[f"{k}_e_0", f"{k}_e_1"]].to_numpy(float).ravel(); ax.plot(np.full(len(v), i) + {"start": -0.2, "hold": 0, "end": 0.2}[k], v, mk,
                                                                                ms=4, alpha=0.6, label=k if i == 0 else None)
    ax.axhline(0, color="0.5", lw=0.8); ax.set_xticks(range(S["cell"].nunique())); ax.set_xticklabels([c.replace("P2_N0900_", "") for c in S["cell"].unique()], rotation=30, fontsize=7, ha="right")
    ax.set_ylabel("strain e [%] per compartment"); ax.grid(alpha=0.25); ax.legend(fontsize=8); ax.set_title("P2: the strain at the start, at the hold's end and at the end", fontsize=10)
    fig.savefig(os.path.join(out, "fig_p2_strain.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
