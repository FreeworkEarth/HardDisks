#!/usr/bin/env python3
"""Independent verification of Report D57-3 task 4 (a)-(c): the reorientation at eta 0.72 (read-only on the repo)."""
import glob, os, sys
import numpy as np, pandas as pd

P2 = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_gen3_stage2_koa_261016/P2/merged/stage2/gen3"
CC = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/experiments_gen3_d57_261016/reorientation/reorientation_p2.csv"
OUT = os.path.dirname(os.path.abspath(__file__))
CELLS = ["P2_N0900_eta0.720000_G1_pin0", "P2_N0900_eta0.720000_G2_pin0"]
MASSES = [900, 1800, 4500, 9000, 18000]

rows = []
for cell in CELLS:
    for M in MASSES:
        d = os.path.join(P2, cell, f"m_{M}")
        fr = pd.read_csv(os.path.join(d, "p2_frames.csv"))
        ru = pd.read_csv(os.path.join(d, "p2_runs.csv"))
        assert sorted(ru["run"].tolist()) == list(range(10)), (cell, M)
        for run in range(10):
            f = fr[fr["run"] == run].sort_values("t").reset_index(drop=True)
            tt = f["t"].to_numpy(float)
            pf = glob.glob(os.path.join(d, f"psi6_t_L0_*_run{run}.csv.gz")); assert len(pf) == 1, pf
            ps = pd.read_csv(pf[0]); ts = ps["t_sigma"].to_numpy(float)
            rr = ru[ru["run"] == run].iloc[0]
            t_first_pos = tt[tt > 0].min()
            for k in (0, 1):
                th = f[f"theta0_{k}"].to_numpy(float)
                m60 = np.mod(th, 60.0)
                isy = np.abs(m60 - 30.0) <= 15.0           # NaN -> False -> 'x' (no y)
                isy_strict = np.abs(m60 - 30.0) < 15.0
                fy = tt[isy][0] if isy.any() else np.nan
                fy_s = tt[isy_strict][0] if isy_strict.any() else np.nan
                o = np.mod(ps[f"psi6_orient_deg_{k}"].to_numpy(float) + 15.0, 60.0) - 15.0
                near = np.abs(o - 30.0) <= 7.5
                near_s = np.abs(o - 30.0) < 7.5
                fp = ts[near][0] if near.any() else np.nan
                fp_s = ts[near_s][0] if near_s.any() else np.nan
                rows.append(dict(cell=cell, M=M, run=run, side=k, first_y=fy, first_y_strict=fy_s, first_pos_snap=t_first_pos,
                                 at_first_snap=bool(fy == t_first_pos), snap_dt_median=float(np.median(np.diff(tt))),
                                 n_nan_theta=int(np.isnan(th).sum()), first_p=fp, first_p_strict=fp_s,
                                 psi_dt_median=float(np.median(np.diff(ts))),
                                 start_line=rr[f"start_line_{k}"], hold_line=rr[f"hold_line_{k}"], end_line=rr[f"end_line_{k}"],
                                 n_frames=len(f), n_psi=len(ps), t_last=tt[-1]))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "v4_abc_per_compartment.csv"), index=False, float_format="%.9g")

print("compartments per cell:", D.groupby("cell").size().to_dict())
print("NaN theta0 frames total:", D["n_nan_theta"].sum())
for cell in CELLS:
    g = D[D["cell"] == cell]; fy = g["first_y"]; fp = g["first_p"]
    print(f"\n== {cell}")
    print(f"(a) first y snapshot: missing {fy.isna().sum()}; median {fy.median():.4f}, range {fy.min():.4f} to {fy.max():.4f}; "
          f"<200 {(fy < 200).sum()}, <300 {(fy < 300).sum()}, <500 {(fy < 500).sum()}; at the first snapshot after 0: {g['at_first_snap'].sum()}")
    print(f"    strict (<15) variant identical: {np.array_equal(g['first_y'].to_numpy(), g['first_y_strict'].to_numpy(), equal_nan=True)}")
    print(f"    snapshot interval median per (M): {g.groupby('M')['snap_dt_median'].median().round(4).to_dict()}")
    print(f"    CC's '<= 1.5 interval' rule: {(g['first_y'] <= 1.5 * g['snap_dt_median']).sum()}")
    print(f"(b) first psi6_t within 7.5 of 30: missing {fp.isna().sum()}; median {fp.median():.4f}, range {fp.min():.4f} to {fp.max():.4f}; "
          f"<200 {(fp < 200).sum()}, <300 {(fp < 300).sum()}, <500 {(fp < 500).sum()}")
    print(f"    strict (<7.5) variant identical: {np.array_equal(g['first_p'].to_numpy(), g['first_p_strict'].to_numpy(), equal_nan=True)}")
    print(f"    psi6_t interval range: {g['psi_dt_median'].min():.4f} to {g['psi_dt_median'].max():.4f}")
    print(f"(c) start_line counts {g['start_line'].value_counts().to_dict()}, hold_line {g['hold_line'].value_counts().to_dict()}, "
          f"end_line {g['end_line'].value_counts().to_dict()}")
print(f"\n(c) all 200: start 0: {(D['start_line'] == 0).all()}, end 90: {(D['end_line'] == 90).all()}, hold 90: {(D['hold_line'] == 90).sum()} of {len(D)}")

# row-by-row against CC's CSV
C = pd.read_csv(CC)
J = D.merge(C, on=["cell", "M", "run", "side"], how="outer", indicator=True)
print("\nmerge with CC's reorientation_p2.csv:", J["_merge"].value_counts().to_dict())
for mine, cc in (("first_y", "first_y_snapshot"), ("first_p", "first_psi6t_near30"), ("snap_dt_median", "snapshot_every"), ("psi_dt_median", "psi6t_every")):
    a, b = J[mine].to_numpy(float), J[cc].to_numpy(float)
    rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-12)
    print(f"  {mine} vs CC {cc}: NaN mismatch {int((np.isnan(a) != np.isnan(b)).sum())}, max rel diff {np.nanmax(rel):.3e}")
print("  start_line equal:", (J["start_line_x"] == J["start_line_y"]).all(), " end_line equal:", (J["end_line_x"] == J["end_line_y"]).all())
