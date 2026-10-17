#!/usr/bin/env python3
"""Extra context check: P2 executed events / the package model (upper and lower Z), per cell and pin; P2 per model x pin counts. Read-only."""
import glob, math, os, sys
import numpy as np, pandas as pd
HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
for p in (os.path.join(HS, "cluster/gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import plot_speed_of_sound_edmd as SOS
OUT = os.path.dirname(os.path.abspath(__file__))
J = pd.read_csv(os.path.join(OUT, "my_p2_join.csv"))
R = pd.concat([pd.read_csv(f, dtype={"pin": str}) for f in sorted(glob.glob(os.path.join(HS, "experiments_gen3_stage2_koa_261016/P2/merged/stage2/gen3/P2_*/m_*/p2_runs.csv")))], ignore_index=True)
R = R.merge(J[["cell", "M", "run", "model"]], on=["cell", "M", "run"])
def zl(eta): return float((SOS.Z_liu_global if eta < 0.72 - 1e-12 else SOS.Z_liu_solid)(np.array([eta]))[0])
ev = []
for r in R.to_dict("records"):
    nup = GEN.nu_pred_binary(900, r["H"], r["L0"], r["M"])
    Trec = math.ceil((200.0 / nup) / GEN.DTF) * GEN.DTF
    f = 1.002 if r["pin"] == "turned" else 1.171
    Z = zl(r["eta"])
    ev.append(dict(nup_mine=nup, nup_csv=r["nu_pred"], low=(1.2004 * (Z - 1) + 0.0112) * 900 * (5000 + Trec) * f,
                   up=(1.2004 * (Z * 1.268 - 1) + 0.0112) * 900 * (5000 + Trec) * f))
E = pd.DataFrame(ev); R = pd.concat([R, E], axis=1)
print("max rel |nu_pred_binary(H, L0 from csv) - csv nu_pred|:", float(((R.nup_mine - R.nup_csv).abs() / R.nup_csv).max()))
R["r_up"] = R.events / R.up; R["r_low"] = R.events / R.low
print(R.groupby(["eta", "G", "pin"]).agg(n=("events", "size"), r_up_min=("r_up", "min"), r_up_med=("r_up", "median"), r_up_max=("r_up", "max"),
                                         r_low_med=("r_low", "median")).round(3).to_string())
print("\nturned cells: executed / upper, range over runs:", round(R[R.pin == "turned"].r_up.min(), 3), "..", round(R[R.pin == "turned"].r_up.max(), 3))
print("turned cells: per-cell medians of executed / upper:", R[R.pin == "turned"].groupby(["eta", "G"]).r_up.median().round(3).to_dict())
print("\nP2 runs per model x pin:"); print(pd.crosstab(R.model, R.pin).to_string())
print("\nP2 runs per model x eta:"); print(pd.crosstab(R.model, R.eta).to_string())
