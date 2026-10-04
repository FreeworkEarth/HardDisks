#!/usr/bin/env python3
"""##CHRIS 2026-10-04 (Task Y, 261012 sec. 3 "A-fixed"): reduction of one held-divider seed -> one row. Fails closed.
As reduce_A.py (dp of the divider events D0 is the PARTICLE's momentum change, edmd.c:1069; the face is the sign of dp:
F_L = -sum_{dp<0} dp / T_w, F_R = sum_{dp>0} dp / T_w), but over the HELD window [t0, t1) only: t0 = 200 (end of the
equilibration), t1 = hold * dt = the release. Extra columns for the sec. 3 gates: u_wall_max = max |divider velocity| over
the window's D0 events (0 if the divider is immovable), W_div = sum of dE over them (0: no work on either gas),
t_last = last event time (must be >= t1). T_L, T_R from the trace, which is written after the release only
(00ALLINONE.c:17069): with the divider held, each compartment is closed and its temperature cannot change.
usage: reduce_AF.py ev_<seed>.csv tr_<seed>.csv red_<seed>.csv <t0> <t1>"""
import sys
import numpy as np, pandas as pd
ev, tr, out, t0, t1 = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
e = pd.read_csv(ev, usecols=["t_sigma", "kind", "u_wall", "dE", "dp"])
t_last = float(e["t_sigma"].max())
if t_last < t1: sys.exit(f"ABORT {ev}: the event log ends at t = {t_last} before the release t1 = {t1}")
d = e[(e["kind"] == "D0") & (e["t_sigma"] >= t0) & (e["t_sigma"] < t1)]
if (d["dp"] == 0).any(): sys.exit(f"ABORT {ev}: {int((d['dp'] == 0).sum())} divider events with dp == 0 (face undefined)")
Tw = t1 - t0; dp = d["dp"].to_numpy()
FL = -dp[dp < 0].sum() / Tw; FR = dp[dp > 0].sum() / Tw
s = pd.read_csv(tr, usecols=["Time", "KE_gas_left", "KE_gas_right", "SegCounts"], low_memory=False)
cnt = s["SegCounts"].astype(str).str.split(";").iloc[0]; NL, NR = int(cnt[0]), int(cnt[-1])
TL = s["KE_gas_left"].mean() / NL; TR = s["KE_gas_right"].mean() / NR
pd.DataFrame([dict(window=Tw, F_L=FL, F_R=FR, n_L=int((dp < 0).sum()), n_R=int((dp > 0).sum()), T_L=TL, T_R=TR,
                   N_L=NL, N_R=NR, u_wall_max=float(np.abs(d["u_wall"]).max()) if len(d) else 0.0,
                   W_div=float(d["dE"].sum()), t_last=t_last)]).to_csv(out, index=False)
