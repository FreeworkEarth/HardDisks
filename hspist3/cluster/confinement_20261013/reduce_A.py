#!/usr/bin/env python3
"""##CHRIS 2026-10-13: method A reduction, one seed -> one row (261012 sec. 1.4). Fails closed.
dp in the event log is the PARTICLE's momentum change (edmd.c:1069). For a held divider a particle arriving from
the left leaves with dp < 0, from the right with dp > 0: F_L = -sum_{dp<0} dp / T_w, F_R = sum_{dp>0} dp / T_w.
Window: release (event t = 200 sigma, after the 12000-step hold) to the last event. T_L, T_R from the trace's
KE_gas_left/right over the same window (2D: KE = N_side kT).
usage: reduce_A.py ev_<seed>.csv tr_<seed>.csv red_<seed>.csv"""
import sys
import numpy as np, pandas as pd
ev, tr, out = sys.argv[1:4]
HOLD = 200.0
e = pd.read_csv(ev, usecols=["t_sigma", "kind", "dp"]); d = e[(e["kind"] == "D0") & (e["t_sigma"] >= HOLD)]
if (d["dp"] == 0).any(): sys.exit(f"ABORT {ev}: {int((d['dp'] == 0).sum())} divider events with dp == 0 (face undefined)")
t0, t1 = HOLD, float(e["t_sigma"].max()); Tw = t1 - t0
dp = d["dp"].to_numpy()
FL = -dp[dp < 0].sum() / Tw; FR = dp[dp > 0].sum() / Tw
s = pd.read_csv(tr, usecols=["Time", "KE_gas_left", "KE_gas_right", "SegCounts"], low_memory=False)
cnt = s["SegCounts"].astype(str).str.split(";").iloc[0]; NL, NR = int(cnt[0]), int(cnt[-1])
TL = s["KE_gas_left"].mean() / NL; TR = s["KE_gas_right"].mean() / NR
pd.DataFrame([dict(window=Tw, F_L=FL, F_R=FR, n_L=int((dp < 0).sum()), n_R=int((dp > 0).sum()), T_L=TL, T_R=TR,
                   N_L=NL, N_R=NR)]).to_csv(out, index=False)
