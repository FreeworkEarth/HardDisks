#!/usr/bin/env python3
"""##CHRIS 2026-10-06 (261012 sec. 4.4.7): diagnostic of the FAILED KOA smoke test of build 73fc07f (minimal rescheduling,
branch engine-divider-resched). Analysis only, no simulation.
Compares the per-mass nu of three one-seed-per-mass pi/8 pilots (Mac 05215ea, KOA 279282b job 14966575, KOA 73fc07f) with the
25-seed distribution of the same cell (epi8_H_H10_L10, 279282b confinement campaign), and bootstraps the one-seed-per-mass
through-origin c_s from those 25 seeds per mass (200000 draws, seed 20261006).
usage (from hspist3/): python3 validation/resched_smoke_diag_261006.py"""
import contextlib, glob, io, math, os, sys
import numpy as np, pandas as pd
HS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(HS, "validation"), HS, os.path.join(HS, "cluster")]
import tests_20260913 as T
import confinement_pilot as CP
from paper1_populate_cs_err_20261002 import cell
L0 = 10.0
P = os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_pilot_20261013")
B = os.path.join(HS, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/confinement_B_20261013/epi8_H_H10_L10")
NEW = {50: 0.07040729, 100: 0.05624957, 200: 0.04478494, 300: 0.03664223, 500: 0.02967332, 750: 0.02414966,
       1000: 0.02112673, 1500: 0.01718229, 2000: 0.01506197}   # [DATA] the per-mass nu table of the KOA smoke-test log
# (conf-smoke, node cn-03-33-01, build 73fc07f, pasted by Chris 2026-10-06); to be replaced by the fetched pilot directory
def pilot(d):
    out = {}
    for M in T.A1_MASSES:
        with contextlib.redirect_stdout(io.StringIO()):
            out[M] = cell((CP.ETA, L0, M, T.cell_runs(os.path.join(d, f"m_{M}"), M)))["nu"]
    return out
mac, koa_old = pilot(os.path.join(P, "mac_pi8_H10_L10")), pilot(os.path.join(P, "koa_pi8_H10_L10"))
x = {M: T.x_of(M, L0) for M in T.A1_MASSES}
camp = {M: pd.read_csv(os.path.join(B, f"m_{M}", "red_nu.csv"))["nu"].to_numpy(float) for M in T.A1_MASSES}
def cs(nu): xs = np.array([x[M] for M in T.A1_MASSES]); y = np.array([nu[M] for M in T.A1_MASSES]); return float((xs * y).sum() / (xs * xs).sum())
print("| M | alpha | weight x^2 share | campaign mean nu (25 seeds) | campaign sd | Mac pilot (z) | KOA 279282b pilot (z) | KOA 73fc07f pilot (z) |")
print("|---|---|---|---|---|---|---|---|")
w = {M: x[M] ** 2 for M in T.A1_MASSES}; W = sum(w.values())
for M in T.A1_MASSES:
    m, s = camp[M].mean(), camp[M].std(ddof=1)
    f = lambda v: f"{v:.6f} ({(v - m) / s:+.2f})"
    print(f"| {M} | {M / 100:g} | {w[M] / W:.3f} | {m:.6f} | {s:.6f} | {f(mac[M])} | {f(koa_old[M])} | {f(NEW[M])} |")
c_camp = cs({M: camp[M].mean() for M in T.A1_MASSES})
print(f"\nthrough-origin c_s: campaign means {c_camp:.5f}; Mac pilot {cs(mac):.5f}; KOA 279282b pilot {cs(koa_old):.5f}; KOA 73fc07f pilot {cs(NEW):.5f}")
rng = np.random.default_rng(20261006); draws = np.array([cs({M: rng.choice(camp[M]) for M in T.A1_MASSES}) for _ in range(200000)])
print(f"bootstrap of a one-seed-per-mass pilot from the 25 campaign seeds (200000 draws): mean {draws.mean():.5f}, sd {draws.std():.5f}")
for lab, v in (("Mac pilot", cs(mac)), ("KOA 279282b pilot", cs(koa_old)), ("KOA 73fc07f pilot", cs(NEW))):
    print(f"  {lab}: {v:.5f}  z = {(v - draws.mean()) / draws.std():+.2f}  fraction of draws <= it: {(draws <= v).mean():.4f}")
d = draws - rng.permutation(draws)
print(f"difference of two independent pilots: sd {d.std():.5f}; P(|diff| > 0.06903) = {(np.abs(d) > 0.06903).mean():.4f}; "
      f"P(diff <= -0.11462) = {(d <= -0.11462).mean():.5f}")
z50 = (NEW[50] - camp[50].mean()) / camp[50].std(ddof=1); z100 = (NEW[100] - camp[100].mean()) / camp[100].std(ddof=1)
print(f"light masses of the 73fc07f pilot: z(M=50) = {z50:+.2f}, z(M=100) = {z100:+.2f} against the 25-seed spread")
