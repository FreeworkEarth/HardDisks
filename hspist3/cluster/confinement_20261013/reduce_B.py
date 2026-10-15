#!/usr/bin/env python3
"""##CHRIS 2026-10-13: method B reduction for one cell directory -> per m_<M>: red_nu.csv and acf_runs.npz.
red_nu.csv: per trajectory, the canonical estimator's frequency, computed with the SAME functions as
paper1_populate_cs_err_20261002.cell (TD = 200, X_EDGE = 2.5, T._load/_prefix/_spectrum, argmax above k_min),
so the Mac analysis can use the summaries instead of the traces. acf_runs.npz: the mean-removed position ACF of
each trajectory to 20 predicted periods (for Gamma = 2/tau_r and P_1 = B, methods sec. 13).
GATE before use: on the full pilot cell copied back, red_nu.csv must equal cell()'s per-run values exactly.
usage (from hspist3/): python3 cluster/confinement_20261013/reduce_B.py <cell_dir> [M ...]
##CHRIS 2026-10-07 (gate v2, Test T): optional masses after the cell directory -- only those m_<M> directories are reduced, so
array tasks that run one mass each never write the same red_nu.csv."""
import glob, os, re, sys
import numpy as np, pandas as pd
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
from paper1_populate_cs_err_20261002 import TD, X_EDGE
cell = sys.argv[1]
only = {f"m_{m}" for m in sys.argv[2:]}
for d in sorted(glob.glob(os.path.join(cell, "m_*"))):
    if only and os.path.basename(d) not in only: continue
    rows, acfs = [], {}
    for p in sorted(glob.glob(os.path.join(d, "wall_x_positions_L0_*_run*.csv"))):
        r = int(re.search(r"_run(\d+)\.csv$", p).group(1))
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD); P, df = T._spectrum(x[:n], dt); k = int(round(TD / X_EDGE))
        nu = (k + int(np.argmax(P[k:]))) * df
        seed = int(pd.read_csv(p, nrows=1)["Seed"].iloc[0])
        y = x[:n] - x[:n].mean(); nl = int(20 / (nup * dt))
        f = np.fft.rfft(y, 2 * len(y)); a = np.fft.irfft(f * np.conj(f))[:nl + 1] / np.arange(len(y), len(y) - nl - 1, -1)
        acfs[f"run{r}"] = a.astype(np.float32)
        rows.append(dict(run=r, seed=seed, nu=nu, nu_pred=nup, dt=dt, n=n))
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(d, "red_nu.csv"), index=False)
        np.savez_compressed(os.path.join(d, "acf_runs.npz"), **acfs)
        print(f"{d}: {len(rows)} trajectories reduced")
