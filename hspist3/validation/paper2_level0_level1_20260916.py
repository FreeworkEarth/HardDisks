#!/usr/bin/env python3
"""##CHRIS 2026-09-16: Paper 2, Levels 0 and 1, from existing data plus 40 short hold runs. Analysis only.

Level 0 -- ledgers of the 2026-09-10 pilot (geometry A, eta = 0.10, N_s = 50, u = 0.02/0.05/0.10, 25 seeds each),
from the per-event log (HD_PISTON_EVENTS: every outer-wall, divider and piston collision) against the per-sample
trace (KE_gas_total, Px_gas, PistonWork). Two analysis corrections, both bookkeeping, neither physics:
  (1) sample times are rebuilt from the step index -- the trace prints Time to 6 decimals, so an event within
      1e-6 sigma of a sample boundary was assigned to the wrong sample (one sample in 17000, |R| = |dp| of one event);
  (2) divider events log dE as the DIVIDER's energy gain, piston events as the GAS's; the gas side of a divider
      event is therefore -dE.
Level 1 -- quasi-static work of the finite box. Z_wall measured in equilibrium at three points of the compression
path (L = 39.25, 37.291667, 35.3125; the last two are grid-exact, 2 L x 24 integer, or the validator aborts),
from outer-wall (WR) and gas-side divider (D0, dp > 0) impulses during the 200-sigma hold. W_qs = N_s kT
(exp INT Z_wall/eta deta - 1), with Z_wall/Z_KR interpolated linearly along the path.
"""
import glob, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos

ET = os.path.join(os.path.dirname(HERE), "experiments_energy_transfer")
DT = 0.0166666669


def ledgers(d, hold_steps=1000):
    hold = hold_steps * DT; worst = dict(raw=0.0, fix=0.0, p=0.0); n = 0
    for tr in sorted(glob.glob(os.path.join(d, "tr_*.csv"))):
        sd = os.path.basename(tr)[3:-4]
        t = pd.read_csv(tr, usecols=["Time", "KE_gas_total", "PistonWork", "Px_gas"], low_memory=False)
        e = pd.read_csv(os.path.join(d, f"ev_{sd}.csv"))
        et = e["t_sigma"].to_numpy(float); o = np.argsort(et, kind="stable"); et = et[o]
        kinds = e["kind"].astype(str).to_numpy()[o]; dE = e["dE"].to_numpy(float)[o]; dp = e["dp"].to_numpy(float)[o]
        isD = np.array([k.startswith("D") for k in kinds]); xk = isD | np.isin(kinds, ["WL", "WR", "PL", "PR"])
        step = np.round(np.diff(t["Time"].to_numpy(float)).mean() / DT)
        ts = hold + np.arange(len(t)) * step * DT + t["Time"].iloc[0]
        idx = np.searchsorted(et, ts, side="right")
        cs = lambda x: (lambda c: c[idx] - c[idx[0]])(np.concatenate([[0.0], np.cumsum(x)]))
        ke = t["KE_gas_total"].to_numpy(float); W = abs(t["PistonWork"].iloc[-1] - t["PistonWork"].iloc[0])
        px = t["Px_gas"].to_numpy(float)
        worst["raw"] = max(worst["raw"], np.max(np.abs((ke - ke[0]) - cs(dE))) / W)
        worst["fix"] = max(worst["fix"], np.max(np.abs((ke - ke[0]) - cs(np.where(isD, -dE, dE)))) / W)
        worst["p"] = max(worst["p"], np.max(np.abs((px - px[0]) - cs(np.where(xk, dp, 0.0)))) / np.abs(np.where(xk, dp, 0.0)).sum())
        n += 1
    return n, worst


def zwall(files, L, hold_steps=12000, Ns=50):
    T = hold_steps * DT; z = []
    for ev in files:
        e = pd.read_csv(ev); h = e[e["t_sigma"] < T]
        zr = h[h["kind"] == "WR"]["dp"].abs().sum() * L / (T * Ns)
        zd = h[(h["kind"] == "D0") & (h["dp"] > 0)]["dp"].sum() * L / (T * Ns)
        z.append(0.5 * (zr + zd))
    z = np.array(z); return z.mean(), z.std(ddof=1) / math.sqrt(len(z)), len(z)


if __name__ == "__main__":
    print("## Level 0 -- ledgers of the 2026-09-10 pilot\n")
    print("| speed | seeds | max|R_E|/W as logged | max|R_E|/W, divider sign consistent | max|R_px|/Σ|dp_x| |")
    print("|---|---|---|---|---|")
    for u in ("u0.02", "u0.05", "u0.10"):
        n, w = ledgers(os.path.join(ET, "level0_pilot_20260910", u))
        print(f"| {u[1:]} | {n} | {w['raw']:.2e} | {w['fix']:.2e} | {w['p']:.2e} |")

    print("\n## Level 1 -- wall pressure of the finite box along the compression path\n")
    Ns = 50
    pts = [(39.25, *zwall(sorted(glob.glob(os.path.join(ET, "level0_Wqs_20260911", "u*", "ev_*.csv"))), 39.25))]
    for L, tag in ((37.291667, "L37p291667"), (35.3125, "L35p3125")):
        pts.append((L, *zwall(sorted(glob.glob(os.path.join(ET, "level1_Zwall_path_20260916", tag, "ev_*.csv"))), L)))
    print("| compartment length | η | seeds | Z_wall | sem | Z_KR | Z_wall/Z_KR |")
    print("|---|---|---|---|---|---|---|")
    rows = []
    for L, z, s, n in pts:
        eta = 3.9269908169872414 / L; kr = float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0])
        rows.append((eta, z, s, kr)); print(f"| {L} | {eta:.6f} | {n} | {z:.4f} | {s:.4f} | {kr:.4f} | {z/kr:.4f} |")
    eta = np.array([r[0] for r in rows]); ratio = np.array([r[1] / r[3] for r in rows]); rsem = np.array([r[2] / r[3] for r in rows])
    e = np.linspace(eta.min(), eta.max(), 4001); Zk = sos.Z_kolafa_rottner_2006(e)
    W_bulk = Ns * (math.exp(np.trapezoid(Zk / e, e)) - 1)
    W_const = Ns * (math.exp(np.trapezoid(ratio[0] * Zk / e, e)) - 1)
    W_meas = Ns * (math.exp(np.trapezoid(np.interp(e, eta, ratio) * Zk / e, e)) - 1)
    rng = np.random.default_rng(1)
    dW = float(np.std([Ns * (math.exp(np.trapezoid(np.interp(e, eta, ratio + rng.normal(0, rsem)) * Zk / e, e)) - 1) for _ in range(4000)]))
    print(f"\npath η {eta.min():.6f} → {eta.max():.6f} (the push ended at 0.111183; the grid-exact hold is 0.02 % further)")
    print(f"W_qs bulk KR adiabat:                           {W_bulk:.4f}")
    print(f"W_qs finite box, start ratio held constant:     {W_const:.4f}")
    print(f"W_qs finite box, ratio measured along the path: {W_meas:.4f} ± {dW:.4f}")
    wq = {}
    for u in ("u0.005", "u0.01", "u0.02"):
        w = [float(r) for r in pd.read_csv(os.path.join(ET, "level0_Wqs_20260911", u, "summary.csv"))["W_in_max"]]
        wq[u] = (np.mean(w), np.std(w, ddof=1) / math.sqrt(len(w)), len(w))
    for u, (m, s, n) in wq.items():
        g = m - W_meas; sg = math.hypot(s, dW)
        print(f"gap  <W>(u = {u[1:]}, {n} seeds) = {m:.4f} ± {s:.4f}  →  W − W_qs^finite = {g:+.4f} ± {sg:.4f} kT ({g/sg:+.1f} σ)")
