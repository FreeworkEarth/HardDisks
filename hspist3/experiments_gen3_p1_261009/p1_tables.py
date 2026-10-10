#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.17; stage D): the tables of PILOT P1 (the structural clock), printed from its data. A pilot.
Per cell (N, eta), over its 6 seeds, from the held-divider psi6(t) series (psi6 global, every 1 sigma-time, 2e4 sigma-time):
  first quarter [0, T/4) and last quarter [3T/4, T]: the mean psi6 per seed; their means over seeds with the SE (SD/sqrt(n)); the drift
  = last - first, its SE from the per-seed differences (paired); "NOT STATIONARY within 2e4" if |drift| > 2 SE;
  the integrated autocorrelation time of psi6 in the second half [T/2, T] per seed (Sokal's automatic window: the smallest M with
  M >= 5 tau(M); tau = dt (1/2 + sum_{k=1..M} rho_k)), its mean over seeds and the seed-to-seed spread (SD, min, max); "TAU NOT
  RESOLVED" if the mean tau > (the second half's length)/20 = 500 sigma-time;
  the proposed T_eq = 10 x the mean tau (and 10 x the largest seed's tau);
  events per second through the driver (events / run_s of the run record), the cost table for the pre-registration;
  health: the runs with a run record clean=1 out of 6 (rule 10).
  INFORMATION column (added after the first cell, not part of the programme's criterion): the paired drift from the SECOND to the
  last quarter, Q4 - Q2, with its SE. The first quarter holds the melting or relaxing of the starting lattice; Q4 - Q2 separates
  that start transient from a slow drift. The programme's flag stays the Q4 - Q1 one.
usage (from hspist3/ of the engine-gen3 worktree): python3 experiments_gen3_p1_261009/p1_tables.py --out <data root> > p1_tables_output.txt
"""
import argparse, gzip, math, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p1_run as P1


def series(d, col=2):
    fs = [f for f in os.listdir(d) if f.startswith("psi6_t_") and (f.endswith(".csv.gz") or f.endswith(".csv"))]
    if not fs: return None
    f = os.path.join(d, fs[0]); op = gzip.open if f.endswith(".gz") else open
    t, v = [], []
    with op(f, "rt") as fh:
        next(fh)
        for l in fh:
            p = l.split(",")
            if p[1] != "hold": continue
            t.append(float(p[0])); v.append(float(p[col]))
    return np.array(t), np.array(v)


def tau_int(x, dt, c=5.0):
    x = np.asarray(x, float) - np.mean(x); n = len(x)
    if n < 10 or np.var(x) == 0: return math.nan
    f = np.fft.rfft(x, 2 * n); acf = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); rho = acf / acf[0]
    s = 0.5
    for M in range(1, n):
        s += rho[M]
        if M >= c * s: return dt * s
    return math.nan                                        # no window found: not resolved


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args(); O = os.path.abspath(a.out)
    T_HOLD = P1.HOLD_STEPS / 60.0
    print("# Pilot P1, the structural clock (261012 sec. 4.7.17), printed by experiments_gen3_p1_261009/p1_tables.py -- A PILOT\n")
    print(f"gen3, M4 seeding, divider held for {T_HOLD:g} sigma-time; psi6 every 1 sigma-time; 6 seeds per cell; SE = SD/sqrt(n) over seeds\n")
    for col, lab in ((2, "psi6 GLOBAL (|<psi6>| over all disks; the primary measure)"),
                     (3, "psi6 LOCAL (the mean of |psi6_j| per disk; orientation-independent, information)")):
        print(f"\n## {lab}\n")
        table(O, T_HOLD, col)


def table(O, T_HOLD, col):
    print("| N | eta | runs clean | psi6 first quarter | psi6 last quarter | drift (paired) | drift / SE | stationary within 2e4? | tau_int, "
          "second half [sigma] (mean; SD; min-max) | resolved? | proposed T_eq = 10 tau (mean; largest seed) [sigma] | events/s (driver) | "
          "information: drift Q4 - Q2 (paired) / SE |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in P1.cells():
        cd = os.path.join(O, c["name"])
        q1, q2, q4, taus, eps, clean = [], [], [], [], [], 0
        for r in range(P1.NSEED):
            d = os.path.join(cd, f"seed{r}")
            if not os.path.isdir(d): continue
            lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
            rec = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) (.*)$", lg, re.M)
            if rec and rec[-1][0] == "1":
                clean += 1
                kv = dict(re.findall(r"(\w+)=(\S+)", rec[-1][1]))
                ev = sum(int(kv[f]) for f in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band"))
                eps.append(ev / float(kv["run_s"]))
            S = series(d, col)
            if S is None: continue
            t, v = S
            q1.append(v[t < T_HOLD / 4].mean()); q4.append(v[t >= 3 * T_HOLD / 4].mean())
            q2.append(v[(t >= T_HOLD / 4) & (t < T_HOLD / 2)].mean())
            h = t >= T_HOLD / 2; dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
            taus.append(tau_int(v[h], dt))
        n = len(q1)
        if n < 2: print(f"| {c['N']} | {c['eta']:.3f} | {clean} of 6 | no data | | | | | | | | | |"); continue
        q1, q2, q4, taus = np.array(q1), np.array(q2), np.array(q4), np.array(taus)
        d42 = q4 - q2; dr42, se42 = d42.mean(), d42.std(ddof=1) / math.sqrt(n)
        se1, se4 = q1.std(ddof=1) / math.sqrt(n), q4.std(ddof=1) / math.sqrt(n)
        dd = q4 - q1; drift, sed = dd.mean(), dd.std(ddof=1) / math.sqrt(n)
        stat = "yes" if abs(drift) <= 2 * sed else "**NOT STATIONARY within 2e4**"
        tf = taus[np.isfinite(taus)]
        tm = tf.mean() if len(tf) else math.nan
        res = "yes" if len(tf) == n and tm <= (T_HOLD / 2) / 20 else "**TAU NOT RESOLVED**"
        tau_s = f"{tm:.3g}; {tf.std(ddof=1) if len(tf) > 1 else math.nan:.2g}; {tf.min() if len(tf) else math.nan:.3g}-{tf.max() if len(tf) else math.nan:.3g}" + \
                (f" ({n - len(tf)} seeds without a window)" if len(tf) < n else "")
        teq = f"{10 * tm:.3g}; {10 * tf.max() if len(tf) else math.nan:.3g}"
        print(f"| {c['N']} | {c['eta']:.3f} | {clean} of 6 | {q1.mean():.4f} +- {se1:.4f} | {q4.mean():.4f} +- {se4:.4f} | {drift:+.4f} +- {sed:.4f} | "
              f"{drift / sed if sed > 0 else math.nan:+.2f} | {stat} | {tau_s} | {res} | {teq} | {np.mean(eps) if eps else math.nan:.3g} | "
              f"{dr42:+.4f} +- {se42:.4f} / {dr42 / se42 if se42 > 0 else math.nan:+.2f} |")


if __name__ == "__main__":
    main()
