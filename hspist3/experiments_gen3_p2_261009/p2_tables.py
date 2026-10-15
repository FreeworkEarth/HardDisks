#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.19; stage F): the table of PILOT P2 next to A1 v2's rows of sec. 4.7.11. A pilot.
Per eta_true (all masses together, as sec. 4.7.11's table 2) and per mass: trajectories, mean dnu = nu2/nu1 - 1 with its SE, the
number with dnu exactly 0 (one argmax bin), mean dA = A2/A1 - 1 with its SE, psi6 global at the release (hold) and at the end (means
over the trajectories with their SE); the estimator is sec. 4.7.11's own (validation/paper1_window_aging_261009.traj2: the canonical
argmax over TD periods and the same estimator on each half of that record; A = RMS displacement about the half's mean). A1 v2's rows of
sec. 4.7.11 (printed there by paper1_window_aging_261009.py) are quoted next to them.
usage (from hspist3/ of the engine-gen3 worktree): python3 experiments_gen3_p2_261009/p2_tables.py --out <root> > p2_tables_output.txt
"""
import argparse, csv, glob, math, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p2_run as P2
sys.path.insert(0, os.path.join(P2.MAIN_HS, "validation"))
import paper1_window_aging_261009 as WA
A1V2_ROWS = {0.6905: "-0.41 (0.22) | 33 (15 %) | +2.85 (1.12)", 0.7060: "+7.45 (0.78) | 15 (7 %) | -9.63 (1.21)",
             0.7113: "+10.37 (0.85) | 18 (8 %) | -13.87 (1.49)", 0.7167: "+7.78 (0.85) | 23 (10 %) | -2.31 (1.76)"}   # sec. 4.7.11, table 2


def mse(x):
    x = np.asarray([v for v in x if v is not None and math.isfinite(v)], float)
    return (x.mean(), x.std(ddof=1) / math.sqrt(len(x)), len(x)) if len(x) > 1 else (math.nan, math.nan, len(x))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args(); O = os.path.abspath(a.out)
    print("# Pilot P2, drift with an equilibrated start (261012 sec. 4.7.19), printed by experiments_gen3_p2_261009/p2_tables.py -- A PILOT\n")
    print("| eta_true | M | trajectories (clean / runs) | T_eq [sigma] | mean dnu [%] (SE) | dnu exactly 0 | mean dA [%] (SE) | psi6 at release | "
          "psi6 at the end | A1 v2 (sec. 4.7.11, all masses): mean dnu [%] (SE) / dnu = 0 / mean dA [%] (SE) |\n|---|---|---|---|---|---|---|---|---|---|")
    for eta_t, L0a, eta_p1 in P2.A1V2:
        allr = []
        for M in P2.MASSES + ("all",):
            rows = []
            if M == "all": rows = allr
            else:
                cdir = os.path.join(O, f"eta{eta_t:.4f}_M{M}")
                for r in range(P2.NSEED):
                    d = os.path.join(cdir, f"seed{r}")
                    if not os.path.isdir(d): continue
                    lg = open(os.path.join(d, "run.log"), errors="ignore").read()
                    clean = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M) == ["1"]
                    tr = glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv"))
                    ps = os.path.join(d, "speed_of_sound_psi6.csv")
                    res = WA.traj2(tr[0]) if tr else None
                    p_h = p_e = math.nan
                    if os.path.exists(ps):
                        row = list(csv.DictReader(open(ps)))[-1]; p_h, p_e = float(row["psi6_global_hold"]), float(row["psi6_global_end"])
                    teq = None
                    m = re.search(r"T_eq=(\S+) sigma-time", lg)
                    if m: teq = float(m.group(1))
                    rows.append(dict(clean=clean, res=res, ph=p_h, pe=p_e, teq=teq))
                allr += rows
            ok = [x for x in rows if x["clean"] and isinstance(x["res"], dict)]      # traj2: dict(nu, nu1, nu2, TD, a1, a2), None, "refused"
            dn = [100 * (x["res"]["nu2"] / x["res"]["nu1"] - 1) for x in ok if x["res"]["nu1"]]
            da = [100 * (x["res"]["a2"] / x["res"]["a1"] - 1) for x in ok if x["res"]["a1"]]
            n0 = sum(1 for x in ok if x["res"]["nu1"] == x["res"]["nu2"])
            m1, s1, _ = mse(dn); m2, s2, _ = mse(da); h1, hs, _ = mse([x["ph"] for x in ok]); e1, es, _ = mse([x["pe"] for x in ok])
            teqs = sorted({x["teq"] for x in rows if x["teq"]})
            a1 = A1V2_ROWS[eta_t] if M == "all" else ""
            print(f"| {eta_t:.4f} | {M} | {len(ok)} / {len(rows)} | {', '.join(f'{t:g}' for t in teqs)} | {m1:+.2f} ({s1:.2f}) | {n0} | {m2:+.2f} ({s2:.2f}) | "
                  f"{h1:.3f} ({hs:.3f}) | {e1:.3f} ({es:.3f}) | {a1} |")
    print("\n(PILOT. dnu = nu2/nu1 - 1 and dA = A2/A1 - 1 between the halves of each record, sec. 4.7.11's estimator; a settled start is "
          "lattice + jitter + a held hold of T_eq; A1 v2 started from the driver's grid seeding with a 2000-step hold.)")


if __name__ == "__main__":
    main()
