#!/usr/bin/env python3
"""##CHRIS 2026-10-04 (Task Y, 261012 sec. 3): checks for the "A-fixed" design, on the Mac.

PRE-CHECK (existing data, method-A pilot, 20 runs): during the 200 sigma-time hold of every method-A run the divider is
immovable -- every divider event D0 with t < 200 has u_wall = 0 and dE = 0 -- which is what A-fixed extends to the whole
record (code: 00ALLINONE.c:16890-16892 mass 0 and velocity 0 while held; edmd.c:1174-1183 the infinite-mass branch).

GATE G1 (after the A-fixed pilot cell is fetched with `bash hspist3/cluster/confinement_20261013/fetch_afix.sh pilot`):
for each of its 20 runs (pi/8 anchor, seeds 9700-9703 at the five positions; the method-A pilot's tasks, held)
  (a) hold-phase identity: the event-log lines with t < 200 are identical to the method-A pilot's (same seed, same
      position; the two command lines differ only in --wall-hold-steps and --steps);
  (b) immobility: every D0 event with t < t1 = 5200 has u_wall == 0 (event log) and reduce_AF.py's u_wall_max == 0;
  (c) no work: the sum of dE over those events is 0 exactly (W_div == 0);
  (d) the record is complete: t_last >= t1, window = 5000;
  (e) health: no EDMD-HEALTH line in run_<seed>.log;
  (f) the drift check of sec. 3: |1 - f| < 0.002, f from the recorded temperatures by the C4 formula
      (paper1_confinement_heldwall_posthoc_261004.drift); the max |T - 1| is printed for information.
usage (from hspist3/): python3 cluster/afix_pilot_check_261004.py
"""
import glob, math, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
A_P = os.path.join(HS, "experiments_energy_transfer/paper1_confinement_A_20261013/pilot_epi8_H_H10_L10")
AF_P = os.path.join(HS, "experiments_energy_transfer/paper1_confinement_Afix_261004/pilot_epi8_H_H10_L10")
POS = (("m2", -2), ("m1", -1), ("0", 0), ("p1", 1), ("p2", 2))
T0, T1 = 200.0, 312000 * 0.4 / 24.0


def hold_lines(ev):
    out = []
    with open(ev) as fh:
        next(fh)
        for l in fh:
            if float(l.split(",", 1)[0]) >= T0: break
            out.append(l)
    return out


def main():
    print("### Pre-check on existing data: the hold phase (t < 200) of the 20 method-A pilot runs\n")
    print("| position | seed | D0 events, t < 200 | max abs u_wall | sum dE |\n|---|---|---|---|---|")
    ok_pre = True
    for lab, j in POS:
        for f in sorted(glob.glob(os.path.join(A_P, f"x_{lab}", "ev_*.csv"))):
            e = pd.read_csv(f, usecols=["t_sigma", "kind", "u_wall", "dE"]); d = e[(e.kind == "D0") & (e.t_sigma < T0)]
            um, w = float(d.u_wall.abs().max()), float(d.dE.sum()); ok_pre &= (um == 0.0 and w == 0.0 and len(d) > 0)
            print(f"| x_{lab} | {f[-8:-4]} | {len(d)} | {um:g} | {w:g} |")
    print(f"\npre-check: {'the divider is immovable during the hold in all 20 runs' if ok_pre else '**FAIL**'}")
    if not os.path.isdir(AF_P):
        print(f"\nGATE G1: not run -- the A-fixed pilot cell is not on the Mac yet ({os.path.relpath(AF_P, HS)})")
        return 0 if ok_pre else 1
    print("\n### GATE G1 -- the A-fixed pilot (sec. 3)\n")
    print("| position | seed | (a) hold lines identical to method A | (b) max abs u_wall, t < t1 | (c) sum dE | (d) t_last | (e) health | max abs(T - 1) |")
    print("|---|---|---|---|---|---|---|---|")
    ok = True; rows = []
    for lab, j in POS:
        for f in sorted(glob.glob(os.path.join(AF_P, f"x_{lab}", "ev_*.csv"))):
            s = f[-8:-4]; fa = os.path.join(A_P, f"x_{lab}", f"ev_{s}.csv")
            same = hold_lines(f) == hold_lines(fa)
            e = pd.read_csv(f, usecols=["t_sigma", "kind", "u_wall", "dE"]); d = e[(e.kind == "D0") & (e.t_sigma < T1)]
            um, w, tl = float(d.u_wall.abs().max()), float(d.dE.sum()), float(e.t_sigma.max())
            r = pd.read_csv(os.path.join(AF_P, f"x_{lab}", f"red_{s}.csv")).iloc[0]
            h = open(os.path.join(AF_P, f"x_{lab}", f"run_{s}.log"), errors="ignore").read().count("EDMD-HEALTH")
            dT = max(abs(r.T_L - 1), abs(r.T_R - 1)); rows.append((j, r))
            good = same and um == 0.0 and w == 0.0 and r.u_wall_max == 0.0 and r.W_div == 0.0 and tl >= T1 and h == 0
            ok &= good
            print(f"| x_{lab} | {s} | {'yes' if same else '**NO**'} | {um:g} | {w:g} | {tl:.2f} | {h} | {dT:.2e} |")
    # (f) drift check by the C4 formula: delta_j = mean of -N_s dT_L/F_L and +N_s dT_R/F_R; 1 - f = -sum delta x / sum x^2
    Ns, dL = 50, 0.125; num = den = 0.0
    for j in (-2, -1, 1, 2):
        R = pd.DataFrame([r for jj, r in rows if jj == j])
        dLft = -Ns * (R.T_L.mean() - 1) / R.F_L.mean(); dRgt = Ns * (R.T_R.mean() - 1) / R.F_R.mean()
        num += 0.5 * (dLft + dRgt) * j * dL; den += (j * dL) ** 2
    one_f = -num / den
    print(f"\n(f) drift check: |1 - f| = {abs(one_f):.2e} (must be < 0.002)")
    ok &= abs(one_f) < 0.002
    print(f"\n**GATE G1: {'PASS' if ok and ok_pre else 'FAIL'}**")
    return 0 if ok and ok_pre else 1


if __name__ == "__main__":
    sys.exit(main())
