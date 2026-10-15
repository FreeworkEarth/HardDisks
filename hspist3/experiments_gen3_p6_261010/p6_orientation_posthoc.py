#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.43; plan-author decision 18, item 3) -- POST HOC, written after P6 ran (p6_run.py and
p6_tables.py were committed before the runs, 36f2b9b, and are not edited): the orientation of every P6 compartment through the hold,
from the runs' own psi6(t) columns (--gen3-psi6-compartments, every sigma-time). Information, no verdict.
Per run and compartment: the orientation at t = 0, 1000, 5000 and 2e4 sigma-time, and the first time it leaves 0 deg -- the earliest
hold row t0 with |orientation| > 15 deg in every row of [t0, t0 + 199] (200 sigma-time) -- or "never" (a turned start reads +-30 from
t = 0); then per cell and start the count of compartments that turned, and their Z in the last quarter (Q4, virial blocks 300-399)
split at P6's tables' cut between the two Z groups.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p6_261010/p6_orientation_posthoc.py [--data <P6 data root>]
"""
import argparse, csv, glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p6_run as R
import p6_tables as PT
TIMES = (0, 1000, 5000, 20000)


def turned_at(rows, c, thr=15.0, span=200):
    run = 0
    for x in rows:
        run = run + 1 if abs(float(x[f"psi6_orient_deg_{c}"])) > thr else 0
        if run >= span: return float(x["t_sigma"]) - (span - 1)
    return None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default=os.path.join(HERE, "data_f72d764")); a = ap.parse_args()
    print("# P6, POST HOC: the orientation of every compartment through the hold -- printed by experiments_gen3_p6_261010/p6_orientation_posthoc.py "
          "(written after the runs; information, no verdict)\n")
    print("leaves 0 deg = the earliest t0 with |orientation| > 15 deg in every psi6(t) row of [t0, t0 + 199] sigma-time; 'from the start' for a "
          "compartment already beyond 15 deg at t = 0; Z = Q4 (blocks 300-399)\n")
    for c in map(R.cell, R.cells()):
        for lab, o, l in R.starts_for(c):
            print(f"## {c['name']} {lab}{' (ADDITION)' if lab == 'zero' else ''}\n")
            print("| seed | compartment | orientation at t = 0 / 1e3 / 5e3 / 2e4 [deg] | leaves 0 deg at [sigma-time] | Z (Q4) |\n|---|---|---|---|---|")
            n_turn, n = 0, 0
            for r in range(6):
                d = os.path.join(os.path.abspath(a.data), c["name"], lab, f"seed{r}")
                f = glob.glob(os.path.join(d, "psi6_t_*_run0.csv"))
                if len(f) != 1: print(f"| {r} | - | no psi6(t) file | - | - |"); continue
                rows = [x for x in csv.DictReader(open(f[0])) if x["phase"] == "hold"]
                pick = {round(float(x["t_sigma"])): x for x in rows}; z = PT.zwin(d, PT.Q4)
                for k in (0, 1):
                    t0 = turned_at(rows, k); n += 1
                    o0 = abs(float(rows[0][f"psi6_orient_deg_{k}"]))
                    when = "from the start" if (t0 is not None and o0 > 15) else ("never" if t0 is None else f"{t0:.0f}")
                    if t0 is not None and o0 <= 15: n_turn += 1
                    ors = " / ".join(f"{float(pick[t][f'psi6_orient_deg_{k}']):+.1f}" if t in pick else "-" for t in TIMES)
                    print(f"| {r} | {k} | {ors} | {when} | {z[k]:.4f} |" if z else f"| {r} | {k} | {ors} | {when} | - |")
            print(f"\ncompartments that started near 0 deg and turned: {n_turn} of {n}\n")
    print("(POST HOC, information only; no verdict.)")


if __name__ == "__main__":
    main()
