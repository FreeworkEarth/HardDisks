#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.38; plan-author decision 17, "additions to P5's report") -- POST HOC: written while P5 ran,
after its registration (p5_run.py, p5_tables.py, engine-gen3 bbaf904), and labelled so; existing data only; a pilot, no paper use.
Z per compartment, since at N = 400 the Z of a compartment tells its orientation (sec. 4.7.37: 14.33 for M4's 0-degree lattice, 13.29
for the turned one). Per seed, in the last quarter of the hold [1.5e4, 2e4) sigma-time (virial blocks 300-399 of 50 sigma-time):
Z_L and Z_R (compartments 0 and 1), Z_L - Z_R and Z_global, each the time average over the quarter's blocks (dt-weighted, as
p4_tables.run() averages Z_global), for every P5 cell (start B at u = 0.001 or 0.0001) and for P4's N = 1600 cells (starts A and B at
u = 0.01, the same quarter of the same 2e4 hold). Then per cell and start: the mean (SE over seeds) of Z_L, Z_R, Z_global, the mean of
abs(Z_L - Z_R), and the smallest and largest compartment value.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p5_261010/p5_compartments_posthoc.py [--data <P5 data root>] [--p4 <P4 data root>]
"""
import argparse, csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import p5_run as R
P5DATA = os.path.join(HERE, "data_06e7248"); P4DATA = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
B0, B1 = 300, 400                                                      # the last quarter's blocks


def quarter(d):
    f = glob.glob(os.path.join(d, "virial_blocks_*_run0.csv"))
    if len(f) != 1: return None
    bl = [x for x in csv.DictReader(open(f[0])) if x["phase"] == "hold" and B0 <= int(x["block"]) < B1]
    if len(bl) != B1 - B0: return None
    w = [float(x["dt_sigma"]) for x in bl]; W = sum(w)
    zl, zr, zg = (sum(float(x[k]) * wi for x, wi in zip(bl, w)) / W for k in ("Z_0", "Z_1", "Z_global"))
    return dict(zl=zl, zr=zr, zg=zg)


def mse(v):
    m = sum(v) / len(v); return m, (math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / math.sqrt(len(v)) if len(v) > 1 else float("nan"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default=P5DATA); ap.add_argument("--p4", default=P4DATA); a = ap.parse_args()
    print("# P5, post hoc (decision 17): Z per compartment in the last quarter of the hold -- printed by "
          "experiments_gen3_p5_261010/p5_compartments_posthoc.py; POST HOC (written while P5 ran), existing data only; a pilot, no paper use\n")
    print("Last quarter [1.5e4, 2e4) sigma-time = virial blocks 300-399; dt-weighted time averages; L = compartment 0 (left of the divider)\n")
    rows = []
    for c in map(R.cell, R.cells()):
        rows.append((c["name"], f"B, u = {c['u']:g} (P5)", [quarter(os.path.join(os.path.abspath(a.data), c["name"], "B", f"seed{r}")) for r in range(6)]))
    for p4name in ("N1600_eta0.740", "N1600_eta0.780"):
        for st, lab in (("A", "A, the M4 lattice (P4)"), ("B", "B, u = 0.01 (P4)")):
            rows.append((p4name, lab, [quarter(os.path.join(os.path.abspath(a.p4), p4name, st, f"seed{r}")) for r in range(6)]))
    print("## Per seed\n\n| cell | start | seed | Z_L | Z_R | Z_L - Z_R | Z_global |\n|---|---|---|---|---|---|---|")
    for cell, lab, vs in rows:
        for r, v in enumerate(vs):
            print(f"| {cell} | {lab} | {r} | " + (f"{v['zl']:.4f} | {v['zr']:.4f} | {v['zl'] - v['zr']:+.4f} | {v['zg']:.4f} |" if v else "- | - | - | - |"))
    print("\n## Per cell and start (mean (SE) over seeds)\n\n| cell | start | n | Z_L | Z_R | Z_global | mean abs(Z_L - Z_R) | smallest / largest "
          "compartment value |\n|---|---|---|---|---|---|---|---|")
    for cell, lab, vs in rows:
        ok = [v for v in vs if v]
        if not ok: print(f"| {cell} | {lab} | 0 | - | - | - | - | - |"); continue
        comp = [x for v in ok for x in (v["zl"], v["zr"])]
        print(f"| {cell} | {lab} | {len(ok)} | " + " | ".join("{:.4f} ({:.4f})".format(*mse([v[k] for v in ok])) for k in ("zl", "zr", "zg"))
              + f" | {sum(abs(v['zl'] - v['zr']) for v in ok) / len(ok):.4f} | {min(comp):.4f} / {max(comp):.4f} |")
    print("\n## N = 1600: the compartment values split at their largest gap (all starts of the cell pooled), and the count per start\n")
    print("| eta | the lower group | the gap | the higher group | start | compartments in the higher group |\n|---|---|---|---|---|---|")
    for eta in ("0.740", "0.780"):
        grp = [(lab, [x for v in vs if v for x in (v["zl"], v["zr"])]) for cell, lab, vs in rows if cell.startswith(f"N1600_eta{eta}")]
        allv = sorted(x for _, xs in grp for x in xs)
        g, k = max((allv[i + 1] - allv[i], i) for i in range(len(allv) - 1)); cut = 0.5 * (allv[k] + allv[k + 1])
        for j, (lab, xs) in enumerate(grp):
            head = (f"| {eta} | {allv[0]:.4f} - {allv[k]:.4f} | {g:.4f} | {allv[k + 1]:.4f} - {allv[-1]:.4f} |" if j == 0 else "| | | | |")
            print(f"{head} {lab} | {sum(x > cut for x in xs)} of {len(xs)} |")
    print("\n(POST HOC, information only; no verdict.)")


if __name__ == "__main__":
    main()
