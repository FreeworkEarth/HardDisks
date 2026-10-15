#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.37; plan-author decision 16, item 2b): P4a's LOCAL psi6 -- information, no verdict, no runs.
From pilot P4's part a (N = 400, eta 0.740, the divider held 1e5 sigma-time; experiments_gen3_p4_261010/data_06e7248/N400_eta0.740),
the existing psi6(t) files (psi6_t_*_run0.csv: t_sigma, phase, psi6_global, psi6_local, mean_neighbors; written every sigma-time):
per run the means over the hold rows of the last tenth, 9e4 <= t < 1e5 sigma-time, of
  psi6_global    |<psi6(j)>| over all disks (long-range orientational order),
  psi6_local     <|psi6(j)|> over the disks with at least one neighbour (local order),
  mean_neighbors the mean number of neighbours within 1.4 sigma (the code's cutoff, PSI6_CUTOFF_SIGMA),
for the three start-A seeds that end at psi6 ~ 0.06 (seeds 1, 3, 4), the three that end at ~ 0.72 (seeds 0, 2, 5) and start B's six
seeds; then each group's mean (SE over its seeds). The grouping is the plan author's; the global psi6 column shows it.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d16_261010/d16_p4a_local_psi6.py [--data <P4 data root>]
"""
import argparse, csv, glob, math, os

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
P4DATA = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
CELL, T0, T1 = "N400_eta0.740", 9.0e4, 1.0e5
GROUPS = (("A, ending at psi6 ~ 0.06", "A", (1, 3, 4)), ("A, ending at psi6 ~ 0.72", "A", (0, 2, 5)), ("B (compressed start)", "B", (0, 1, 2, 3, 4, 5)))
COLS = ("psi6_global", "psi6_local", "mean_neighbors")


def mse(v):
    m = sum(v) / len(v); return m, (math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / math.sqrt(len(v)) if len(v) > 1 else float("nan"))


def last_tenth(d):
    f = glob.glob(os.path.join(d, "psi6_t_*_run0.csv"))
    if len(f) != 1: return None
    rows = [r for r in csv.DictReader(open(f[0])) if r["phase"] == "hold" and T0 <= float(r["t_sigma"]) < T1]
    return dict(n=len(rows), **{c: sum(float(r[c]) for r in rows) / len(rows) for c in COLS}) if rows else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default=P4DATA); a = ap.parse_args()
    O = os.path.join(os.path.abspath(a.data), CELL)
    print("# Item 2b (261012 sec. 4.7.37; decision 16): P4a's local psi6 and mean neighbour number in the last tenth of the hold -- printed "
          "by experiments_gen3_d16_261010/d16_p4a_local_psi6.py from the existing psi6(t) files; information, no verdict, no runs\n")
    print(f"cell {CELL} (N = 400, H = 20, eta 0.740, the divider held 1e5 sigma-time); the hold rows with {T0:.0f} <= t < {T1:.0f} sigma-time; "
          "neighbours within 1.4 sigma\n")
    print("| group | start | seed | rows | psi6 global | psi6 local | mean neighbours |\n|---|---|---|---|---|---|---|")
    G = {}
    for lab, st, seeds in GROUPS:
        for r in seeds:
            v = last_tenth(os.path.join(O, st, f"seed{r}"))
            if v is None: print(f"| {lab} | {st} | {r} | **missing** | - | - | - |"); continue
            G.setdefault(lab, []).append(v)
            print(f"| {lab} | {st} | {r} | {v['n']} | {v['psi6_global']:.4f} | {v['psi6_local']:.4f} | {v['mean_neighbors']:.4f} |")
    print("\n### Group means (SE over the group's seeds)\n")
    print("| group | seeds | psi6 global | psi6 local | mean neighbours |\n|---|---|---|---|---|")
    for lab, st, seeds in GROUPS:
        vs = G.get(lab, [])
        if not vs: continue
        print(f"| {lab} | {len(vs)} | " + " | ".join("{:.4f} ({:.4f})".format(*mse([v[c] for v in vs])) for c in COLS) + " |")
    print("\n(Information only; no verdict.)")


if __name__ == "__main__":
    main()
