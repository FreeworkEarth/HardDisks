#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.26; decision 12, part 3a): POST HOC, INFORMATION ONLY, no verdict -- written after
p3a_tables.py printed "NOT MET" for the held-cell criterion (2 of 6 runs with abs(z_LR) >= 2). The criterion and its result stand;
this script only describes the six held-cell runs further, from their existing virial-block files (nothing is run):
  1. the left - right difference of Z per block: lag-1 autocorrelations, and the z with blocks merged 2, 5 and 10 at a time (an
     autocorrelation makes 50-sigma-time blocks underestimate the SE), and the two halves of the hold;
  2. the pair-collision COUNTS per compartment and block (integers, assigned by the same per-disk compartment index as the virial,
     but with no virial arithmetic): left - right with its SE, whether left + right = the global count in every block, and the
     compartments' kinetic energies (each is conserved exactly by the held divider).
usage (from hspist3/ of the engine-gen3 worktree): python3 experiments_gen3_p3a_261009/p3a_heldcell_diagnostic.py --out <evidence dir>
"""
import argparse, csv, glob, math, os


def mse(v):
    m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)); return m, sd / math.sqrt(len(v))


def ac1(v):
    m = sum(v) / len(v); den = sum((x - m) ** 2 for x in v)
    return sum((v[i] - m) * (v[i + 1] - m) for i in range(len(v) - 1)) / den if den > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    print("# Part 3a, the held cells: post-hoc description (INFORMATION, no verdict), printed by experiments_gen3_p3a_261009/p3a_heldcell_diagnostic.py\n")
    print("## 1. Z_L - Z_R per block: autocorrelation, merged blocks, halves\n")
    print("| run | blocks | mean (SE) | z | lag-1 autocorr. Z_L / Z_R / difference | z, blocks merged x2 / x5 / x10 | first half (SE) | second half (SE) |\n|---|---|---|---|---|---|---|---|")
    runs = sorted(glob.glob(os.path.join(os.path.abspath(a.out), "e4", "*")))
    data = {}
    for d in runs:
        rows = [r for r in csv.DictReader(open(glob.glob(os.path.join(d, "virial_blocks_*.csv"))[0])) if r["phase"] == "hold"]
        data[d] = rows
        zl = [float(r["Z_0"]) for r in rows]; zr = [float(r["Z_1"]) for r in rows]; dd = [x - y for x, y in zip(zl, zr)]
        m, s = mse(dd); zs = []
        for k in (2, 5, 10):
            db = [sum(dd[i:i + k]) / k for i in range(0, len(dd) - k + 1, k)]; mb, sb = mse(db); zs.append(f"{mb / sb:+.2f}")
        h = len(dd) // 2; m1, s1 = mse(dd[:h]); m2, s2 = mse(dd[h:])
        print(f"| {os.path.basename(d)} | {len(dd)} | {m:+.5f} ({s:.5f}) | {m / s:+.2f} | {ac1(zl):+.2f} / {ac1(zr):+.2f} / {ac1(dd):+.2f} | {' / '.join(zs)} | "
              f"{m1:+.5f} ({s1:.5f}) | {m2:+.5f} ({s2:.5f}) |")
    print("\n## 2. Pair-collision counts per compartment and block (integers; no virial arithmetic)\n")
    print("| run | mean count L | mean count R | L - R (SE) | z | L + R = global in every block | KE_L values | KE_R values |\n|---|---|---|---|---|---|---|---|")
    for d in runs:
        rows = data[d]
        nl = [int(r["pair_events_0"]) for r in rows]; nr = [int(r["pair_events_1"]) for r in rows]; ng = [int(r["pair_events"]) for r in rows]
        dd = [x - y for x, y in zip(nl, nr)]; m, s = mse(dd)
        print(f"| {os.path.basename(d)} | {sum(nl) / len(nl):.1f} | {sum(nr) / len(nr):.1f} | {m:+.1f} ({s:.1f}) | {m / s:+.2f} | "
              f"{'yes' if all(x + y == g for x, y, g in zip(nl, nr, ng)) else '**no**'} | {', '.join(sorted({r['KE_0'] for r in rows}))} | "
              f"{', '.join(sorted({r['KE_1'] for r in rows}))} |")
    print("\n(POST HOC. The counts are assigned by the same per-disk compartment index as the virial; they agree with the Z differences "
          "in sign and size, so the bookkeeping of the virial is not their cause.)")


if __name__ == "__main__":
    main()
