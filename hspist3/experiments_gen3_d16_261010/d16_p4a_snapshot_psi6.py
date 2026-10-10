#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.37; plan-author decision 16, item 2b -- AN ADDITION, not ordered): psi6 PER COMPARTMENT from
P4a's existing snapshots -- information, no verdict, no runs. P4 wrote snapshots every 2000 sigma-time for seeds 0 and 1 of each start
(not committed; SHA-256 in P4_SHA256_FILES.txt): start A seed 1 is one of the three that end at global psi6 ~ 0.06, A seed 0 one of
the three that end at ~ 0.72; B seeds 0 and 1. From the positions, with the code's definition (00ALLINONE.c compute_psi6_from_edmd:
psi6(j) = the mean of exp(6 i theta_jk) over the neighbours k within 1.4 sigma, all disks; the global value |sum_j psi6(j)| / N):
  per compartment c (disks left / right of the divider): Psi_c = sum_{j in c} psi6(j) / N_c, its modulus, and the angle between the
  two compartments' orientations, arg(Psi_L / Psi_R) / 6 in degrees (in (-30, 30]: 30 means orientations 30 degrees apart, i.e. the
  same lattice turned by 90 degrees, so that Psi_L ~ -Psi_R and the global value cancels);
  each compartment's orientation against the x axis, arg(Psi_c) / 6 in degrees (in (-30, 30]; 0 for M4's 'parallel' rows);
  the global |psi6| from the same positions, next to the psi6(t) file's row at that time (a check of the reading);
and for all twelve runs of the cell (no snapshots needed) Z per compartment in the last tenth [9e4, 1e5) from the existing virial
blocks (blocks 1800-1999 of 50 sigma-time, as P4's tables; the plain mean of the blocks' Z_0 and Z_1, as part 3a's table), with
the global Z.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d16_261010/d16_p4a_snapshot_psi6.py [--data <P4 data root>]
"""
import argparse, cmath, csv, glob, math, os, struct, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
P4DATA = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
CELL, CUT, TIMES = "N400_eta0.740", 1.4, (0.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0, 60000.0, 70000.0, 80000.0, 90000.0, 100000.0)
RUNS = (("A", 1, "ends at ~ 0.06"), ("A", 0, "ends at ~ 0.72"), ("B", 0, "compressed start"), ("B", 1, "compressed start"))


def read_snaps(p):
    """the G3SNAP1 file (as experiments_gen3_p3a_261009/p3a_tables.py reads it): header N, dividers, box W, H, radius [sigma], interval;
    then per snapshot t [sigma-time], phase, the divider positions [sigma], the positions x, y [sigma]"""
    b = open(p, "rb").read()
    assert b[:8] == b"G3SNAP1\n", "magic"
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, T = struct.unpack_from("<4d", b, 16); o = 48; out = []
    while o < len(b):
        t, = struct.unpack_from("<d", b, o); ph, k = struct.unpack_from("<ii", b, o + 8); o += 16
        dx = struct.unpack_from(f"<{k}d", b, o); o += 8 * k
        xy = struct.unpack_from(f"<{2 * N}d", b, o); o += 16 * N
        out.append((t, ph, dx, xy))
    return (N, nd, W, H, R, T), out


def psi6(xy, N, R, xdiv):
    cut2 = (CUT * 2.0 * R) ** 2; ps = []
    for i in range(N):
        xi, yi = xy[2 * i], xy[2 * i + 1]; s, n = 0j, 0
        for j in range(N):
            if j == i: continue
            dx, dy = xy[2 * j] - xi, xy[2 * j + 1] - yi; d2 = dx * dx + dy * dy
            if 0.0 < d2 <= cut2: s += cmath.exp(6j * math.atan2(dy, dx)); n += 1
        ps.append(s / n if n else 0j)
    L = [p for i, p in enumerate(ps) if xy[2 * i] < xdiv]; Rr = [p for i, p in enumerate(ps) if xy[2 * i] >= xdiv]
    PL, PR = sum(L) / len(L), sum(Rr) / len(Rr)
    ang = math.degrees(cmath.phase(PL / PR)) / 6.0 if abs(PL) > 0 and abs(PR) > 0 else float("nan")
    return dict(nl=len(L), nr=len(Rr), PL=abs(PL), PR=abs(PR), ang=ang, glob=abs(sum(ps)) / N,
                oL=math.degrees(cmath.phase(PL)) / 6.0, oR=math.degrees(cmath.phase(PR)) / 6.0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default=P4DATA); a = ap.parse_args()
    print("# Item 2b, an addition (261012 sec. 4.7.37; decision 16; not ordered): psi6 per compartment from P4a's existing snapshots -- "
          "printed by experiments_gen3_d16_261010/d16_p4a_snapshot_psi6.py; information, no verdict, no runs\n")
    print(f"cell {CELL}; neighbours within {CUT} sigma (all disks, as the code); the snapshots every 2000 sigma-time, shown every 1e4\n")
    print("| start | seed | group | t [sigma-time] | N_L / N_R | abs(Psi_L) | abs(Psi_R) | orientation L [deg] | orientation R [deg] | angle L vs R [deg] | global abs(psi6) (snapshot) | global (psi6(t) file) |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
    for st, r, grp in RUNS:
        d = os.path.join(os.path.abspath(a.data), CELL, st, f"seed{r}")
        sp = glob.glob(os.path.join(d, "snapshots_*_run0.bin")); pf = glob.glob(os.path.join(d, "psi6_t_*_run0.csv"))
        if len(sp) != 1: print(f"| {st} | {r} | {grp} | no snapshot file |"); continue
        (N, nd, W, H, R, T), snaps = read_snaps(sp[0])
        rows = {round(float(x["t_sigma"])): float(x["psi6_global"]) for x in csv.DictReader(open(pf[0])) if x["phase"] == "hold"} if len(pf) == 1 else {}
        for t, ph, dx, xy in snaps:
            if not any(abs(t - tt) < 1e-6 * max(1.0, tt) + 1e-6 for tt in TIMES) and t != snaps[-1][0]: continue
            v = psi6(xy, N, R, dx[0])
            f = rows.get(round(t)); fs = f"{f:.4f}" if f is not None else "-"
            print(f"| {st} | {r} | {grp} | {t:.0f} | {v['nl']} / {v['nr']} | {v['PL']:.4f} | {v['PR']:.4f} | {v['oL']:+.1f} | {v['oR']:+.1f} | {v['ang']:+.1f} | {v['glob']:.4f} | {fs} |")
    print("\n### Z per compartment in the last tenth [9e4, 1e5) sigma-time (the existing virial blocks; plain block means)\n")
    print("| start | seed | blocks | Z_L | Z_R | Z_L - Z_R | Z_global |\n|---|---|---|---|---|---|---|")
    for st in ("A", "B"):
        for r in range(6):
            vb = glob.glob(os.path.join(os.path.abspath(a.data), CELL, st, f"seed{r}", "virial_blocks_*_run0.csv"))
            if len(vb) != 1: print(f"| {st} | {r} | no file |"); continue
            bl = [x for x in csv.DictReader(open(vb[0])) if x["phase"] == "hold" and 1800 <= int(x["block"]) < 2000]
            zl, zr, zg = (sum(float(x[k]) for x in bl) / len(bl) for k in ("Z_0", "Z_1", "Z_global"))
            print(f"| {st} | {r} | {len(bl)} | {zl:.4f} | {zr:.4f} | {zl - zr:+.4f} | {zg:.4f} |")
    print("\n(Information only; no verdict.)")


if __name__ == "__main__":
    main()
