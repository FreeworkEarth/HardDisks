#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.43; plan-author decision 18, item 3): the tables of PILOT P6 -- a pilot, no paper use, no
verdict; the two questions are answered with the numbers only. Per compartment (c = 0 left of the divider, 1 right), in the last
quarter of the 2e4-sigma-time hold (Q4 = [1.5e4, 2e4) sigma-time; virial blocks 300-399 of 50 sigma-time; for P4a, which held 1e5,
the same window, and also its first quarter [0, 5000), before any compartment turned):
  Z_c        the dt-weighted time average of the block values Z_c (as p5_compartments_posthoc.py)
  |Psi_c|    P6 (--gen3-psi6-compartments): the mean of the psi6_c column over the Q4 rows; P4 and P5 (no such columns): from their
             snapshots (seeds 0 and 1 only, every 2000 sigma-time) read by --gen3-psi6-snapshots of the same binary, the snapshots at
             1.5e4 <= t <= 2e4 (for P4a's first quarter, 0 <= t <= 5000); other seeds: Z only
  orientation_c  arg of the mean of |Psi_c| exp(6 i orient_c) over the same rows or snapshots, / 6, in degrees in (-30, 30] (0: the
             close-packed lines along x; +-30: along y, turned)
Rows per cell: P6's turned start (ordered); P6's 0-deg start (ADDITION, N = 1600 only); P4's lattice start A (the present M4 family:
0 deg at N = 400, turned at N = 1600 -- sec. 4.7.42); the compressed starts P4 B (u = 0.01) and P5 B (u = 0.001, 0.0001).
Tables: the inventory; per cell and start, per seed: Z_0, Z_1 and, where known, |Psi| and the orientation; per cell and start the mean
and SE over its compartments (each compartment one sample; the divider is held, so the two compartments of a run are independent);
for the compressed starts at N = 1600 also their two Z groups (split at the largest gap, as p5_compartments_posthoc.py).
THE TWO QUESTIONS, numbers only:
  Q1 is the compressed start's state the turned crystal? The turned start's mean Z_c against the compressed starts' (all compartments,
     and at N = 1600 each Z group), difference (SE) and its ratio to the SE; and the compressed compartments' orientations from
     their snapshots.
  Q2 which orientation has the lower Z? The 0-deg start's mean Z_c against the turned start's, difference (SE), ratio, and the lower
     one named. N = 400: 0 deg = P4a's start A (Q1 window, all compartments at 0; and Q4) against P6's turned start (Q4). N = 1600:
     0 deg = P6's 0-deg start (the addition) against P6's turned start, both Q4.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p6_261010/p6_tables.py --data <P6 data root> --bin <00ALLINONE with --gen3-psi6-snapshots>
  python3 experiments_gen3_p6_261010/p6_tables.py --selftest <empty folder>
"""
import argparse, cmath, csv, glob, io, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import p6_run as R
P4D = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
P5D = os.path.join(HS, "experiments_gen3_p5_261010", "data_06e7248")
Q4, Q1 = (300, 400, 15000.0, 20000.0), (0, 100, 0.0, 5000.0)    # (first block, end block, t0, t1)


def one(d, pat):
    g = sorted(glob.glob(os.path.join(d, pat))); return g[0] if len(g) == 1 else None


def zwin(d, w):
    """the dt-weighted Z_0, Z_1 over the blocks [w0, w1)"""
    f = one(d, "virial_blocks_*_run0.csv")
    if not f: return None
    bl = [x for x in csv.DictReader(open(f)) if x["phase"] == "hold" and w[0] <= int(x["block"]) < w[1]]
    if len(bl) != w[1] - w[0]: return None
    W = sum(float(x["dt_sigma"]) for x in bl)
    return tuple(sum(float(x[k]) * float(x["dt_sigma"]) for x in bl) / W for k in ("Z_0", "Z_1"))


def cmean(rows):
    """[(|Psi|, orientation deg)] -> (mean |Psi|, orientation of the mean complex Psi)"""
    if not rows: return None
    z = sum(a * cmath.exp(6j * math.radians(o)) for a, o in rows) / len(rows)
    return sum(a for a, _ in rows) / len(rows), math.degrees(cmath.phase(z)) / 6.0


def psi_cols(d, w):
    """P6: per compartment (mean |Psi_c|, orientation) over the psi6(t) rows of the hold with t0 <= t < t1"""
    f = one(d, "psi6_t_*_run0.csv")
    if not f: return None
    rows = [x for x in csv.DictReader(open(f)) if x["phase"] == "hold" and w[2] <= float(x["t_sigma"]) < w[3]]
    if not rows or "psi6_0" not in rows[0]: return None
    return tuple(cmean([(float(x[f"psi6_{c}"]), float(x[f"psi6_orient_deg_{c}"])) for x in rows]) for c in (0, 1))


def psi_snaps(d, w, binp):
    """P4, P5: per compartment (mean |Psi_c|, orientation) over the snapshots with t0 <= t <= t1 (hold), by --gen3-psi6-snapshots"""
    f = one(d, "snapshots_*_run0.bin")
    if not f: return None
    out = subprocess.run([binp, f"--gen3-psi6-snapshots={f}"], capture_output=True, text=True).stdout
    rows = [x for x in csv.DictReader(io.StringIO("\n".join(l for l in out.splitlines() if not l.startswith("#"))))
            if x["phase"] == "hold" and w[2] <= float(x["t_sigma"]) <= w[3] + 0.5]
    if not rows: return None
    return tuple(cmean([(float(x[f"psi6_{c}"]), float(x[f"psi6_orient_deg_{c}"])) for x in rows]) for c in (0, 1))


def mse(v):
    v = [x for x in v if x is not None]
    if not v: return float("nan"), float("nan"), 0
    m = sum(v) / len(v)
    return m, (math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / math.sqrt(len(v)) if len(v) > 1 else float("nan")), len(v)


def groups(vals):
    """split at the largest gap of the sorted values: (cut, gap)"""
    s = sorted(vals)
    if len(s) < 2: return None, None
    g, k = max((s[i + 1] - s[i], i) for i in range(len(s) - 1)); return 0.5 * (s[k] + s[k + 1]), g


def starts(c, P6):
    """(label, kind, run folders, window, psi source) per start of a cell"""
    n = c["name"]; out = [("P6 turned (ordered)", "turned", [os.path.join(P6, n, "turned", f"seed{r}") for r in range(6)], Q4, "cols")]
    if c["N"] == 1600: out.append(("P6 0-deg (ADDITION)", "zero", [os.path.join(P6, n, "zero", f"seed{r}") for r in range(6)], Q4, "cols"))
    out.append((f"P4 A, the present family ({'0 deg' if c['N'] == 400 else 'turned'})", "A", [os.path.join(P4D, n, "A", f"seed{r}") for r in range(6)], Q4, "snaps"))
    if c["N"] == 400:
        out.append(("P4a A, first quarter [0, 5000) (all at 0 deg)", "A_Q1", [os.path.join(P4D, n, "A", f"seed{r}") for r in range(6)], Q1, "snaps"))
    out.append(("P4 B, compressed u = 0.01", "B", [os.path.join(P4D, n, "B", f"seed{r}") for r in range(6)], Q4, "snaps"))
    for u in ("0.001", "0.0001"):
        dd = os.path.join(P5D, f"{n}_u{u}", "B")
        if os.path.isdir(dd): out.append((f"P5 B, compressed u = {u}", f"B{u}", [os.path.join(dd, f"seed{r}") for r in range(6)], Q4, "snaps"))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data"); ap.add_argument("--bin"); ap.add_argument("--selftest"); a = ap.parse_args()
    if a.selftest: return selftest(os.path.abspath(a.selftest))
    P6, binp = os.path.abspath(a.data), os.path.abspath(a.bin)
    print("# Pilot P6, the turned-lattice start (261012 sec. 4.7.43; decision 18, item 3), printed by experiments_gen3_p6_261010/p6_tables.py -- "
          "A PILOT, no paper use, no verdict; the questions are answered with the numbers only\n")
    print("Q4 = [1.5e4, 2e4) sigma-time of the hold (virial blocks 300-399); per compartment Z_c, |Psi_c| and the orientation (0: close-packed lines "
          "along x; +-30: along y, turned); |Psi| and orientation from P6's psi6(t) columns, for P4/P5 from their snapshots (seeds 0, 1 only)\n")
    print("## Inventory (P6)\n\n| cell | start | runs | clean | infeasible | findings |\n|---|---|---|---|---|---|")
    for c in map(R.cell, R.cells()):
        for lab, o, l in R.starts_for(c):
            st = []
            for r in range(6):
                d = os.path.join(P6, c["name"], lab, f"seed{r}")
                lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
                e = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
                st.append(("infeasible" if "[EDMD3-M4] INFEASIBLE" in lg else ("clean" if e == "0" and "clean=1" in lg else f"exit {e}")))
            print(f"| {c['name']} | {lab}{' (ADDITION)' if lab == 'zero' else ''} | 6 | {st.count('clean')} | {st.count('infeasible')} | "
                  f"{'; '.join(f'seed{r}: {s}' for r, s in enumerate(st) if s not in ('clean', 'infeasible')) or '-'} |")
    summary = {}
    for c in map(R.cell, R.cells()):
        print(f"\n## {c['name']}\n\n### Per seed (Q4; for 'P4a A, first quarter' the window [0, 5000))\n")
        print("| start | seed | Z_0 | Z_1 | abs(Psi_0) | orientation_0 [deg] | abs(Psi_1) | orientation_1 [deg] |\n|---|---|---|---|---|---|---|---|")
        S = {}
        for lab, kind, ds, w, src in starts(c, P6):
            zs, ps = [], []
            for r, d in enumerate(ds):
                z = zwin(d, w); p = psi_cols(d, w) if src == "cols" else (psi_snaps(d, w, binp) if r in (0, 1) else None)
                zs.append(z); ps.append(p)
                f = lambda q, k, j: f"{q[k][j]:.4f}" if q and q[k] else "-"
                fo = lambda q, k: f"{q[k][1]:+.2f}" if q and q[k] else "-"
                print(f"| {lab} | {r} | " + (f"{z[0]:.4f} | {z[1]:.4f}" if z else "- | -") + f" | {f(p, 0, 0)} | {fo(p, 0)} | {f(p, 1, 0)} | {fo(p, 1)} |")
            S[kind] = dict(lab=lab, z=[x for z in zs if z for x in z], psi=[p for p in ps if p])
        print("\n### Per start: mean (SE) over the compartments; orientations where known\n")
        print("| start | compartments | mean Z_c (SE) | smallest / largest Z_c | mean abs(Psi_c) | orientations known [deg] |\n|---|---|---|---|---|---|")
        for kind, v in S.items():
            m, se, n = mse(v["z"]); pa = [q[k][0] for q in v["psi"] for k in (0, 1) if q[k]]
            ors = [f"{q[k][1]:+.1f}" for q in v["psi"] for k in (0, 1) if q[k]]
            print(f"| {v['lab']} | {n} | {m:.4f} ({se:.4f}) | " + (f"{min(v['z']):.4f} / {max(v['z']):.4f}" if v["z"] else "-") +
                  f" | {(sum(pa) / len(pa)) if pa else float('nan'):.4f} | {' '.join(ors) or '-'} |")
        comp = [x for k in S if k.startswith("B") for x in S[k]["z"]]
        cut, gap = groups(comp) if c["N"] == 1600 else (None, None)
        if cut is not None:
            lo = [x for x in comp if x < cut]; hi = [x for x in comp if x > cut]
            print(f"\nthe compressed starts' compartments (all {len(comp)}): split at their largest gap {gap:.4f} (cut {cut:.4f}): the lower group "
                  f"{len(lo)} at {mse(lo)[0]:.4f} ({mse(lo)[1]:.4f}), the higher group {len(hi)} at {mse(hi)[0]:.4f} ({mse(hi)[1]:.4f})")
        summary[c["name"]] = (c, S, comp, cut)
    print("\n## Q1: is the compressed start's state the turned crystal? (the turned start's mean Z_c against the compressed starts')\n")
    print("| cell | turned start: mean Z_c (SE), n | against | compressed: mean Z_c (SE), n | difference turned - compressed (SE) | / SE |\n|---|---|---|---|---|---|")
    for name, (c, S, comp, cut) in summary.items():
        t = mse(S["turned"]["z"]) if "turned" in S else (float("nan"),) * 3
        sets = [("all compressed compartments", comp)]
        if cut is not None: sets += [("the compressed lower Z group", [x for x in comp if x < cut]), ("the compressed higher Z group", [x for x in comp if x > cut])]
        for lab, v in sets:
            m, se, n = mse(v); d, sd = t[0] - m, math.hypot(t[1], se)
            print(f"| {name} | {t[0]:.4f} ({t[1]:.4f}), {t[2]} | {lab} | {m:.4f} ({se:.4f}), {n} | {d:+.4f} ({sd:.4f}) | {d / sd:+.2f} |")
    print("\nthe compressed compartments' orientations from their snapshots (seeds 0 and 1 of each compressed start, Q4) are in the per-start tables above")
    print("\n## Q2: which orientation has the lower Z? (the 0-deg start's mean Z_c against the turned start's)\n")
    print("| cell | 0 deg: which start | mean Z_c (SE), n | turned: P6's turned start, mean Z_c (SE), n | difference 0 - turned (SE) | / SE | lower |\n|---|---|---|---|---|---|---|")
    for name, (c, S, comp, cut) in summary.items():
        t = mse(S["turned"]["z"]) if "turned" in S else (float("nan"),) * 3
        zero = [("A_Q1", S.get("A_Q1")), ("A", S.get("A"))] if c["N"] == 400 else [("zero", S.get("zero"))]
        for kind, v in zero:
            if not v: continue
            m, se, n = mse(v["z"]); d, sd = m - t[0], math.hypot(se, t[1])
            low = "0 deg" if d < 0 else ("turned" if d > 0 else "equal")
            print(f"| {name} | {v['lab']} | {m:.4f} ({se:.4f}), {n} | {t[0]:.4f} ({t[1]:.4f}), {t[2]} | {d:+.4f} ({sd:.4f}) | {d / sd:+.2f} | {low} |")
    print("\n(PILOT: no verdict, no paper use; the questions are answered with the numbers above.)")


def selftest(root):
    """the window averages and the circular orientation mean on synthetic files"""
    if os.listdir(root): sys.exit("--selftest needs an empty folder")
    ok = True
    d = os.path.join(root, "run"); os.makedirs(d)
    with open(os.path.join(d, "virial_blocks_L0_1_wallmassfactor_300_run0.csv"), "w") as f:
        f.write("t_sigma,phase,block,dt_sigma,N_0,KE_0,pair_events_0,Z_0,N_1,KE_1,pair_events_1,Z_1,pair_events,Z_global,Z_KE_weighted_mean,weighted_minus_global\n")
        for b in range(400): f.write(f"{50.0 * (b + 1):.6f},hold,{b},{50.0 if b % 2 else 50.0000014},200,200,1,{13.0 + (0.5 if b >= 300 else 0.0) + 0.001 * (b % 3)},200,200,1,{14.0},2,0,0,0\n")
    z = zwin(d, Q4); want0 = sum((13.5 + 0.001 * (b % 3)) * (50.0 if b % 2 else 50.0000014) for b in range(300, 400)) / sum((50.0 if b % 2 else 50.0000014) for b in range(300, 400))
    e = abs(z[0] - want0) < 1e-12 and abs(z[1] - 14.0) < 1e-12; ok &= e
    print(f"selftest Z window: {z} (direct {want0:.12f}, 14): {'ok' if e else '**WRONG**'}")
    with open(os.path.join(d, "psi6_t_L0_1_wallmassfactor_300_run0.csv"), "w") as f:
        f.write("t_sigma,phase,psi6_global,psi6_local,mean_neighbors,N_0,psi6_0,psi6_orient_deg_0,psi6_local_0,N_1,psi6_1,psi6_orient_deg_1,psi6_local_1\n")
        for t in range(14990, 20010):
            o0 = 29.5 if t % 2 else -29.5                       # +-29.5 alternate: the same orientation, mod 60
            f.write(f"{t:.6f},hold,0.1,0.8,5.5,200,0.8,{o0},0.85,200,0.6,{0.5 if t % 2 else -0.5},0.8\n")
    p = psi_cols(d, Q4)
    e = abs(p[0][0] - 0.8) < 1e-12 and abs(abs(p[0][1]) - 30.0) < 0.6 and abs(p[1][0] - 0.6) < 1e-12 and abs(p[1][1]) < 1e-9; ok &= e
    print(f"selftest psi columns (+-29.5 alternating -> +-30 mod 60; +-0.5 -> 0): {p}: {'ok' if e else '**WRONG**'}")
    cg = groups([1.0, 1.1, 2.0, 2.05]); e = abs(cg[0] - 1.55) < 1e-12 and abs(cg[1] - 0.9) < 1e-12; ok &= e
    print(f"selftest groups: {cg} (cut 1.55, gap 0.9): {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
