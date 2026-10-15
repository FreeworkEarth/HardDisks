#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.37; plan-author decision 16, item 2a): the print of THE PART-3A OUTLIER RUN -- information, no
verdict, no interpretation. For the run of d16_outlier_run.py (pi/8, seed 1, the hold extended to 2e4 sigma-time) it prints:
  0. the inventory: exit code, one build, one run record and clean=1;
  1. whether the first 100 hold virial blocks are identical to the original run's (experiments_gen3_p3a_261009/evidence_5d19aa7/e4/
     pi8_seed1): the CSV rows compared as text, byte for byte;
  2. d = Z_L - Z_R per 2500 sigma-time segment (8 segments of 50 blocks), with the pair-collision counts of each compartment. Z_L and
     Z_R are the plain means of the blocks' Z_0 and Z_1 (equal block lengths), their SEs the block scatter, as part 3a's table
     (p3a_tables.py, section 4) computed them; also [0, 5000) (the original's span) and [0, 2e4);
  3. the initial state per compartment, from --gen3-dump-initial: N, KE (= 1/2 sum v^2, m = 1, internal units = kT units; next to it
     the first virial block's KE_0 / KE_1), the total momentum (x, y), the lattice form (the run header [EDMD3-RUN]), and the smallest
     gaps: disk-disk (surface to surface, pairs in the same compartment), disk-outer wall (which wall), disk-divider.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d16_261010/d16_outlier_tables.py --data <data root>
  python3 experiments_gen3_d16_261010/d16_outlier_tables.py --selftest <empty folder>
"""
import argparse, csv, glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import d16_outlier_run as RUN
SEG_BLOCKS, NSEG, NORIG = 50, 8, 100


def one(d, pat):
    f = sorted(glob.glob(os.path.join(d, pat))); return f[0] if len(f) == 1 else None


def hold_rows_text(p):
    """the hold rows of a virial-block CSV as their exact text lines (without the line end)"""
    lines = open(p, "rb").read().decode().split("\n")
    return [l for l in lines[1:] if l and l.split(",")[1] == "hold"], lines[0]


def mse(v):
    m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else float("nan")
    return m, sd / math.sqrt(len(v)), sd


def seg_stats(bl):
    zl, sl, _ = mse([float(b["Z_0"]) for b in bl]); zr, sr, _ = mse([float(b["Z_1"]) for b in bl])
    return dict(t0=float(bl[0]["t_sigma"]) - float(bl[0]["dt_sigma"]), t1=float(bl[-1]["t_sigma"]), n=len(bl), zl=zl, sl=sl, zr=zr, sr=sr,
                d=zl - zr, sd=math.hypot(sl, sr), pl=sum(int(b["pair_events_0"]) for b in bl), pr=sum(int(b["pair_events_1"]) for b in bl))


def read_dump(p):
    L = open(p).read().split("\n")
    if L[0] != "gen3-initial v1": raise ValueError("dump format")
    N, NL = int(L[1].split()[1]), int(L[1].split()[3])
    W, H, R = (float.fromhex(L[k].split()[1]) for k in (2, 3, 4))
    c, th = (float.fromhex(v) for v in L[5].split()[1:3])
    P = [tuple(float.fromhex(v) for v in l.split()) for l in L[6:6 + N]]
    return dict(N=N, NL=NL, W=W, H=H, R=R, c=c, th=th, P=P)


def compartments(D):
    """per compartment (left: x < divider centre): N, KE, momentum, smallest gaps [px]"""
    out = []
    for side in ("L", "R"):
        Q = [p for p in D["P"] if (p[0] < D["c"]) == (side == "L")]
        R2, best = 2.0 * D["R"], (float("inf"), None)
        for i in range(len(Q)):
            for j in range(i + 1, len(Q)):
                g = math.hypot(Q[i][0] - Q[j][0], Q[i][1] - Q[j][1]) - R2
                if g < best[0]: best = (g, (i, j))
        walls = []
        for x, y, _, _ in Q:
            if side == "L": walls.append((x - D["R"], "outer left (x = 0)"))
            else: walls.append((D["W"] - D["R"] - x, "outer right (x = boxW)"))
            walls += [(y - D["R"], "bottom (y = 0)"), (D["H"] - D["R"] - y, "top (y = boxH)")]
        div = [((D["c"] - D["th"] / 2) - (x + D["R"])) if side == "L" else ((x - D["R"]) - (D["c"] + D["th"] / 2)) for x, _, _, _ in Q]
        out.append(dict(side=side, N=len(Q), KE=0.5 * sum(vx * vx + vy * vy for _, _, vx, vy in Q), px=sum(q[2] for q in Q), py=sum(q[3] for q in Q),
                        pair=best[0], wall=min(walls), div=min(div), xs=(min(q[0] for q in Q), max(q[0] for q in Q))))
    return out


def lattice(lg):
    m = re.search(r"^\[EDMD3-RUN\] .*? left=\[(.*?)\] right=\[(.*?)\]", lg, re.M)
    return (m.group(1), m.group(2)) if m else ("-", "-")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data"); ap.add_argument("--selftest"); a = ap.parse_args()
    if a.selftest: return selftest(os.path.abspath(a.selftest))
    d = os.path.join(os.path.abspath(a.data), RUN.NAME)
    print("# Item 2a (261012 sec. 4.7.37; decision 16): the part-3a outlier pi/8 seed 1, the hold extended to 2e4 sigma-time -- printed by "
          "experiments_gen3_d16_261010/d16_outlier_tables.py; information, no verdict\n")
    lg = open(os.path.join(d, "run.log"), errors="ignore").read()
    ec = open(os.path.join(d, "exit_code.txt")).read().strip()
    recs = re.findall(r"^\[EDMD3-HEALTH\] (.*)$", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
    kv = dict(re.findall(r" (\w+)=(\S+)", " " + recs[0].split(": ", 1)[-1])) if len(recs) == 1 else {}
    print(f"## 0. Inventory\n\nexit {ec}; builds {nb}; run records {len(recs)}; clean={kv.get('clean')}; build {kv.get('build')} target {kv.get('target')}; "
          f"hash {kv.get('hash')}; t_end {kv.get('t_end')}; run_s {kv.get('run_s')}\n")
    rows_new, hdr_new = hold_rows_text(one(d, "virial_blocks_*_run0.csv"))
    rows_old, hdr_old = hold_rows_text(one(RUN.ORIG, "virial_blocks_*_run0.csv"))
    same = [i for i in range(NORIG) if i < len(rows_new) and rows_new[i] == rows_old[i]]
    first_diff = next((i for i in range(NORIG) if i >= len(rows_new) or rows_new[i] != rows_old[i]), None)
    print("## 1. The first 100 hold virial blocks against the original run's\n")
    print(f"original: {len(rows_old)} hold blocks; this run: {len(rows_new)} hold blocks; header identical: {hdr_new == hdr_old}")
    print(f"rows byte-identical: {len(same)} of {NORIG}" + (" -- IDENTICAL" if len(same) == NORIG else f" -- **NOT IDENTICAL**, first differing block {first_diff}"))
    if first_diff is not None and first_diff < len(rows_new):
        print(f"  original block {first_diff}: {rows_old[first_diff]}\n  this run block {first_diff}: {rows_new[first_diff]}")
    bl = [b for b in csv.DictReader(open(one(d, "virial_blocks_*_run0.csv"))) if b["phase"] == "hold"]
    print(f"\n## 2. d = Z_L - Z_R per 2500 sigma-time segment (50 blocks of 50 sigma-time each); SE from the block scatter (as part 3a)\n")
    print("| segment | sigma-time | blocks | Z_L (SE) | Z_R (SE) | d = Z_L - Z_R (SE) | d / SE | pair collisions L | pair collisions R | R - L |\n"
          "|---|---|---|---|---|---|---|---|---|---|")
    segs = [seg_stats(bl[SEG_BLOCKS * j: SEG_BLOCKS * (j + 1)]) for j in range(NSEG) if len(bl) >= SEG_BLOCKS * (j + 1)]
    for j, s in enumerate(segs):
        print(f"| {j + 1} | [{s['t0']:.0f}, {s['t1']:.0f}) | {s['n']} | {s['zl']:.5f} ({s['sl']:.5f}) | {s['zr']:.5f} ({s['sr']:.5f}) | "
              f"{s['d']:+.5f} ({s['sd']:.5f}) | {s['d'] / s['sd']:+.2f} | {s['pl']} | {s['pr']} | {s['pr'] - s['pl']:+d} |")
    for lab, sub in (("[0, 5000) -- the original's span", bl[:NORIG]), ("[0, 2e4) -- the whole hold", bl)):
        s = seg_stats(sub)
        print(f"| {lab} | [{s['t0']:.0f}, {s['t1']:.0f}) | {s['n']} | {s['zl']:.5f} ({s['sl']:.5f}) | {s['zr']:.5f} ({s['sr']:.5f}) | "
              f"{s['d']:+.5f} ({s['sd']:.5f}) | {s['d'] / s['sd']:+.2f} | {s['pl']} | {s['pr']} | {s['pr'] - s['pl']:+d} |")
    if segs:
        m, se, sd = mse([s["d"] for s in segs])
        print(f"\nthe {len(segs)} segment values of d: mean {m:+.5f}, SD {sd:.5f}, SD / sqrt({len(segs)}) {se:.5f}")
    D = read_dump(os.path.join(d, "initial_state.txt")); sg = 2.0 * D["R"]
    C = compartments(D); latl, latr = lattice(lg)
    print(f"\n## 3. The initial state per compartment (--gen3-dump-initial; px, sigma = 2 x radius = {sg:g} px)\n")
    print(f"box {D['W']:.10g} x {D['H']:.10g} px ({D['W'] / sg:.10g} x {D['H'] / sg:.10g} sigma); divider centre {D['c']:.10g} px, thickness "
          f"{D['th']:.6g} px; N {D['N']}, N_L (header) {D['NL']}\n")
    ke0 = (float(bl[0]["KE_0"]), float(bl[0]["KE_1"]))
    print("| compartment | N | KE (dump) | KE (first virial block) | momentum x | momentum y | smallest disk-disk gap [px] ([sigma]) | "
          "smallest disk-wall gap [px] ([sigma]), wall | smallest disk-divider gap [px] ([sigma]) | disk centres x from / to [px] |\n|---|---|---|---|---|---|---|---|---|---|")
    for k, c in enumerate(C):
        print(f"| {c['side']} | {c['N']} | {c['KE']:.12g} | {ke0[k]:.12g} | {c['px']:+.6e} | {c['py']:+.6e} | {c['pair']:.6g} ({c['pair'] / sg:.6g}) | "
              f"{c['wall'][0]:.6g} ({c['wall'][0] / sg:.6g}), {c['wall'][1]} | {c['div']:.6g} ({c['div'] / sg:.6g}) | {c['xs'][0]:.6g} / {c['xs'][1]:.6g} |")
    print(f"\nthe lattice form (run header [EDMD3-RUN]):\n  left:  [{latl}]\n  right: [{latr}]")
    print("\n(Information only; no verdict, no interpretation.)")


def selftest(root):
    """a synthetic dump and block file: the gap, KE, momentum and segment arithmetic against direct values"""
    if os.listdir(root): sys.exit("--selftest needs an empty folder")
    ok = True
    D = dict(N=4, NL=2, W=100.0, H=50.0, R=2.0, c=50.0, th=1.0,
             P=[(10.0, 12.0, 1.0, 0.0), (15.0, 12.0, -1.0, 2.0), (60.0, 40.0, 0.5, 0.5), (90.0, 47.0, 0.0, -1.0)])
    C = compartments(D)
    want = [dict(N=2, KE=0.5 * (1 + 1 + 4), px=0.0, py=2.0, pair=1.0, wall=(8.0, "outer left (x = 0)"), div=(49.5 - 17.0)),
            dict(N=2, KE=0.5 * (0.5 + 1.0), px=0.5, py=-0.5, pair=math.hypot(30, 7) - 4.0, wall=(1.0, "top (y = boxH)"), div=(58.0 - 50.5))]
    for c, w in zip(C, want):
        e = all(abs(c[k] - w[k]) < 1e-12 for k in ("N", "KE", "px", "py", "pair", "div")) and c["wall"] == w["wall"]; ok &= e
        print(f"selftest compartment {c['side']}: N {c['N']} KE {c['KE']} p ({c['px']}, {c['py']}) pair {c['pair']:.6f} wall {c['wall']} div {c['div']}: "
              f"{'ok' if e else '**WRONG**'}")
    p = os.path.join(root, "vb.csv")
    with open(p, "w") as f:
        f.write("t_sigma,phase,block,dt_sigma,N_0,KE_0,pair_events_0,Z_0,N_1,KE_1,pair_events_1,Z_1,pair_events,Z_global,Z_KE_weighted_mean,weighted_minus_global\n")
        for b in range(8): f.write(f"{50.0 * (b + 1):.6f},hold,{b},50.0000007,200,200,{100 + b},{2.8 + 0.01 * b},200,200,{90 + b},{2.79 + 0.02 * (b % 2)},0,0,0,0\n")
        f.write("401.000000,record,8,1,200,200,1,9,200,200,1,9,0,0,0,0\n")
    rows, _ = hold_rows_text(p); bl = [b for b in csv.DictReader(open(p)) if b["phase"] == "hold"]
    s = seg_stats(bl[:4]); zl = sum(2.8 + 0.01 * b for b in range(4)) / 4; zr = sum(2.79 + 0.02 * (b % 2) for b in range(4)) / 4
    e = len(rows) == 8 and abs(s["d"] - (zl - zr)) < 1e-12 and s["pl"] == 406 and s["pr"] == 366 and abs(s["t0"]) < 1e-6 and abs(s["t1"] - 200) < 1e-6; ok &= e
    print(f"selftest segment: rows {len(rows)}, d {s['d']:.6f} (direct {zl - zr:.6f}), pairs {s['pl']} / {s['pr']} (direct 406 / 366): {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
