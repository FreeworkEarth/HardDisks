#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.38; plan-author decision 16, item 3): the tables of PILOT P5 (the compression speed) -- a pilot,
no paper use, no verdict. SE = SD / sqrt(n) over seeds (the seed-to-seed scatter); differences carry sqrt(SE_1^2 + SE_2^2).
Per cell (p5_run.cells()), the plan author's list:
  - the inventory (every run present, exit 0, one run record clean=1, one build, one [EDMD3-M4B] line) and the builds;
  - the compression: psi6 at its end (global, local, neighbours; [EDMD3-M4B]), its length in sigma-time, its wall time, and the KE
    per compartment before the rescale (as KE / N_c = kT, the compartments' mean), next to P4's start B at u = 0.01;
  - psi6 (global) and Z_global in the four quarters of the hold (mean (SE) over seeds), computed as P4 computed them
    (p4_tables.run(): per seed the mean of the psi6(t) rows of the hold in the quarter and the time average of Z_global over its 100
    virial blocks), next to P4's rows for start A and for start B (u = 0.01) of the same N and eta. For N = 400, P4a ran 1e5
    sigma-time: its rows are its first 2e4 sigma-time in the same quarters (5000 sigma-time each) and its last tenth [9e4, 1e5);
  - the differences in the last quarter, P5 - A and P5 - B(0.01), with their SEs and the ratio to the SE (for N = 400 also against
    P4a's last tenth, as information);
  - information: the share of P4's A - B(0.01) difference that P5 closes in the last quarter, (P5 - B) / (A - B), with a first-order
    SE from the three independent seed sets.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p5_261010/p5_tables.py --data <P5 data root> [--p4 <P4 data root>]
  python3 experiments_gen3_p5_261010/p5_tables.py --selftest <empty folder>
"""
import argparse, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HS, "experiments_gen3_p4_261010"))
import p5_run as R
import p4_tables as P4T
P4DATA = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248")
mse = P4T.mse


def m4b(d):
    """the [EDMD3-M4B] numbers: compression length [sigma-time], wall time [s], KE / N_c before the rescale (both compartments),
    psi6 at the end of the compression (global, local, neighbours), the piston speed"""
    p = os.path.join(d, "run.log")
    lg = open(p, errors="ignore").read() if os.path.exists(p) else ""
    out = {}
    m = re.search(r"^\[EDMD3-M4B\] L0=.*? u=(\S+) sigma/sigma-time;.*?before the rescale (\S+) (\S+), after (\S+) (\S+) .*?"
                  r"\(melt (\S+) \+ compression (\S+)\);.*?N_L=(\d+) N_R=(\d+) .*; ([0-9.]+) s$", lg, re.M)
    if m:
        u, kl, kr, al, ar, melt, comp, nl, nr, ws = m.groups()
        out.update(u=float(u), kt_before=(float(kl) / int(nl), float(kr) / int(nr)), kt_after=(float(al) / int(nl), float(ar) / int(nr)),
                   melt=float(melt), comp_t=float(comp), wall_s=float(ws))
    m = re.search(r"^\[EDMD3-M4B\] psi6 at the end of the compression \(the carried state\): global (\S+) local (\S+) neighbours (\S+)", lg, re.M)
    if m: out.update(psi_g=float(m.group(1)), psi_l=float(m.group(2)), nb=float(m.group(3)))
    return out


def ratio(x, sx, a, sa, b, sb):
    """(x - b) / (a - b) and its first-order SE for independent x, a, b"""
    g = a - b
    if g == 0: return float("nan"), float("nan")
    f = (x - b) / g
    return f, math.sqrt((sx / g) ** 2 + ((x - b) * sa / g ** 2) ** 2 + ((x - a) * sb / g ** 2) ** 2)


def p5_runs(O, c):
    out = []
    for r in range(R.NSEED):
        d = os.path.join(O, c["name"], "B", f"seed{r}")
        v = P4T.run(d, c["hold"] / 60.0, 4); v.update(m4b(d)); out.append(v)
    return out


def p4_rows(P4O, c):
    """P4's start A and B (u = 0.01) rows for the cell's N and eta: four quarters of 2e4 sigma-time (N = 1600: the whole hold; N = 400:
    P4a's first 2e4 sigma-time, as quarters 1-4 of its twenty 5000-sigma-time segments), and for N = 400 also P4a's last tenth"""
    rows = {}
    for lab in ("A", "B"):
        rs = []
        for r in range(R.NSEED):
            d = os.path.join(P4O, c["p4name"], lab, f"seed{r}")
            if c["N"] == 1600: v = P4T.run(d, 1200000 / 60.0, 4)
            else:
                v = P4T.run(d, 6000000 / 60.0, 20)
                if v["ok"]:
                    t10 = P4T.run(d, 6000000 / 60.0, 10)
                    v["psi"], v["z"], v["psi_t10"], v["z_t10"] = v["psi"][:4], v["z"][:4], t10["psi"][9], t10["z"][9]
            v.update(m4b(d)); rs.append(v)
        rows[lab] = rs
    return rows


def fmt(m, s): return f"{m:.4f} ({s:.4f})"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data"); ap.add_argument("--p4", default=P4DATA); ap.add_argument("--selftest")
    a = ap.parse_args()
    if a.selftest: return selftest(os.path.abspath(a.selftest))
    O, P4O = os.path.abspath(a.data), os.path.abspath(a.p4)
    C = [R.cell(c) for c in R.cells()]
    print("# Pilot P5, the compression speed (261012 sec. 4.7.38; decision 16, item 3), printed by experiments_gen3_p5_261010/p5_tables.py "
          "-- A PILOT, no paper use, no verdict\n")
    print("Start B (M4b: the lattice at eta_0 = 0.60 in a longer box, 500 sigma-time melt, the pistons at u, KE rescaled per compartment) "
          "with u = 0.001 and 0.0001, then the divider held 2e4 sigma-time; 6 independent seeds per cell (base 20261226). Next to it P4 "
          "(base 20261225): start A (the M4 lattice at the target eta) and start B at u = 0.01. SE = SD / sqrt(n) over seeds.\n")
    V, P, builds = {}, {}, set()
    print("## Inventory\n\n| cell | runs | clean (exit 0, one record clean=1, one build) | with an [EDMD3-M4B] line | findings |\n|---|---|---|---|---|")
    for c in C:
        rs = p5_runs(O, c); V[c["name"]] = rs
        for v in rs:
            if v["ok"]: builds.add(v["build"])
        bad = [f"seed{r}: {v['why']}" for r, v in enumerate(rs) if not v["ok"]]
        print(f"| {c['name']} | {R.NSEED} | {sum(v['ok'] for v in rs)} | {sum('comp_t' in v for v in rs)} | {'; '.join(bad) if bad else '-'} |")
        P[c["name"]] = p4_rows(P4O, c)
    print(f"\nP5 builds: {sorted(builds)}; P4 builds: {sorted({v['build'] for c in C for lab in 'AB' for v in P[c['name']][lab] if v['ok']})}")
    print("\n## The compression (start B): psi6 at its end, its length, its wall time, kT before the rescale (mean (SE) over seeds)\n")
    print("| cell | u [sigma/sigma-time] | n | compression [sigma-time] | psi6 global at its end | psi6 local | neighbours | wall time [s] | "
          "kT before the rescale (KE / N_c, both compartments) | events/s (hold) |\n|---|---|---|---|---|---|---|---|---|---|")
    for c in C:
        for lab, rs in ((f"P5 B", V[c["name"]]), ("P4 B", P[c["name"]]["B"])):
            ok = [v for v in rs if v["ok"] and "comp_t" in v]
            if not ok: print(f"| {c['name']} ({lab}) | - | 0 | - | - | - | - | - | - | - |"); continue
            kt = [x for v in ok for x in v["kt_before"]]; ws, sws, _ = mse([v["wall_s"] for v in ok])
            print(f"| {c['name'] if lab == 'P5 B' else c['p4name']} ({lab}) | {ok[0]['u']:g} | {len(ok)} | {mse([v['comp_t'] for v in ok])[0]:.1f} | "
                  f"{fmt(*mse([v['psi_g'] for v in ok])[:2])} | {fmt(*mse([v['psi_l'] for v in ok])[:2])} | {fmt(*mse([v['nb'] for v in ok])[:2])} | "
                  f"{ws:.1f} ({sws:.1f}) | {fmt(*mse(kt)[:2])} | {mse([v['evs'] for v in ok])[0]:.3g} |")
    print("\n## psi6 (global) and Z_global in the four quarters of the hold (mean (SE) over seeds)\n")
    print("Q1-Q4 = [0, 5000), [5000, 1e4), [1e4, 1.5e4), [1.5e4, 2e4) sigma-time of the hold. For N = 400 P4's rows are P4a's first 2e4 "
          "sigma-time; its last tenth [9e4, 1e5) is the column 'P4a T10'.\n")
    print("| cell | start | n | psi6 Q1 | psi6 Q2 | psi6 Q3 | psi6 Q4 | Z Q1 | Z Q2 | Z Q3 | Z Q4 | psi6 P4a T10 | Z P4a T10 |\n" + "|---" * 13 + "|")
    for c in C:
        for lab, rs in ((f"B, u = {c['u']:g} (P5)", V[c["name"]]), ("A (P4)", P[c["name"]]["A"]), ("B, u = 0.01 (P4)", P[c["name"]]["B"])):
            ok = [v for v in rs if v["ok"]]
            t10 = (f" | {fmt(*mse([v['psi_t10'] for v in ok])[:2])} | {fmt(*mse([v['z_t10'] for v in ok])[:2])}" if "P4" in lab and c["N"] == 400
                   else " | - | -")
            print(f"| {c['name']} | {lab} | {len(ok)} | " + " | ".join(fmt(*mse([v['psi'][j] for v in ok])[:2]) for j in range(4)) + " | "
                  + " | ".join(fmt(*mse([v['z'][j] for v in ok])[:2]) for j in range(4)) + t10 + " |")
    print("\n## The differences in the last quarter (Q4): P5 - A and P5 - B(0.01), with the seed-scatter SE\n")
    print("| cell | against | psi6 difference (SE) | / SE | Z difference (SE) | / SE |\n|---|---|---|---|---|---|")
    for c in C:
        ok5 = [v for v in V[c["name"]] if v["ok"]]
        p5, sp5, _ = mse([v["psi"][3] for v in ok5]); z5, sz5, _ = mse([v["z"][3] for v in ok5])
        refs = [("A (P4), Q4", "A", "psi", 3), ("B u = 0.01 (P4), Q4", "B", "psi", 3)]
        if c["N"] == 400: refs += [("A (P4a), last tenth", "A", "t10", None), ("B u = 0.01 (P4a), last tenth", "B", "t10", None)]
        for lab, st, kind, j in refs:
            ok = [v for v in P[c["name"]][st] if v["ok"]]
            pr, spr, _ = mse([v["psi"][j] if kind == "psi" else v["psi_t10"] for v in ok]); zr, szr, _ = mse([v["z"][j] if kind == "psi" else v["z_t10"] for v in ok])
            dp, sdp, dz, sdz = p5 - pr, math.hypot(sp5, spr), z5 - zr, math.hypot(sz5, szr)
            print(f"| {c['name']} | {lab} | {dp:+.4f} ({sdp:.4f}) | {dp / sdp:+.2f} | {dz:+.4f} ({sdz:.4f}) | {dz / sdz:+.2f} |")
    print("\n## Information: the share of P4's A - B(0.01) difference that P5 closes in Q4, (P5 - B) / (A - B) (0: at B(0.01); 1: at A)\n")
    print("| cell | psi6: A - B(0.01) (SE) | psi6: share (SE) | Z: A - B(0.01) (SE) | Z: share (SE) |\n|---|---|---|---|---|")
    for c in C:
        ok5 = [v for v in V[c["name"]] if v["ok"]]; okA = [v for v in P[c["name"]]["A"] if v["ok"]]; okB = [v for v in P[c["name"]]["B"] if v["ok"]]
        cols = []
        for key in ("psi", "z"):
            x, sx, _ = mse([v[key][3] for v in ok5]); aa, sa, _ = mse([v[key][3] for v in okA]); bb, sb, _ = mse([v[key][3] for v in okB])
            f, sf = ratio(x, sx, aa, sa, bb, sb); cols += [f"{aa - bb:+.4f} ({math.hypot(sa, sb):.4f})", f"{f:+.2f} ({sf:.2f})"]
        print(f"| {c['name']} | " + " | ".join(cols) + " |")
    print("\n### Per seed (P5): psi6 at the end of the compression and at the end of the hold; Z in Q4\n")
    print("| cell | seed | psi6 at the end of the compression | psi6 Q4 | psi6 at the end of the hold | Z Q4 | compression wall time [s] |\n|---|---|---|---|---|---|---|")
    for c in C:
        for r, v in enumerate(V[c["name"]]):
            if not v["ok"]: print(f"| {c['name']} | {r} | - | - | - | - | - |"); continue
            print(f"| {c['name']} | {r} | {v.get('psi_g', float('nan')):.4f} | {v['psi'][3]:.4f} | {v['psi_end']:.4f} | {v['z'][3]:.4f} | {v.get('wall_s', float('nan')):.1f} |")
    print("\n(PILOT: no verdict, no paper use.)")


def selftest(root):
    """the [EDMD3-M4B] parser on a real P4 line; the share and its SE against direct values; P4's per-run parser is P4's (self-tested)"""
    if os.listdir(root): sys.exit("--selftest needs an empty folder")
    ok = True
    d = os.path.join(root, "run"); os.makedirs(d)
    open(os.path.join(d, "run.log"), "w").write(
        "[EDMD3-M4B] L0=21.2 M=300 run=0 seed=2628733938: eta_0=0.59999999999999998 (seeded 0.59999999999999998; target 0.73999999999999999) "
        "T_melt=500 sigma-time u=0.01 sigma/sigma-time; long box 1304.6370614359173 x 960 px, divider 652.31853071795865 px, pistons from 24 and "
        "1280.6370614359173 px; L0 requested 21.226977389120222 sigma, achieved left 21.226977389120325 right 21.226977389121569 sigma (minus "
        "requested 1.03e-13, 1.35e-12); pistons stopped at 142.87107337907088 and 1161.7659880568763 px; piston work left 5144.2227064016761 right "
        "5101.9783298209732; KE left/right at the start 799.99999999999807 799.99999999999955, before the rescale 5944.2227064016424 "
        "5901.9783298209531, after 799.99999999999989 799.9999999999992 (after/N_s kT - 1: -1.11e-16 -9.99e-16); KE gain 10246.201036222597 minus "
        "piston work 10246.20103622265 = -5.28e-11; ledger E - E0 - W -6.18e-11, scale 3.22e-07, ratio 0.000192; end of the compression "
        "t=995.29613907947214 sigma-time (melt 500 + compression 495.29613907947214); events 20464302 (pair 17162941 wall 352667 cross 2342404 div "
        "303504 piston 302780 band 6); clean=1 overlap_repair=0 wall_overdue=0 past_event=0 clamp_repair=0 cell_repair=0 grid_escape=0 stagnation=0 "
        "obj_overlap_repair=0 body_findings=0 local_findings=0 full_findings=0 hash=7a9ca3ef047221e7; contact audit 0 events, max gap [px] pair 0 "
        "wall 0 divider 0 piston 0; carry S=142.87107337907332 px; N_L=800 N_R=800 (moved 0); smallest gap after the carry [px]: walls 0.0148, "
        "divider 0.0263; 30.2 s\n"
        "[EDMD3-M4B] psi6 at the end of the compression (the carried state): global 0.733714 local 0.861015 neighbours 5.6988\n")
    v = m4b(d)
    want = dict(u=0.01, comp_t=495.29613907947214, melt=500.0, wall_s=30.2, psi_g=0.733714, psi_l=0.861015, nb=5.6988)
    e = all(abs(v[k] - w) < 1e-12 for k, w in want.items()) and abs(v["kt_before"][0] - 5944.2227064016424 / 800) < 1e-12 \
        and abs(v["kt_before"][1] - 5901.9783298209531 / 800) < 1e-12 and abs(v["kt_after"][0] - 1.0) < 1e-12; ok &= e
    print(f"selftest [EDMD3-M4B] parser: {v}: {'ok' if e else '**WRONG**'}")
    x, sx, a, sa, b, sb = 0.73, 0.004, 0.75, 0.003, 0.71, 0.005
    f, sf = ratio(x, sx, a, sa, b, sb); h = 1e-7
    fx = (ratio(x + h, sx, a, sa, b, sb)[0] - f) / h; fa = (ratio(x, sx, a + h, sa, b, sb)[0] - f) / h; fb = (ratio(x, sx, a, sa, b + h, sb)[0] - f) / h
    sd = math.sqrt((fx * sx) ** 2 + (fa * sa) ** 2 + (fb * sb) ** 2)
    e = abs(f - 0.5) < 1e-12 and abs(sf - sd) < 1e-6; ok &= e
    print(f"selftest share: {f:.6f} (direct 0.5), SE {sf:.6f} (numerical derivatives {sd:.6f}): {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
