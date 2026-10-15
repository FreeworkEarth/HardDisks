#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (261012 sec. 4.7.34; plan-author decision 15, item 2): the tables of PILOT P4 -- a pilot, no paper use, no verdict.
Starts A (the M4 lattice at the target density) and B (M4b, the compressed fluid); 6 seeds per start and cell; SE = SD / sqrt(n) over
seeds (the seed-to-seed scatter).
PART a (N = 400, eta 0.740, the divider held 1e5 sigma-time), the plan author's list:
  - psi6 (global) and Z_global in the ten tenths of the hold (tenth j = [1e4 (j - 1), 1e4 j) sigma-time), per start: per seed the mean of
    the psi6(t) rows of the hold in the tenth and the time average of Z_global over the tenth's 200 virial blocks; mean (SE) over seeds;
  - per seed psi6 and Z in the last tenth;
  - A - B in the last tenth with its seed-scatter SE (sqrt(SE_A^2 + SE_B^2)), for psi6 and Z, and the label (information, as in P3):
    "converged" if A - B is within 2 SE in both;
  - for each start-A seed the time at which psi6 first stays below 0.8 for 500 sigma-time -- the earliest hold row t0 with psi6 < 0.8 in
    every row of [t0, t0 + 500] -- or "never";
  - A's lattice form (orientation, rows, n_s x n_r, vacancies, spacings, gap) from the run header [EDMD3-RUN].
PART b (N = 1600 at eta 0.740 and 0.780, held 2e4 sigma-time): P3's tables (p3_tables.py): psi6 and Z_global in the four quarters (mean
(SE) over seeds), A - B in the last quarter with its SE and the label, the per-seed psi6 at the end of the hold; with A's lattice form.
Both parts: the inventory (every run present, exit 0, one record clean=1, one build) first, events/s (the run proper), and for B the
compression's wall time and psi6 at its end ([EDMD3-M4B]).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p4_261010/p4_tables.py --data <data root>
  python3 experiments_gen3_p4_261010/p4_tables.py --selftest <empty folder>    synthetic run folders: the parsers against direct sums
"""
import argparse, csv, glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p4_run as R
BLOCK = 50.0


def mse(v):
    v = [x for x in v if x is not None and math.isfinite(x)]
    if len(v) < 2: return (v[0] if v else float("nan")), float("nan"), len(v)
    m = sum(v) / len(v); return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) / math.sqrt(len(v)), len(v)


def first_below(rows, thr=0.8, span=500.0):
    """the earliest hold time t0 with psi6 < thr in every psi6(t) row of [t0, t0 + span]; None if never"""
    start = None
    for t, p in rows:
        if p < thr:
            if start is None: start = t
            if t - start >= span - 1e-6: return start
        else: start = None
    return None


def lattice(lg):
    m = re.search(r"^\[EDMD3-RUN\] .*? left=\[(.*?)\] right=\[(.*?)\]", lg, re.M)
    if not m: return None
    def f(s):
        g = re.search(r"orient (\S+) rows (\S+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) g (\S+) s (\S+) n_vac (\d+)", s)
        return f"{g.group(1)}, {g.group(2)} rows, {g.group(3)} x {g.group(4)}, {g.group(9)} vacancies, a_s {g.group(5)} a_r {g.group(6)} px, gap s {g.group(8)} px" if g else s
    return f"left: {f(m.group(1))}; right: {f(m.group(2))}"


def run(d, hold, nseg):
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    ec = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
    recs = re.findall(r"^\[EDMD3-HEALTH\] (.*)$", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
    kv = dict(re.findall(r" (\w+)=(\S+)", " " + recs[0].split(": ", 1)[-1])) if len(recs) == 1 else {}
    ok = ec == "0" and nb == 1 and kv.get("clean") == "1"
    v = dict(ok=ok, why="" if ok else f"exit {ec}, builds {nb}, records {len(recs)}, clean {kv.get('clean')}", build=(kv.get("build"), kv.get("target")),
             lattice=lattice(lg))
    if not ok: return v
    v["evs"] = sum(int(kv[k]) for k in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")) / float(kv["run_s"])
    ps = [(float(r["t_sigma"]), float(r["psi6_global"])) for r in csv.DictReader(open(glob.glob(os.path.join(d, "psi6_t_*_run0.csv"))[0])) if r["phase"] == "hold"]
    vb = [r for r in csv.DictReader(open(glob.glob(os.path.join(d, "virial_blocks_*_run0.csv"))[0])) if r["phase"] == "hold"]
    seg, nb_seg = hold / nseg, int(round(hold / nseg / BLOCK))
    v["psi"], v["z"] = [], []
    for j in range(nseg):
        x = [p for t, p in ps if seg * j <= t < seg * (j + 1)]
        v["psi"].append(sum(x) / len(x) if x else float("nan"))
        b = [r for r in vb if nb_seg * j <= int(r["block"]) < nb_seg * (j + 1)]
        dt = [float(r["dt_sigma"]) for r in b]
        v["z"].append(sum(float(r["Z_global"]) * w for r, w in zip(b, dt)) / sum(dt) if len(b) == nb_seg else float("nan"))
    v["psi_end"], v["t_end"] = ps[-1][1], ps[-1][0]
    v["below"] = first_below(ps)
    m = re.search(r"^\[EDMD3-M4B\] .*; ([0-9.]+) s$", lg, re.M); v["comp_s"] = float(m.group(1)) if m else None
    m = re.search(r"^\[EDMD3-M4B\] psi6 at the end of the compression \(the carried state\): global (\S+)", lg, re.M)
    v["psi_comp"] = float(m.group(1)) if m else None
    return v


def tables(O, cells_, nseg, seglab):
    V, builds = {}, set()
    print("| cell | start | runs | clean | findings |\n|---|---|---|---|---|")
    for c in cells_:
        for lab, _ in R.STARTS:
            rs = [run(os.path.join(O, c["name"], lab, f"seed{r}"), c["hold"] / 60.0, nseg) for r in range(R.NSEED)]
            V[(c["name"], lab)] = rs
            for v in rs:
                if v["ok"]: builds.add(v["build"])
            bad = [f"seed{r}: {v['why']}" for r, v in enumerate(rs) if not v["ok"]]
            print(f"| {c['name']} | {lab} | {R.NSEED} | {sum(v['ok'] for v in rs)} | {'; '.join(bad) if bad else '-'} |")
    print(f"\nbuilds: {sorted(builds)}")
    print("\n### A's lattice form (from the run header [EDMD3-RUN])\n")
    for c in cells_:
        L = {v["lattice"] for v in V[(c["name"], "A")] if v["lattice"]}
        print(f"- {c['name']}, start A: " + (" | ".join(sorted(L)) if L else "-") + (f" (the same in all {R.NSEED} seeds)" if len(L) == 1 else " (**differs between seeds**)"))
    print(f"\n### psi6 (global) and Z_global in the {seglab}s of the hold (mean (SE) over seeds)\n")
    hdr = [f"psi6 {seglab[0].upper()}{j + 1}" for j in range(nseg)] + [f"Z {seglab[0].upper()}{j + 1}" for j in range(nseg)]
    print("| cell | start | n | " + " | ".join(hdr) + " | events/s (run) | B: compression wall time [s]; psi6 at its end (mean) |\n" + "|---" * (len(hdr) + 5) + "|")
    for c in cells_:
        for lab, _ in R.STARTS:
            ok = [v for v in V[(c["name"], lab)] if v["ok"]]
            ps = [mse([v["psi"][j] for v in ok]) for j in range(nseg)]; zs = [mse([v["z"][j] for v in ok]) for j in range(nseg)]
            comp = "-" if lab == "A" else f"{mse([v['comp_s'] for v in ok])[0]:.1f}; {mse([v['psi_comp'] for v in ok])[0]:.4f}"
            print(f"| {c['name']} | {lab} | {len(ok)} | " + " | ".join(f"{m:.4f} ({s:.4f})" for m, s, _ in ps) + " | "
                  + " | ".join(f"{m:.4f} ({s:.4f})" for m, s, _ in zs) + f" | {mse([v['evs'] for v in ok])[0]:.3g} | {comp} |")
    last = f"{seglab[0].upper()}{nseg}"
    print(f"\n### A - B in the last {seglab} ({last}), and the label (information: 'converged' if within 2 SE in both psi6 and Z)\n")
    print("| cell | psi6 A - B (SE) | / SE | Z A - B (SE) | / SE | label |\n|---|---|---|---|---|---|")
    for c in cells_:
        okA = [v for v in V[(c["name"], "A")] if v["ok"]]; okB = [v for v in V[(c["name"], "B")] if v["ok"]]
        pa, sa, _ = mse([v["psi"][-1] for v in okA]); pb, sb, _ = mse([v["psi"][-1] for v in okB])
        za, sza, _ = mse([v["z"][-1] for v in okA]); zb, szb, _ = mse([v["z"][-1] for v in okB])
        dp, sp = pa - pb, math.hypot(sa, sb); dz, sz = za - zb, math.hypot(sza, szb)
        lab = "converged" if abs(dp) < 2 * sp and abs(dz) < 2 * sz else "not converged"
        print(f"| {c['name']} | {dp:+.4f} ({sp:.4f}) | {dp / sp:+.2f} | {dz:+.4f} ({sz:.4f}) | {dz / sz:+.2f} | {lab} |")
    return V


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data"); ap.add_argument("--selftest"); a = ap.parse_args()
    if a.selftest: return selftest(os.path.abspath(a.selftest))
    O = os.path.abspath(a.data); C = [R.cell(c) for c in R.cells()]
    print("# Pilot P4 (261012 sec. 4.7.34; decision 15, item 2), printed by experiments_gen3_p4_261010/p4_tables.py -- A PILOT, no paper use\n")
    print("A = the M4 lattice at the target density; B = M4b (the lattice at eta_0 = 0.60 in a longer box, 500 sigma-time held, the pistons at "
          "0.01 sigma per sigma-time, KE rescaled); 6 seeds per start and cell (independent seeds); SE = SD / sqrt(n) over seeds\n")
    print("## Part a: N = 400, eta 0.740, the divider held 1e5 sigma-time\n")
    Ca = [c for c in C if c["part"] == "a"]
    V = tables(O, Ca, 10, "tenth")
    print("\n### Per seed: psi6 and Z_global in the last tenth; psi6 at the end of the hold; for start A the first time psi6 stays below 0.8 for "
          "500 sigma-time\n")
    print("| cell | start | seed | psi6 T10 | Z T10 | psi6 at the end | first stays below 0.8 for 500 sigma-time [sigma-time] |\n|---|---|---|---|---|---|---|")
    for c in Ca:
        for lab, _ in R.STARTS:
            for r, v in enumerate(V[(c["name"], lab)]):
                if not v["ok"]: print(f"| {c['name']} | {lab} | {r} | - | - | - | - |"); continue
                fb = ("never" if v["below"] is None else f"{v['below']:.0f}") if lab == "A" else "(not asked for start B)"
                print(f"| {c['name']} | {lab} | {r} | {v['psi'][-1]:.4f} | {v['z'][-1]:.4f} | {v['psi_end']:.4f} | {fb} |")
    print("\n## Part b: N = 1600, eta 0.740 and 0.780, the divider held 2e4 sigma-time (P3's tables)\n")
    Cb = [c for c in C if c["part"] == "b"]
    Vb = tables(O, Cb, 4, "quarter")
    print("\n### The per-seed psi6 at the end of the hold (its last psi6(t) row)\n")
    print("| cell | start | seed 0 | seed 1 | seed 2 | seed 3 | seed 4 | seed 5 |\n|---|---|---|---|---|---|---|---|")
    for c in Cb:
        for lab, _ in R.STARTS:
            print(f"| {c['name']} | {lab} | " + " | ".join(f"{v['psi_end']:.4f}" if v["ok"] else "-" for v in Vb[(c["name"], lab)]) + " |")
    print("\n(PILOT: no verdict, no paper use.)")


def selftest(root):
    """synthetic run folders: the per-run parser against direct sums; the 'stays below' rule on constructed series"""
    if os.listdir(root): sys.exit("--selftest needs an empty folder")
    ok = True
    for name, series, want in (("never", [(t, 0.85) for t in range(2001)], None),
                               ("dip shorter than 500", [(t, 0.7 if 100 <= t < 599 else 0.85) for t in range(2001)], None),
                               ("first stay at 1000", [(t, 0.7 if (300 <= t < 600 or t >= 1000) else 0.85) for t in range(2001)], 1000),
                               ("exactly 500 long", [(t, 0.7 if 1200 <= t <= 1700 else 0.85) for t in range(2001)], 1200),
                               ("at the start", [(t, 0.5) for t in range(600)], 0)):
        got = first_below([(float(t), p) for t, p in series]); ok &= got == want
        print(f"selftest stays-below, {name}: {got} (expected {want})")
    d = os.path.join(root, "run"); os.makedirs(d)
    import random
    rnd = random.Random(5); hold, nseg = 2000.0, 4
    psi = [(float(t), float(f"{0.6 + 0.1 * rnd.random():.6f}")) for t in range(int(hold) + 1)]   # as the file holds them (%.6f)
    with open(os.path.join(d, "psi6_t_L0_100_wallmassfactor_300_run0.csv"), "w") as f:
        f.write("t_sigma,phase,psi6_global,psi6_local,mean_neighbors\n")
        for t, p in psi: f.write(f"{t:.6f},hold,{p:.6f},0.7,5.5\n")
        f.write(f"{hold + 1:.6f},record,0.1,0.7,5.5\n")
    zs = []
    with open(os.path.join(d, "virial_blocks_L0_100_wallmassfactor_300_run0.csv"), "w") as f:
        f.write("t_sigma,phase,block,dt_sigma,N_0,KE_0,pair_events_0,Z_0,N_1,KE_1,pair_events_1,Z_1,pair_events,Z_global,Z_KE_weighted_mean,weighted_minus_global\n")
        for b in range(int(hold / BLOCK)):
            z = float(f"{13.5 + rnd.random():.12g}"); dt = 50.0000007 if b % 2 else 49.9999993; zs.append((z, dt))
            f.write(f"{BLOCK * (b + 1):.6f},hold,{b},{dt},200,200,1,{z},200,200,1,{z},2,{z},{z},0\n")
    open(os.path.join(d, "exit_code.txt"), "w").write("0\n")
    open(os.path.join(d, "run.log"), "w").write(
        "[EDMD3-RUN] L0=10.6 M=300 run=0 seed=1: engine=gen3 left=[orient parallel rows alternate n_s 21 n_r 10 a_s 14.25 a_r 50.17 g 2.25 s 4.5 n_vac 3 "
        "jitter 1.1 px] right=[orient parallel rows alternate n_s 21 n_r 10 a_s 14.25 a_r 50.17 g 2.25 s 4.5 n_vac 3 jitter 1.1 px] jitter_f=0.25\n"
        "[EDMD3] built #1 at t = 0\n[EDMD3-HEALTH] L0=10.6 M=300 run=0 seed=1: clean=1 ev_pair=100 ev_wall=1 ev_cross=1 ev_div=1 ev_piston=0 ev_band=0 "
        "hash=00 t_end=1 engine=gen3 build=\"x\" target=y run_s=2 engine_s=1 driver_share=0.5\n")
    v = run(d, hold, nseg); seg = hold / nseg; nbs = int(seg / BLOCK)
    for j in range(nseg):
        x = [p for t, p in psi if seg * j <= t < seg * (j + 1)]
        zz = zs[nbs * j: nbs * (j + 1)]; zt = sum(z * w for z, w in zz) / sum(w for _, w in zz)
        e = abs(v["psi"][j] - sum(x) / len(x)) < 1e-12 and abs(v["z"][j] - zt) < 1e-12; ok &= e
        print(f"selftest segment {j + 1}: psi6 {v['psi'][j]:.6f} (direct {sum(x) / len(x):.6f}), Z {v['z'][j]:.6f} (direct {zt:.6f}): {'ok' if e else '**DIFFERENT**'}")
    e = v["lattice"] and "3 vacancies" in v["lattice"] and abs(v["evs"] - 51.5) < 1e-12 and v["psi_end"] == psi[-1][1]; ok &= bool(e)
    print(f"selftest header and record: lattice '{v['lattice']}', events/s {v['evs']}, psi6 at the end {v['psi_end']}: {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


if __name__ == "__main__":
    main()
