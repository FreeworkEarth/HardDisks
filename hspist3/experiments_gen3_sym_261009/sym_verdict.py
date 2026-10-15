#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.28; plan-author decision 13, item 1): THE SYMMETRY TEST, REGISTERED -- this script, the seed
list sym_seeds.txt (its SHA-256 is SEEDS_SHA below) and the runner sym_run.py were committed together before any of the test's runs;
none of them is edited afterwards.
The plan author's text (261012 sec. 4.7.27, decision 13, item 1), verbatim:
  1. SYMMETRY TEST (registered: commit this text, the seed list and the
     verdict script before any run; script not edited afterwards)
   Cells: N = 400, H = 20, M4 seeding, divider held: eta = pi/8 and
     eta = 0.70. 40 new seeds per cell (not the six of part 3a).
   Run: 5500 sigma-time held, --gen3-virial-blocks=50; the first 500
     sigma-time discarded.
   Per seed: d = time average of Z_L - Z_R over the remaining 5000.
   Statistic per cell: t = mean(d) / (SD of d over seeds / sqrt(n)).
     The error is the seed-to-seed scatter; block SEs are not used for
     the verdict.
   Rule per cell: SYMMETRIC if |t| < 2; ASYMMETRIC if |t| >= 3;
     otherwise one extension of 40 more seeds, judged on 80 with the
     same thresholds; still between -> UNRESOLVED.
   PASS = SYMMETRIC in both cells. Inventory: all runs present and
     clean=1, one build; else INVALID.
   The script also prints: the false-alarm rates; the smallest
     asymmetry that would give t = 3, as % of Z; and as information:
     SD of d over seeds against the mean block-based SE (is the block
     SE calibrated?), sum of z_i^2 against chi-square, the same t for
     the pair-collision counts, and Z_global per cell.
   If ASYMMETRIC or INVALID: stop everything below and report.
IMPLEMENTATION (stated before any data):
  - The runs are sym_run.py's: <data>/<cell>/seed<r>/ with cell pi8 (eta = pi/8, L0 = 20) or 070 (eta = 0.70, P1's L0), r = 0..39 the
    main set and r = 40..79 the extension, the seeds of sym_seeds.txt.
  - Per run: the hold rows of virial_blocks_*_run0.csv must be blocks 0..109 (110 blocks of 50 sigma-time = 5500 sigma-time). Blocks
    0..9 (the first 500 sigma-time) are discarded; on blocks 10..109 (5000 sigma-time):
      d = sum_b (Z_0 - Z_1) dt_b / sum_b dt_b      the time average of Z_L - Z_R (Z_0 the left compartment, Z_1 the right);
      SE_block = SD_b(Z_0 - Z_1) / sqrt(100)       (information only);  z = d / SE_block;
      c = sum_b (pair_events_0 - pair_events_1) / sum_b dt_b   the left minus right pair-collision rate (information);
      Zg = sum_b Z_global dt_b / sum_b dt_b        (information).
  - Per cell: t = mean(d) / (SD(d) / sqrt(n)), SD with n - 1, n = 40 (80 after an extension).
  - RULE per cell: SYMMETRIC if |t| < 2; ASYMMETRIC if |t| >= 3; otherwise the extension (all 40 runs r = 40..79 of that cell),
    judged on all 80 with the same thresholds; still 2 <= |t| < 3 -> UNRESOLVED. PASS = SYMMETRIC in both cells.
  - INVENTORY (else INVALID, no verdict), for every run of every judged set (and of an extension set that has any folder): the folder
    exists; exit code 0; exactly one "[EDMD3] built #" line and exactly one [EDMD3-HEALTH] record, with clean=1, the seed of the
    seed list, build="5d19aa7" and target=mac-O3-e0pre (one build: BUILD below); the hold blocks are exactly 0..109; N_0 = N_1 = 200
    and KE_0, KE_1 each constant over the hold.
  - The false-alarm rates under no asymmetry: the first stage exactly (Student's t, n - 1 = 39 degrees of freedom); the whole two-stage
    rule by Monte Carlo with normal d (400000 replicates, numpy PCG64 seed 20261222); for two independent cells, P(PASS) and the
    false-alarm rate 1 - P(PASS).
  - The smallest asymmetry that gives t = 3: 3 SD(d) / sqrt(n), printed absolute and as % of the cell's mean Z_global.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_sym_261009/sym_verdict.py --data <data root>
  python3 experiments_gen3_sym_261009/sym_verdict.py --selftest <empty folder>    synthetic files: parser, inventory, statistic, rule
"""
import argparse, csv, glob, hashlib, math, os, re, sys
import numpy as np
from scipy.stats import t as student_t, chi2

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = os.path.join(HERE, "sym_seeds.txt")
SEEDS_SHA = "77d8018bd9dfdd752a8b33d2506150e9a0cd42b65e52aed1b271ebd48338c5f2"
BUILD = ('"5d19aa7"', "mac-O3-e0pre")
CELLS = (("pi8", "eta = pi/8, L0 = 20"), ("070", "eta = 0.70, L0 = 11.219974"))
NMAIN, NBLK, DISCARD, T_SYM, T_ASYM = 40, 110, 10, 2.0, 3.0


def seed_list():
    txt = open(SEEDS).read()
    rows = [l.split() for l in txt.splitlines() if l.strip()]
    return txt, {(c, int(r)): (s, int(seed)) for c, r, s, seed in rows}


def read_run(d, seed):
    """the per-run numbers and the inventory findings of one run folder (an empty findings list = clean)"""
    bad = []
    if not os.path.isdir(d): return None, ["missing"]
    ec = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
    if ec != "0": bad.append(f"exit {ec}")
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M)); recs = re.findall(r"^\[EDMD3-HEALTH\] (.*)$", lg, re.M)
    if nb != 1: bad.append(f"{nb} build lines")
    if len(recs) != 1: bad.append(f"{len(recs)} run records")
    else:
        kv = dict(re.findall(r" (\w+)=(\S+)", " " + recs[0].split(": ", 1)[-1]))
        sd = re.search(r" seed=(\d+):", " " + recs[0])
        if kv.get("clean") != "1": bad.append(f"clean={kv.get('clean')}")
        if not sd or int(sd.group(1)) != seed: bad.append(f"seed {sd.group(1) if sd else None} != {seed}")
        if (kv.get("build"), kv.get("target")) != BUILD: bad.append(f"build {kv.get('build')} {kv.get('target')}")
    p = glob.glob(os.path.join(d, "virial_blocks_*_run0.csv"))
    if len(p) != 1: return None, bad + [f"{len(p)} virial-block files"]
    rows = [r for r in csv.DictReader(open(p[0])) if r["phase"] == "hold"]
    if [int(r["block"]) for r in rows] != list(range(NBLK)): bad.append("hold blocks are not 0..109")
    if any(r["N_0"] != "200" or r["N_1"] != "200" for r in rows): bad.append("N per compartment")
    if len({r["KE_0"] for r in rows}) != 1 or len({r["KE_1"] for r in rows}) != 1: bad.append("KE not constant")
    if bad: return None, bad
    use = rows[DISCARD:NBLK]
    dt = np.array([float(r["dt_sigma"]) for r in use]); dz = np.array([float(r["Z_0"]) - float(r["Z_1"]) for r in use])
    dn = np.array([int(r["pair_events_0"]) - int(r["pair_events_1"]) for r in use], float)
    zg = np.array([float(r["Z_global"]) for r in use])
    d_ = float((dz * dt).sum() / dt.sum()); se = float(dz.std(ddof=1) / math.sqrt(len(dz)))
    return dict(d=d_, se=se, z=d_ / se, c=float(dn.sum() / dt.sum()), zg=float((zg * dt).sum() / dt.sum())), []


def tstat(v):
    v = np.asarray(v, float); return float(v.mean() / (v.std(ddof=1) / math.sqrt(len(v))))


def outcome(t):
    return "SYMMETRIC" if abs(t) < T_SYM else ("ASYMMETRIC" if abs(t) >= T_ASYM else "BETWEEN")


def judge(dm, de):
    """the rule for one cell: (outcome, t on the main 40, t on all 80 or None)"""
    t1 = tstat(dm); o = outcome(t1)
    if o != "BETWEEN": return o, t1, None
    if de is None or len(de) != NMAIN: return "EXTENSION REQUIRED (not run yet)", t1, None
    t2 = tstat(list(dm) + list(de)); o2 = outcome(t2)
    return ("UNRESOLVED" if o2 == "BETWEEN" else o2), t1, t2


def rates(n=NMAIN, reps=400000, chunk=50000):
    """under no asymmetry: exact first stage (Student's t, n - 1 dof) and the two-stage rule by Monte Carlo (normal d)"""
    p3 = 2 * student_t.sf(T_ASYM, n - 1); p2 = 2 * student_t.sf(T_SYM, n - 1)
    g = np.random.Generator(np.random.PCG64(20261222)); cnt = dict(SYMMETRIC=0, ASYMMETRIC=0, UNRESOLVED=0, ext=0); done = 0
    while done < reps:
        x = g.standard_normal((chunk, 2 * n))
        def tt(a): return a.mean(1) / (a.std(1, ddof=1) / math.sqrt(a.shape[1]))
        t1 = np.abs(tt(x[:, :n])); t2 = np.abs(tt(x))
        mid = (t1 >= T_SYM) & (t1 < T_ASYM)
        cnt["ext"] += int(mid.sum())
        cnt["SYMMETRIC"] += int((t1 < T_SYM).sum() + (mid & (t2 < T_SYM)).sum())
        cnt["ASYMMETRIC"] += int((t1 >= T_ASYM).sum() + (mid & (t2 >= T_ASYM)).sum())
        cnt["UNRESOLVED"] += int((mid & (t2 >= T_SYM) & (t2 < T_ASYM)).sum())
        done += chunk
    return dict(p3=p3, pmid=p2 - p3, p2=p2, **{k: v / done for k, v in cnt.items()}, reps=done)


def print_rates():
    R = rates(); se = lambda p: math.sqrt(p * (1 - p) / R["reps"])
    print(f"| quantity (no asymmetry) | one cell | two cells (independent) |\n|---|---|---|")
    print(f"| first stage, exact (t, {NMAIN - 1} dof): P(abs(t) >= {T_ASYM:g}) | {R['p3']:.5f} | {1 - (1 - R['p3']) ** 2:.5f} (at least one) |")
    print(f"| first stage, exact: P({T_SYM:g} <= abs(t) < {T_ASYM:g}) = P(extension) | {R['pmid']:.5f} | {1 - (1 - R['pmid']) ** 2:.5f} (at least one) |")
    print(f"| whole rule (Monte Carlo, {R['reps']} replicates): P(extension) | {R['ext']:.5f} ({se(R['ext']):.1e}) | |")
    print(f"| whole rule: P(ASYMMETRIC) | {R['ASYMMETRIC']:.5f} ({se(R['ASYMMETRIC']):.1e}) | {1 - (1 - R['ASYMMETRIC']) ** 2:.5f} (at least one) |")
    print(f"| whole rule: P(UNRESOLVED) | {R['UNRESOLVED']:.5f} ({se(R['UNRESOLVED']):.1e}) | {1 - (1 - R['UNRESOLVED']) ** 2:.5f} (at least one) |")
    print(f"| whole rule: P(SYMMETRIC) | {R['SYMMETRIC']:.5f} ({se(R['SYMMETRIC']):.1e}) | P(PASS) = {R['SYMMETRIC'] ** 2:.5f}; false alarm (not PASS) = {1 - R['SYMMETRIC'] ** 2:.5f} |")
    print("(Monte Carlo standard errors in brackets)")


def analyse(data, check_seeds=True):
    txt, S = seed_list()
    if check_seeds and hashlib.sha256(txt.encode()).hexdigest() != SEEDS_SHA:
        print("**the seed list is not the registered one (SHA-256)**"); return "INVALID"
    print("## Inventory\n")
    print("| cell | set | runs present | clean (exit 0, one build line, one record clean=1, seed, build) | blocks and compartments | findings |\n|---|---|---|---|---|---|")
    ok, res = True, {}
    for cell, _ in CELLS:
        for st in ("main", "ext"):
            rs = [r for (c, r), (s, _) in S.items() if c == cell and s == st]
            dirs = [os.path.join(data, cell, f"seed{r}") for r in rs]
            if st == "ext" and not any(os.path.isdir(d) for d in dirs): continue
            got, bad = {}, []
            for r, d in zip(rs, dirs):
                v, b = read_run(d, S[(cell, r)][1])
                if b: bad.append(f"seed{r}: {', '.join(b)}")
                else: got[r] = v
            present = sum(os.path.isdir(d) for d in dirs)
            ok &= not bad and len(got) == len(rs)
            res[(cell, st)] = got
            print(f"| {cell} | {st} | {present} of {len(rs)} | {len(got)} of {len(rs)} | {'yes' if not bad else '**no**'} | {'; '.join(bad[:6]) + (' ...' if len(bad) > 6 else '') if bad else '-'} |")
    print(f"\nbuild required: build={BUILD[0]} target={BUILD[1]}; inventory {'clean' if ok else '**NOT CLEAN -> INVALID**'}")
    if not ok: print("\nSYMMETRY TEST: INVALID -- the inventory is not clean; no verdict"); return "INVALID"
    print("\n## The registered statistic per cell (blocks 10..109 of the hold, 5000 sigma-time; d = time average of Z_L - Z_R)\n")
    print("| cell | judged on | mean d | SD of d over seeds | SE = SD / sqrt(n) | t | outcome by the rule |\n|---|---|---|---|---|---|---|")
    outs, judged = {}, {}
    for cell, lab in CELLS:
        m = res[(cell, "main")]; dm = [m[r]["d"] for r in sorted(m)]; e = res.get((cell, "ext"))
        de = [e[r]["d"] for r in sorted(e)] if e else None
        o, t1, t2 = judge(dm, de)
        print(f"| {cell} ({lab}) | {len(dm)} | {np.mean(dm):+.6f} | {np.std(dm, ddof=1):.6f} | {np.std(dm, ddof=1) / math.sqrt(len(dm)):.6f} | {t1:+.2f} | "
              f"**{o if t2 is None else 'EXTENSION'}** |")
        use = m
        if t2 is not None:
            use = {**m, **e}; da = dm + de
            print(f"| {cell} ({lab}) | {len(da)} | {np.mean(da):+.6f} | {np.std(da, ddof=1):.6f} | {np.std(da, ddof=1) / math.sqrt(len(da)):.6f} | {t2:+.2f} | **{o}** |")
        outs[cell] = o; judged[cell] = use
    print("\n## The thresholds' false-alarm rates (no asymmetry)\n"); print_rates()
    print("\n## The smallest asymmetry that would give t = 3 (with the observed SD of d)\n")
    print("| cell | judged on | 3 SD(d) / sqrt(n) | mean Z_global | as % of Z_global |\n|---|---|---|---|---|")
    for cell, lab in CELLS:
        u = judged[cell]; dd = np.array([u[r]["d"] for r in u]); zg = np.mean([u[r]["zg"] for r in u])
        a3 = 3 * dd.std(ddof=1) / math.sqrt(len(dd))
        print(f"| {cell} | {len(dd)} | {a3:.6f} | {zg:.5f} | {100 * a3 / zg:.4f} |")
    print("\n## Information (no verdict)\n")
    print("| cell | n | SD of d over seeds | mean block-based SE | ratio SD / mean SE | sum of z_i^2 | chi-square, n dof: mean, P(>= sum) | "
          "t of the pair-collision rate difference (mean L - R per sigma-time) | Z_global: mean (SE over seeds) |\n|---|---|---|---|---|---|---|---|---|")
    for cell, lab in CELLS:
        u = judged[cell]; rs = sorted(u); n = len(rs)
        dd = np.array([u[r]["d"] for r in rs]); se = np.array([u[r]["se"] for r in rs]); zz = np.array([u[r]["z"] for r in rs])
        cc = np.array([u[r]["c"] for r in rs]); zg = np.array([u[r]["zg"] for r in rs]); s2 = float((zz ** 2).sum())
        print(f"| {cell} | {n} | {dd.std(ddof=1):.6f} | {se.mean():.6f} | {dd.std(ddof=1) / se.mean():.3f} | {s2:.1f} | {n}, {chi2.sf(s2, n):.3g} | "
              f"{tstat(cc):+.2f} ({cc.mean():+.3f}) | {zg.mean():.5f} ({zg.std(ddof=1) / math.sqrt(n):.5f}) |")
    print("\n## Per seed (information)\n")
    print("| cell | r | seed | d | block SE | z | pair-rate L - R [1/sigma-time] | Z_global |\n|---|---|---|---|---|---|---|---|")
    _, S2 = seed_list()
    for cell, _ in CELLS:
        u = judged[cell]
        for r in sorted(u):
            v = u[r]
            print(f"| {cell} | {r} | {S2[(cell, r)][1]} | {v['d']:+.6f} | {v['se']:.6f} | {v['z']:+.2f} | {v['c']:+.3f} | {v['zg']:.5f} |")
    if any(o == "ASYMMETRIC" for o in outs.values()): verdict = "FAIL (ASYMMETRIC)"   # decides even if the other cell needs an extension
    elif any(o.startswith("EXTENSION REQUIRED") for o in outs.values()):
        verdict = "EXTENSION REQUIRED at " + ", ".join(c for c, o in outs.items() if o.startswith("EXTENSION REQUIRED"))
    elif all(o == "SYMMETRIC" for o in outs.values()): verdict = "PASS"
    else: verdict = "FAIL (UNRESOLVED)"
    ext = [c for c in outs if (c, "ext") in res]
    print("\n" + "; ".join(f"{c}: {o}" for c, o in outs.items()) + (f"; the extension was run for {', '.join(ext)} (disclosed)" if ext else "; no extension"))
    print(f"SYMMETRY TEST: {verdict}")
    return verdict


def selftest(root):
    """synthetic run folders in an empty folder: the parser and the statistic against an independent computation; the rule's paths"""
    if os.listdir(root): sys.exit("--selftest needs an empty folder")
    txt, S = seed_list(); g = np.random.Generator(np.random.PCG64(7)); hdr = ("t_sigma,phase,block,dt_sigma,N_0,KE_0,pair_events_0,Z_0,N_1,KE_1,"
                                                                              "pair_events_1,Z_1,pair_events,Z_global,Z_KE_weighted_mean,weighted_minus_global\n")
    def make(base, shift, sd=0.003, bad=None, ext=False):
        want = {}
        for (cell, r), (st, seed) in S.items():
            if st == "ext" and not ext: continue
            d = os.path.join(base, cell, f"seed{r}"); os.makedirs(d)
            open(os.path.join(d, "exit_code.txt"), "w").write("0\n")
            clean = "0" if bad == (cell, r) else "1"
            open(os.path.join(d, "run.log"), "w").write(f"[EDMD3] built #1 at t = 0\n[EDMD3-HEALTH] L0=20.0 M=300 run=0 seed={seed}: clean={clean} "
                                                         f"hash=0123 engine=gen3 build=\"5d19aa7\" target=mac-O3-e0pre run_s=1 engine_s=1 driver_share=0\n")
            mu = shift[cell] + g.normal(0, sd); x = []
            with open(os.path.join(d, "virial_blocks_L0_200_wallmassfactor_300_run0.csv"), "w") as f:
                f.write(hdr)
                for b in range(NBLK):
                    dt = 50.0000007
                    z0 = float(f"{2.82 + mu / 2 + g.normal(0, 0.01):.12g}"); z1 = float(f"{2.82 - mu / 2 + g.normal(0, 0.01):.12g}")
                    if b >= DISCARD: x.append((z0 - z1, dt))
                    f.write(f"{50.000001 * (b + 1):.6f},hold,{b},{dt},200,200,20500,{z0:.12g},200,200,20400,{z1:.12g},40900,{(z0 + z1) / 2:.12g},{(z0 + z1) / 2:.12g},0\n")
                f.write(f"5501,record,{NBLK},1,200,200,1,2.8,200,200,1,2.8,2,2.8,2.8,0\n")
            want[(cell, r)] = sum(a * w for a, w in x) / sum(w for _, w in x)
        return want
    ok = True
    base_ = np.array([1.0, -1.0] * (NMAIN // 2)); sdm = base_.std(ddof=1) / math.sqrt(NMAIN)
    for lab, t_main, ext_mu, expect in (("t 1.5", 1.5, None, "SYMMETRIC"), ("t -3.5", -3.5, None, "ASYMMETRIC"),
                                         ("t 2.5, no extension", 2.5, None, "EXTENSION REQUIRED (not run yet)"),
                                         ("t 2.5, extension at the same mean", 2.5, 1.0, "ASYMMETRIC"),
                                         ("t 2.5, extension at mean 0", 2.5, 0.0, "SYMMETRIC"),
                                         ("t 2.5, extension at 0.6 of the mean", 2.5, 0.6, "UNRESOLVED")):
        dm = base_ + t_main * sdm; de = None if ext_mu is None else base_ + ext_mu * t_main * sdm
        o, t1, t2 = judge(dm, de); ok &= o == expect
        print(f"selftest rule, {lab}: t1 {t1:+.2f}, t2 {'-' if t2 is None else f'{t2:+.2f}'} -> {o} (expected {expect})")
    for name, shift, bad, ext, expect in (("null", {"pi8": 0.0, "070": 0.0}, None, False, None),
                                          ("asymmetric", {"pi8": 0.01, "070": 0.0}, None, False, "FAIL (ASYMMETRIC)"),
                                          ("one run clean=0", {"pi8": 0.0, "070": 0.0}, ("070", 7), False, "INVALID")):
        base = os.path.join(root, name.replace(" ", "_").replace("=", "")); os.makedirs(base)
        want = make(base, shift, bad=bad, ext=ext)
        for cell, _ in CELLS:
            if bad: break
            got = [read_run(os.path.join(base, cell, f"seed{r}"), S[(cell, r)][1])[0]["d"] for r in range(NMAIN)]
            ref = [want[(cell, r)] for r in range(NMAIN)]
            ok &= max(abs(a - b) for a, b in zip(got, ref)) < 1e-12
            print(f"selftest {name} {cell}: d parsed = d computed independently (max diff {max(abs(a - b) for a, b in zip(got, ref)):.1e}); "
                  f"t = {tstat(got):+.2f} (independent: {np.mean(ref) / (np.std(ref, ddof=1) / math.sqrt(NMAIN)):+.2f})")
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf): v = analyse(base, check_seeds=False)
        print(f"selftest {name}: verdict {v}" + (f" (expected {expect})" if expect else ""))
        if expect: ok &= v == expect
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data"); ap.add_argument("--selftest"); a = ap.parse_args()
    if a.selftest: return selftest(os.path.abspath(a.selftest))
    if SEEDS_SHA.startswith("<"): sys.exit("STOP: not registered")
    print("# The symmetry test (261012 sec. 4.7.28; plan-author decision 13, item 1), printed by experiments_gen3_sym_261009/sym_verdict.py\n")
    print(f"RULE per cell: SYMMETRIC if abs(t) < {T_SYM:g}; ASYMMETRIC if abs(t) >= {T_ASYM:g}; otherwise one extension (40 more seeds), judged on 80 "
          f"with the same thresholds, still between -> UNRESOLVED. PASS = SYMMETRIC in both cells. t = mean(d) / (SD(d) / sqrt(n)), the seed-to-seed scatter.\n")
    analyse(os.path.abspath(a.data))


if __name__ == "__main__":
    main()
