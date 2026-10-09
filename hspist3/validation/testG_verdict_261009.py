#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.15, stage B of the plan-author programme of sec. 4.7.12): TEST G, the Mac part of the gen-3
gate -- written and committed before any of its trajectories; not edited after its registration commit.
The logic of validation/resched_testTprime_261007.py with "policy" replaced by "engine" (gen3 in minimal's place, gen2 in
legacy's), and the rule of the programme:
Cell epi8_H_H10_L10 (eta = pi/8, H = L0 = 10, N = 100), the speed-of-sound protocol and command of Test T-prime except the seeds and
--engine; engines gen2 (default policy) and gen3, one binary, one build; M = 300 and M = 1500, 400 trajectories per engine per mass,
the same seeds for both engines, interleaved (cluster/gen3_gate_261009/gen_testG_261009.py, run_testG_mac.py, testG_worker.sh).
Registered estimator: the mean nu per mass (red_nu.csv, the canonical argmax estimator of reduce_B.py, unchanged), SE = SD/sqrt(n)
per engine, z = (mean_gen3 - mean_gen2) / sqrt(SE_gen3^2 + SE_gen2^2); 95 % interval of the relative difference
100 (d -+ 1.96 s) / mean_gen2 (T-prime's formula), the stated bound.
INVENTORY (else INVALID, no verdict): all 1600 trajectories present (and the extension's, if run); every gen3 run record clean=1
(exactly one record and one build line per trajectory, engine=gen3); no gen2 health line (and no gen3 line in a gen2 trajectory);
no failed trajectory; one build (every cell's .build_git = the registered binary's version line).
RULE per mass (the programme's, verbatim in 261012 sec. 4.7.12): NO DIFFERENCE if |z| < 2; DIFFERENCE if |z| >= 3; otherwise ONE
extension: 400 more trajectories per engine at that mass (r = 400..799 of the same seed stream), judged on all 800 with the same two
thresholds; still 2 <= |z| < 3 -> UNRESOLVED. TEST G (Mac) PASS = NO DIFFERENCE at both masses. The extension is disclosed.
Also printed: the false-alarm rates of the thresholds for one and for two masses (with the extension path), and the smallest true
difference that gives z = 3 with the observed SEs.
INFORMATION rows (no verdict): nu_d (sec. 4.4.13 item 4a) on the main masses; M = 50 and M = 2000 (100 per engine); the dense cell
(the gate's dense geometry, eta = 0.70073, M = 300, 200 per engine); N = 400 at pi/8 (H = L0 = 20, M = 300, 100 per engine); the
static method: the A-fixed cell at x_0, 200 runs per engine, the mean force per face F_L, F_R of the registered Method A
estimator (reduce_AF.py).
usage (from hspist3/):
  python3 validation/testG_verdict_261009.py --design        thresholds, false-alarm rates, expected precision (T-prime's SDs)
  python3 validation/testG_verdict_261009.py --dry-run-Tp    the rule on T-prime's data (gen3 <- minimal, gen2 <- legacy): must
                                                              reproduce T-prime's z = +0.49 (M = 300) and -0.42 (M = 1500)
  python3 validation/testG_verdict_261009.py                 the analysis of the Test G data (experiments_gen3_gate_261009/)
"""
import glob, hashlib, math, os, re, sys
from collections import Counter
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GG = os.path.join(HS, "cluster", "gen3_gate_261009"); GATE = os.path.join(HS, "cluster", "resched_gate_261005")
sys.path.insert(0, HERE); sys.path.insert(0, HS); sys.path.insert(0, GG); sys.path.insert(0, GATE)
import resched_testT_261007 as TT
import gen_testG_261009 as G

CID, MASSES, NSEED = G.CID, G.MASSES, G.NSEED
LOC = os.environ.get("HD_TESTG_LOC", G.LOC)                     # env: test hook only
TASKS = os.path.join(GG, "tasks_testG_261009.txt")
ENG = ("gen3", "gen2")                                           # gen3 in minimal's place (the new one first, as T-prime)
Z_ND, Z_D, Z95 = 2.0, 3.0, 1.96
BUILD = "00ALLINONE  git f42befb  target mac-O3-e0pre"           # the registered binary (sec. 4.7.15)
DESIGN, DRY = "--design" in sys.argv, "--dry-run-Tp" in sys.argv
HEALTH = re.compile(r"EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue")


def cell_dir(eng, cid, M):
    return os.path.join(LOC, G.REL, eng, cid, f"m_{M}")


def load_cell(eng, cid, M):
    d = cell_dir(eng, cid, M)
    p = os.path.join(d, "red_nu.csv")
    return dict(dir=d, red=(pd.read_csv(p).sort_values("run").reset_index(drop=True) if os.path.exists(p) else None),
                log=(open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""),
                failed=[open(f, errors="ignore").read() for f in glob.glob(os.path.join(d, ".failed_run*", "stdout.log"))],
                build=(open(os.path.join(d, ".build_git")).read().strip() if os.path.exists(os.path.join(d, ".build_git")) else None))


def load_Tp():
    """T-prime's data, gen3 <- minimal, gen2 <- legacy (the dry run)."""
    import resched_testTprime_261007 as TP
    D = TP.load("Tprime")
    return {("gen3", M): D[("minimal", M)] for M in MASSES} | {("gen2", M): D[("legacy", M)] for M in MASSES}


def stats(a, b):
    sa, sb = a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b))
    d = a.mean() - b.mean(); s = math.hypot(sa, sb)
    return dict(na=len(a), nb=len(b), ma=a.mean(), sa=sa, mb=b.mean(), sb=sb, d=d, s=s, rel=100 * d / b.mean(),
                lo=100 * (d - Z95 * s) / b.mean(), hi=100 * (d + Z95 * s) / b.mean(), z=d / s)


def outcome(z):
    if not math.isfinite(z): return "INVALID (no z)"
    return "NO DIFFERENCE" if abs(z) < Z_ND else ("DIFFERENCE" if abs(z) >= Z_D else "EXTENSION")


def compare(D):
    """Per mass: the 400-trajectory numbers; if they ask for the extension and runs 400..799 exist for both engines, the 800."""
    rows = {}
    for M in MASSES:
        A, B = D[("gen3", M)]["red"], D[("gen2", M)]["red"]
        a4, b4 = A[A["run"] < NSEED]["nu"].to_numpy(float), B[B["run"] < NSEED]["nu"].to_numpy(float)
        r = stats(a4, b4); r["out"] = outcome(r["z"]); r["judged_on"] = NSEED; r["ext"] = None
        if r["out"] == "EXTENSION":
            a8, b8 = A["nu"].to_numpy(float), B["nu"].to_numpy(float)
            if len(a8) == len(b8) == 2 * NSEED:
                e = stats(a8, b8); o = outcome(e["z"]); e["out"] = "UNRESOLVED" if o == "EXTENSION" else o
                e["judged_on"] = 2 * NSEED; r["ext"] = e
            else:
                r["out"] = "EXTENSION REQUIRED (not run yet)"
        rows[M] = r
    return rows


def final(r):
    return r["ext"]["out"] if r["ext"] else r["out"]


def print_rule(rows):
    print("| mass | judged on (per engine) | n gen3 | n gen2 | gen3 mean nu | SE | gen2 mean nu | SE | difference | relative [%] | "
          "95 % interval [%] | z | outcome by the rule |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for M, r0 in rows.items():
        for r in [r0] + ([r0["ext"]] if r0["ext"] else []):
            print(f"| {M} | {r['judged_on']} | {r['na']} | {r['nb']} | {r['ma']:.7f} | {r['sa']:.2e} | {r['mb']:.7f} | {r['sb']:.2e} | {r['d']:+.3e} | "
                  f"{r['rel']:+.3f} | [{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | **{r['out']}** |")


def rates():
    """Under the null at one mass: P(DIFFERENCE), P(extension), P(UNRESOLVED), P(NO DIFFERENCE); z800 | z400 ~ N(z400/sqrt2, 1/2)."""
    def cond(z, lo, hi):   # P(lo <= |z800| < hi | z400 = z)
        m, s = z / math.sqrt(2), math.sqrt(0.5)
        f = lambda x: norm.cdf((x - m) / s)
        return (f(hi) - f(lo)) + (f(-lo) - f(-hi))
    ext = 2 * (norm.cdf(Z_D) - norm.cdf(Z_ND))
    pd_ext = 2 * quad(lambda z: norm.pdf(z) * cond(z, Z_D, np.inf), Z_ND, Z_D)[0]
    pu = 2 * quad(lambda z: norm.pdf(z) * cond(z, Z_ND, Z_D), Z_ND, Z_D)[0]
    pn_ext = 2 * quad(lambda z: norm.pdf(z) * cond(z, 0.0, Z_ND), Z_ND, Z_D)[0]
    p_d = 2 * norm.sf(Z_D) + pd_ext; p_n = (1 - 2 * norm.sf(Z_ND)) + pn_ext
    return dict(diff=p_d, ext=ext, unres=pu, nd=p_n, d1=2 * norm.sf(Z_D), n1=1 - 2 * norm.sf(Z_ND))


def print_rates():
    R = rates()
    print("| quantity (no true difference) | one mass | two masses (independent) |\n|---|---|---|")
    print(f"| P(abs(z) >= {Z_D:g} at 400) | {R['d1']:.5f} | {1 - (1 - R['d1']) ** 2:.5f} |")
    print(f"| P(extension: {Z_ND:g} <= abs(z) < {Z_D:g} at 400) | {R['ext']:.5f} | {1 - (1 - R['ext']) ** 2:.5f} (at least one) |")
    print(f"| P(DIFFERENCE), extension path included | {R['diff']:.5f} | {1 - (1 - R['diff']) ** 2:.5f} (at least one) |")
    print(f"| P(UNRESOLVED) | {R['unres']:.5f} | {1 - (1 - R['unres']) ** 2:.5f} (at least one) |")
    print(f"| P(NO DIFFERENCE) | {R['nd']:.5f} | P(PASS) = {R['nd'] ** 2:.5f}; false alarm (not PASS) = {1 - R['nd'] ** 2:.5f} |")


def design():
    print("# Test G -- design numbers (261012 sec. 4.7.15), before any of its data\n")
    print("### The thresholds' false-alarm rates under the null\n"); print_rates()
    D = load_Tp(); sd = {(e, M): D[(e, M)]["red"]["nu"].std(ddof=1) for e in ENG for M in MASSES}
    print("\n### Expected precision with 400 + 400 per mass (per-seed SDs of T-prime, gen3 <- minimal, gen2 <- legacy; an assumption "
          "for gen3)\n")
    print("| M | SD (T-prime minimal) | SD (T-prime legacy) | mean nu (T-prime legacy) | SE of the difference | smallest true difference "
          "with z = 3 (expected) [%] |\n|---|---|---|---|---|---|")
    for M in MASSES:
        s = math.sqrt(sd[("gen3", M)] ** 2 + sd[("gen2", M)] ** 2) / math.sqrt(NSEED); mb = D[("gen2", M)]["red"]["nu"].mean()
        print(f"| {M} | {sd[('gen3', M)]:.3e} | {sd[('gen2', M)]:.3e} | {mb:.7f} | {s:.2e} | {100 * 3 * s / mb:.3f} |")
    print("\n### Seeds and task list\n")
    S = G.seed_text(); print(f"seed list SHA-256 {hashlib.sha256(S.encode()).hexdigest()}; task list {os.path.relpath(TASKS, HS)} "
                             f"SHA-256 {hashlib.sha256(open(TASKS, 'rb').read()).hexdigest()} ({sum(1 for _ in open(TASKS))} lines)")


def sections(log):
    return log.split("##RUN")[1:]


def inventory(D):
    want = {}
    lists = [TASKS] + sorted(glob.glob(os.path.join(GG, "tasks_testG_ext_M*_261009.txt")))
    for p in lists:
        for l in open(p):
            f = l.split()
            if f[0] == "B" and f[-1] == "main": want.setdefault((f[10], int(f[2])), set()).add((int(f[3]), int(f[4])))
    print("| engine | M | trajectories matched (expected) | n finite > 0 | log sections | health findings (failed runs incl.) | engine line wrong | "
          "max contact gap [px] | contact lines missing | nodes | ok |\n|---|---|---|---|---|---|---|---|---|---|---|")
    ok_all = True; builds = set(); nodes = {e: Counter() for e in ENG}
    for eng in ENG:
        for M in MASSES:
            e = D[(eng, M)]; red = e["red"]
            if red is None:
                ok_all = False; print(f"| {eng} | {M} | **no red_nu.csv** | | | | | | | | **NO** |"); continue
            have = {(int(a), int(b)) for a, b in zip(red["run"], red["seed"])}
            exp = want.get((eng, M), set()); matched = len(exp & have)
            nfin = int((np.isfinite(red["n"].astype(float)) & (red["n"].astype(float) > 0)).sum())
            secs = sections(e["log"])
            if eng == "gen2":
                nh = sum(len(HEALTH.findall(s)) for s in secs) + sum(len(HEALTH.findall(f)) for f in e["failed"]) + len(e["failed"])
                bad = sum(("[EDMD-RESCHED] divider events: minimal (default)" not in s) or ("[EDMD3" in s) for s in secs)
            else:
                nh = len(e["failed"])
                for s in secs:
                    recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", s, re.M)
                    other = [l for l in s.splitlines() if HEALTH.search(l) and not l.startswith("[EDMD3-HEALTH]")]
                    nh += int(recs != ["1"]) + len(other)
                bad = sum(len(re.findall(r"^\[EDMD3\] built #", s, re.M)) != 1 or " engine=gen3 " not in s for s in secs)
            cg, nmiss = 0.0, 0
            for s in secs:
                m = TT.CONTACT.search(s)
                if not m: nmiss += 1; continue
                cg = max(cg, max(float(v) for v in m.groups()[1:]))
                n = TT.NODE.search(s.split("\n", 1)[0])
                if n: nodes[eng][n.group(3)] += 1
            if e["build"]: builds.add(e["build"])
            ok = matched == len(exp) and len(exp) in (NSEED, 2 * NSEED) and len(red) == len(exp) and nfin == len(red) and len(secs) == len(red) \
                and nh == 0 and bad == 0 and nmiss == 0 and cg <= 1e-6
            ok_all &= ok
            print(f"| {eng} | {M} | {matched} ({len(exp)}) | {nfin} | {len(secs)} | {nh} | {bad} | {cg:.2e} | {nmiss} | {len(nodes[eng])} | {'yes' if ok else '**NO**'} |")
    ok_all &= builds == {BUILD}
    print(f"\nnodes: gen3 {dict(nodes['gen3'])}; gen2 {dict(nodes['gen2'])}")
    print(f"builds (.build_git): {sorted(builds)} (required: {BUILD}); inventory {'clean' if ok_all else '**NOT CLEAN -> INVALID**'}")
    return ok_all


def information_nud(D):
    import warnings; warnings.filterwarnings("ignore")
    import resched_testT_followup_261007 as F
    print("| M | estimator | gen3 mean | gen2 mean | relative [%] | z | same-seed correlation |\n|---|---|---|---|---|---|---|")
    for M in MASSES:
        ref = 0.5 * (D[("gen3", M)]["red"]["nu"].mean() + D[("gen2", M)]["red"]["nu"].mean()); v = {}
        for eng in ENG:
            dd = D[(eng, M)]["dir"]; red = D[(eng, M)]["red"]; z = np.load(os.path.join(dd, "acf_runs.npz")); dt = float(red["dt"].iloc[0]); w = []
            for r in red["run"]:
                try: w.append(F.fit_acf(z[f"run{r}"].astype(float), dt, ref)[4] / (2 * math.pi))
                except Exception: w.append(np.nan)
            v[eng] = np.array(w)
        for lab, a, b in (("argmax (registered)", D[("gen3", M)]["red"]["nu"].to_numpy(float), D[("gen2", M)]["red"]["nu"].to_numpy(float)),
                          ("nu_d (sec. 4.4.13 item 4a)", v["gen3"], v["gen2"])):
            ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
            s = math.hypot(a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b)))
            print(f"| {M} | {lab} | {a.mean():.7f} | {b.mean():.7f} | {100 * (a.mean() - b.mean()) / b.mean():+.3f} | {(a.mean() - b.mean()) / s:+.2f} | "
                  f"{np.corrcoef(a, b)[0, 1]:+.3f} |")


def information_rows():
    print("| row | n gen3 | n gen2 | gen3 mean | SE | gen2 mean | SE | relative [%] | 95 % interval [%] | z | clean (both engines) |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|")
    for lab, cid, M in (("nu, M = 50", CID, 50), ("nu, M = 2000", CID, 2000), ("nu, dense cell, M = 300", "dense_H10_L5.604167", 300),
                        ("nu, N = 400 at pi/8, M = 300", "epi8_N400_H20_L20", 300)):
        c = {e: load_cell(e, cid, M) for e in ENG}
        if any(c[e]["red"] is None for e in ENG): print(f"| {lab} | missing | | | | | | | | | |"); continue
        a, b = c["gen3"]["red"]["nu"].to_numpy(float), c["gen2"]["red"]["nu"].to_numpy(float); r = stats(a, b)
        fl = sum(len(c[e]["failed"]) for e in ENG)
        recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", c["gen3"]["log"], re.M)
        g2h = sum(len(HEALTH.findall(s)) for s in sections(c["gen2"]["log"]))
        cl = f"failed {fl}; gen3 records {len(recs)}, clean=1 {recs.count('1')}; gen2 health lines {g2h}"
        print(f"| {lab} | {r['na']} | {r['nb']} | {r['ma']:.7f} | {r['sa']:.2e} | {r['mb']:.7f} | {r['sb']:.2e} | {r['rel']:+.3f} | "
              f"[{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | {cl} |")
    for face in ("F_L", "F_R"):
        v = {}
        for e in ENG:
            d = os.path.join(LOC, G.REL, e, CID, "x_0")
            fs = sorted(glob.glob(os.path.join(d, "red_*.csv")))
            v[e] = np.array([float(pd.read_csv(f)[face].iloc[0]) for f in fs]) if fs else np.array([])
        if min(len(v[e]) for e in ENG) < 2: print(f"| A-fixed {face} | missing | | | | | | | | | |"); continue
        r = stats(v["gen3"], v["gen2"])
        nf = {e: len(glob.glob(os.path.join(LOC, G.REL, e, CID, "x_0", ".failed_*"))) for e in ENG}
        print(f"| static method, A-fixed at x_0, mean {face} | {r['na']} | {r['nb']} | {r['ma']:.6f} | {r['sa']:.2e} | {r['mb']:.6f} | {r['sb']:.2e} | "
              f"{r['rel']:+.3f} | [{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | failed gen3 {nf['gen3']}, gen2 {nf['gen2']} |")


def main():
    if DESIGN: return design()
    print(f"# Test G, the Mac part of the gen-3 gate (261012 sec. 4.7.15){'  [DRY RUN on T-prime data: a test of this script, no verdict]' if DRY else ''}\n")
    print(f"RULE per mass: NO DIFFERENCE if abs(z) < {Z_ND:g}; DIFFERENCE if abs(z) >= {Z_D:g}; otherwise one extension (400 more per engine), "
          f"judged on 800 with the same thresholds, still between -> UNRESOLVED. PASS = NO DIFFERENCE at both masses.\n")
    if DRY:
        D = load_Tp()
    else:
        D = {(e, M): load_cell(e, CID, M) for e in ENG for M in MASSES}
        print("## Inventory\n"); inv = inventory(D)
        if not inv:
            print("\nTEST G (Mac): INVALID -- the inventory is not clean; no verdict"); return
    rows = compare(D)
    print("\n## The registered numbers\n"); print_rule(rows)
    if DRY:
        ok = round(rows[300]["z"], 2) == 0.49 and round(rows[1500]["z"], 2) == -0.42
        print(f"\ndry run: z(M = 300) = {rows[300]['z']:+.2f}, z(M = 1500) = {rows[1500]['z']:+.2f}; T-prime printed +0.49 and -0.42: "
              f"{'REPRODUCED' if ok else '**NOT REPRODUCED**'}")
    print("\n## The thresholds' false-alarm rates (no true difference)\n"); print_rates()
    print("\n## The smallest true difference that would give z = 3 with the observed SEs\n")
    print("| mass | judged on | SE of the difference | 3 x SE | relative to gen2's mean [%] |\n|---|---|---|---|---|")
    for M, r0 in rows.items():
        r = r0["ext"] or r0
        print(f"| {M} | {r['judged_on']} | {r['s']:.2e} | {3 * r['s']:.2e} | {100 * 3 * r['s'] / r['mb']:.3f} |")
    outs = {M: final(r) for M, r in rows.items()}
    ext = [M for M, r in rows.items() if r["ext"] or r["out"].startswith("EXTENSION")]
    if any(o.startswith("EXTENSION REQUIRED") for o in outs.values()): verdict = "EXTENSION REQUIRED at M = " + ", ".join(str(M) for M in ext)
    elif all(o == "NO DIFFERENCE" for o in outs.values()): verdict = "PASS"
    elif any(o == "DIFFERENCE" for o in outs.values()): verdict = "FAIL (DIFFERENCE)"
    else: verdict = "FAIL (UNRESOLVED)"
    print(f"\nM = 300: {outs[300]}; M = 1500: {outs[1500]}" + (f"; the extension was run at M = {', '.join(map(str, ext))} (disclosed)" if ext else "; no extension"))
    print(f"95 % intervals of the relative difference (the stated bound): " + "; ".join(
        f"M = {M}: [{(r['ext'] or r)['lo']:+.3f}, {(r['ext'] or r)['hi']:+.3f}] %" for M, r in rows.items()))
    print(f"TEST G (Mac){' [dry run]' if DRY else ''}: {verdict}")
    if DRY: return
    print("\n## Information (no verdict): nu_d on the main masses\n"); information_nud(D)
    print("\n## Information (no verdict): the other rows\n"); information_rows()


if __name__ == "__main__":
    main()
