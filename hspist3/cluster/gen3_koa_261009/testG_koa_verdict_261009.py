#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G), REGISTRATION of TEST G ON KOA by plan-author decision 13, item 2g (261012
sec. 4.7.29): the verdict of Test G on KOA. Registered with the task list (gen_testG_koa_261009.py, tasks_testG_koa_261009.txt and its
SHA-256) before any KOA trajectory exists; BUILD_KOA_COMMIT (the commit the KOA binary is built from) is set in the commit that follows
the registration commit and names it; nothing is edited afterwards.
THE RULE (sec. 4.7.15, the registered Mac verdict validation/testG_verdict_261009.py, imported UNCHANGED and used with three masses):
per mass M = 50, 300 and 1500 the mean nu per engine (red_nu.csv, reduce_B.py's canonical argmax estimator), SE = SD/sqrt(n),
z = (mean_gen3 - mean_gen2) / sqrt(SE_gen3^2 + SE_gen2^2): NO DIFFERENCE if |z| < 2; DIFFERENCE if |z| >= 3; otherwise ONE extension of
400 more trajectories per engine at that mass (r = 400..799 of the same stream), judged on all 800 with the same thresholds; still
2 <= |z| < 3 -> UNRESOLVED. TEST G (KOA) PASS = NO DIFFERENCE at all three masses.
INVENTORY (else INVALID, no verdict): the Mac script's inventory (all trajectories of the main list and of every extension that was run;
every gen3 run record clean=1, one build line and one record per trajectory; no gen2 health line; no failed trajectory; one build),
with the masses set to the three and the build line set to the KOA binary's: "00ALLINONE  git <h>  target koa" where <h> is the
abbreviated hash the KOA clone printed (git rev-parse --short: at least 7 hex digits, a prefix of BUILD_KOA_COMMIT).
AN EXTENSION is counted for a mass if the fetched data hold any of its runs r >= 400; its task list (written by
gen_testG_koa_261009.py --extension M from the registered seeds) must then be in this folder, and all its runs must be present.
ALSO PRINTED: the false-alarm rates of the thresholds for one and for three independent masses (the Mac script's rates(): exact normal
integrals, the extension path included); the smallest true difference with z = 3; the 95 % intervals. INFORMATION rows (no verdict):
nu_d on the three verdict masses; the Mac's other rows WITHOUT the M = 50 row (M = 2000, the dense cell, N = 400 at pi/8, the static
method).
usage (from hspist3/, on the Mac after the fetch):
  python3 cluster/gen3_koa_261009/testG_koa_verdict_261009.py               the analysis of the fetched KOA data
  python3 cluster/gen3_koa_261009/testG_koa_verdict_261009.py --design      the false-alarm rates for three masses (no data)
  python3 cluster/gen3_koa_261009/testG_koa_verdict_261009.py --dry-run-mac a test of this script on the Mac's Test G data: the rule
                                                                            table and the verdict line with three masses (M = 50 has
                                                                            100 per engine there); no inventory, no verdict
"""
import glob, math, os, re, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation"))
BUILD_KOA_COMMIT = "<SET IN THE COMMIT AFTER THE REGISTRATION>"
MASSES_KOA = (50, 300, 1500)
import testG_verdict_261009 as V
import numpy as np, pandas as pd
V.MASSES = MASSES_KOA                                   # the Mac script's loops read its module-level MASSES at call time
LOC_KOA = os.path.join(HS, "experiments_gen3_gate_koa_261009")
TASKS_KOA = os.path.join(HERE, "tasks_testG_koa_261009.txt")


def print_rates3():
    R = V.rates()
    print("| quantity (no true difference) | one mass | three masses (independent) |\n|---|---|---|")
    print(f"| P(abs(z) >= {V.Z_D:g} at 400) | {R['d1']:.5f} | {1 - (1 - R['d1']) ** 3:.5f} (at least one) |")
    print(f"| P(extension: {V.Z_ND:g} <= abs(z) < {V.Z_D:g} at 400) | {R['ext']:.5f} | {1 - (1 - R['ext']) ** 3:.5f} (at least one) |")
    print(f"| P(DIFFERENCE), extension path included | {R['diff']:.5f} | {1 - (1 - R['diff']) ** 3:.5f} (at least one) |")
    print(f"| P(UNRESOLVED) | {R['unres']:.5f} | {1 - (1 - R['unres']) ** 3:.5f} (at least one) |")
    print(f"| P(NO DIFFERENCE) | {R['nd']:.5f} | P(PASS) = {R['nd'] ** 3:.5f}; false alarm (not PASS) = {1 - R['nd'] ** 3:.5f} |")


def information_rows3():
    """the Mac script's information_rows() without its M = 50 row (that mass is a verdict mass here); otherwise its code"""
    print("| row | n gen3 | n gen2 | gen3 mean | SE | gen2 mean | SE | relative [%] | 95 % interval [%] | z | clean (both engines) |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|")
    for lab, cid, M in (("nu, M = 2000", V.CID, 2000), ("nu, dense cell, M = 300", "dense_H10_L5.604167", 300),
                        ("nu, N = 400 at pi/8, M = 300", "epi8_N400_H20_L20", 300)):
        c = {e: V.load_cell(e, cid, M) for e in V.ENG}
        if any(c[e]["red"] is None for e in V.ENG): print(f"| {lab} | missing | | | | | | | | | |"); continue
        a, b = c["gen3"]["red"]["nu"].to_numpy(float), c["gen2"]["red"]["nu"].to_numpy(float); r = V.stats(a, b)
        fl = sum(len(c[e]["failed"]) for e in V.ENG)
        recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", c["gen3"]["log"], re.M)
        g2h = sum(len(V.HEALTH.findall(s)) for s in V.sections(c["gen2"]["log"]))
        cl = f"failed {fl}; gen3 records {len(recs)}, clean=1 {recs.count('1')}; gen2 health lines {g2h}"
        print(f"| {lab} | {r['na']} | {r['nb']} | {r['ma']:.7f} | {r['sa']:.2e} | {r['mb']:.7f} | {r['sb']:.2e} | {r['rel']:+.3f} | "
              f"[{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | {cl} |")
    for face in ("F_L", "F_R"):
        v = {}
        for e in V.ENG:
            fs = sorted(glob.glob(os.path.join(V.LOC, V.G.REL, e, V.CID, "x_0", "red_*.csv")))
            v[e] = np.array([float(pd.read_csv(f)[face].iloc[0]) for f in fs]) if fs else np.array([])
        if min(len(v[e]) for e in V.ENG) < 2: print(f"| A-fixed {face} | missing | | | | | | | | | |"); continue
        r = V.stats(v["gen3"], v["gen2"])
        nf = {e: len(glob.glob(os.path.join(V.LOC, V.G.REL, e, V.CID, "x_0", ".failed_*"))) for e in V.ENG}
        print(f"| static method, A-fixed at x_0, mean {face} | {r['na']} | {r['nb']} | {r['ma']:.6f} | {r['sa']:.2e} | {r['mb']:.6f} | {r['sb']:.2e} | "
              f"{r['rel']:+.3f} | [{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | failed gen3 {nf['gen3']}, gen2 {nf['gen2']} |")


def verdict_line(rows, dry):
    outs = {M: V.final(r) for M, r in rows.items()}
    ext = [M for M, r in rows.items() if r["ext"] or r["out"].startswith("EXTENSION")]
    if any(o.startswith("EXTENSION REQUIRED") for o in outs.values()): verdict = "EXTENSION REQUIRED at M = " + ", ".join(str(M) for M in ext)
    elif all(o == "NO DIFFERENCE" for o in outs.values()): verdict = "PASS"
    elif any(o == "DIFFERENCE" for o in outs.values()): verdict = "FAIL (DIFFERENCE)"
    else: verdict = "FAIL (UNRESOLVED)"
    print("\n" + "; ".join(f"M = {M}: {o}" for M, o in outs.items()) + (f"; the extension was run at M = {', '.join(map(str, ext))} (disclosed)" if ext else "; no extension"))
    print("95 % intervals of the relative difference (the stated bound): " + "; ".join(
        f"M = {M}: [{(r['ext'] or r)['lo']:+.3f}, {(r['ext'] or r)['hi']:+.3f}] %" for M, r in rows.items()))
    print(f"TEST G (KOA){' [dry run on the Mac data: a test of this script, no verdict]' if dry else ''}: {verdict}")


def main():
    if "--design" in sys.argv:
        print("# Test G on KOA -- the false-alarm rates of the thresholds for three masses (261012 sec. 4.7.29)\n"); print_rates3(); return
    dry = "--dry-run-mac" in sys.argv
    print(f"# Test G on KOA (261012 sec. 4.7.29; decision 13, item 2g){'  [DRY RUN on the Mac Test G data: a test of this script]' if dry else ''}\n")
    print(f"RULE per mass (M = 50, 300, 1500): NO DIFFERENCE if abs(z) < {V.Z_ND:g}; DIFFERENCE if abs(z) >= {V.Z_D:g}; otherwise one extension (400 more per "
          f"engine), judged on 800 with the same thresholds, still between -> UNRESOLVED. PASS = NO DIFFERENCE at all three masses.\n")
    if dry:
        D = {(e, M): V.load_cell(e, V.CID, M) for e in V.ENG for M in MASSES_KOA}
    else:
        if BUILD_KOA_COMMIT.startswith("<"): sys.exit("STOP: BUILD_KOA_COMMIT is set in the commit after the registration")
        V.LOC, V.TASKS = LOC_KOA, TASKS_KOA
        D = {(e, M): V.load_cell(e, V.CID, M) for e in V.ENG for M in MASSES_KOA}
        builds = {D[k]["build"] for k in D if D[k]["build"]}
        m = re.fullmatch(r"00ALLINONE  git ([0-9a-f]{7,40})  target koa", sorted(builds)[0]) if len(builds) == 1 else None
        V.BUILD = sorted(builds)[0] if (m and BUILD_KOA_COMMIT.startswith(m.group(1))) else f"00ALLINONE  git {BUILD_KOA_COMMIT[:7]}  target koa"
        ran = [M for M in MASSES_KOA if any(D[(e, M)]["red"] is not None and (D[(e, M)]["red"]["run"] >= V.NSEED).any() for e in V.ENG)]
        use = os.path.join(LOC_KOA, ".extension_lists_counted"); os.makedirs(use, exist_ok=True)
        for M in ran:
            src = os.path.join(HERE, f"tasks_testG_ext_M{M}_261009.txt")
            if not os.path.exists(src): sys.exit(f"STOP: runs r >= 400 at M = {M} but no extension list {os.path.basename(src)}")
            dst = os.path.join(use, os.path.basename(src))
            if not os.path.exists(dst): shutil.copyfile(src, dst)
            if open(dst, "rb").read() != open(src, "rb").read(): sys.exit(f"STOP: {dst} differs from {src}")
        V.GG = use                                      # the inventory reads exactly the extension lists that were run
        print(f"extensions counted (runs r >= 400 in the data): {ran if ran else 'none'}\n")
        print("## Inventory\n"); inv = V.inventory(D)
        if not inv:
            print("\nTEST G (KOA): INVALID -- the inventory is not clean; no verdict"); return
    rows = V.compare(D)
    print("\n## The registered numbers\n"); V.print_rule(rows)
    print("\n## The thresholds' false-alarm rates (no true difference)\n"); print_rates3()
    print("\n## The smallest true difference that would give z = 3 with the observed SEs\n")
    print("| mass | judged on | SE of the difference | 3 x SE | relative to gen2's mean [%] |\n|---|---|---|---|---|")
    for M, r0 in rows.items():
        r = r0["ext"] or r0
        print(f"| {M} | {r['judged_on']} | {r['s']:.2e} | {3 * r['s']:.2e} | {100 * 3 * r['s'] / r['mb']:.3f} |")
    verdict_line(rows, dry)
    if dry: return
    print("\n## Information (no verdict): nu_d on the verdict masses\n"); V.information_nud(D)
    print("\n## Information (no verdict): the other rows (the Mac's, without M = 50)\n"); information_rows3()


if __name__ == "__main__":
    main()
