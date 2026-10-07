#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.13, GATE V3, third plan-author decision of 2026-10-07, item 2): TEST T-PRIME, written and
committed before any of its data. A new test, not an extension of Test T.
Test T-prime: cell epi8_H_H10_L10, method B, the campaign protocol; masses M = 300 and M = 1500 only; 400 fresh seeds per policy
per mass (cluster/resched_gate_261005/tasks_Tprime_epi8_H_H10_L10.txt, gen_testTprime_261007.py); the same seeds on both policies
(minimal, --legacy-resched) of the binary 7b08827 that ran Test T (testTprime.sbatch checks its recorded hash and version line);
same partition and interleaving as Test T; HD_CONTACT_AUDIT=1; node recorded per trajectory.
Registered estimator: the mean nu per mass (red_nu.csv, the canonical argmax estimator of reduce_B.py), SE = SD/sqrt(n),
z = (mean_minimal - mean_legacy)/sqrt(SE_min^2 + SE_leg^2); 95 % interval of the relative difference 100 (d -+ 1.96 s)/mean_legacy
(Test T's formula, 1.96 = Phi^-1(0.975) rounded as there).
RULE (plan author, verbatim in 261012 sec. 4.4.13):
  M = 300 (the verdict number): BIAS CONFIRMED if z >= 3; NO BIAS if |z| < 2 AND the 95 % interval of the relative difference
  excludes +0.475 %; anything else, including z <= -2, = FAIL.
  M = 1500: FAIL if |z| >= 3; otherwise information (expected z if the Test T size -0.239 % were real: printed by --design).
  Test T-prime PASS = M = 300 NO BIAS, M = 1500 not FAIL, inventory clean (every task-list trajectory present, health 0, policy
  lines right, contact <= 1e-6 px, one build). Gate v3 also needs ASan CLEAN (asan.sbatch: the fetched report's last line).
STOPPING RULE (recorded before the data): T-prime is the last statistical test for this fix. FAIL or unresolved = the fix is
shelved, 279282b stays, no further test. PASS = the build is ACCEPTED FOR THE FLUID REGIME on the combined record (Test T for ten
numbers, T-prime for the eleventh, nu at M = 300; the second chance disclosed), with the T and T-prime 95 % intervals as the stated
bias bound. Use in the melting window still needs the same-binary A/B at N = 100 inside the window (melting stage 1).
usage (from hspist3/):
  python3 validation/resched_testTprime_261007.py --design      thresholds and false-fail rates, power, seeds, cost, --time
  python3 validation/resched_testTprime_261007.py --dry-run-T   the rule on Test T's M = 300 and 1500 data (must give +3.02, -2.56)
  python3 validation/resched_testTprime_261007.py               the analysis of the T-prime data (experiments_resched_gate2_261007/)
"""
import contextlib, glob, hashlib, io, math, os, re, sys
from collections import Counter
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = os.path.join(HS, "cluster", "resched_gate_261005")
sys.path.insert(0, HERE); sys.path.insert(0, HS); sys.path.insert(0, GATE)
import paper1_confinement_results_261004 as R
import tests_20260913 as T
import resched_testT_261007 as TT
import gen_testTprime_261007 as G

CID, MASSES, NSEED = G.CID, G.MASSES, G.NSEED
LOC = TT.LOC
REL_TP = G.REL
TASKS_TP = os.path.join(GATE, f"tasks_Tprime_{CID}.txt")
TASKS_T = os.path.join(GATE, f"tasks_T_{CID}.txt")
Z_CONF, Z_NOBIAS, SIZE_300, Z_1500, SIZE_1500, Z95 = 3.0, 2.0, 0.475, 3.0, -0.239, 1.96
SHIFTS = (0.0, 0.2, 0.3, 0.475)
POL = ("minimal", "legacy")
DESIGN, DRY = "--design" in sys.argv, "--dry-run-T" in sys.argv


def cell_dir(pol, M, which):
    rel = TT.REL_T if which == "T" else REL_TP
    return os.path.join(LOC, rel, pol, CID, f"m_{M}")


def load(which):
    out = {}
    for pol in POL:
        for M in MASSES:
            d = cell_dir(pol, M, which)
            out[(pol, M)] = dict(dir=d, red=pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True),
                                 log=open(os.path.join(d, "run.log"), errors="ignore").read(),
                                 failed=[open(f, errors="ignore").read() for f in glob.glob(os.path.join(d, ".failed_run*", "stdout.log"))],
                                 build=(open(os.path.join(d, ".build_git")).read().strip() if os.path.exists(os.path.join(d, ".build_git")) else None))
    return out


def rule_300(z, lo, hi):
    if z >= Z_CONF: return "BIAS CONFIRMED"
    if abs(z) < Z_NOBIAS and not (lo <= SIZE_300 <= hi): return "NO BIAS"
    return "FAIL"


def compare(D):
    """The registered numbers per mass: means, SEs, difference, relative difference with its 95 % interval, z, outcome."""
    rows = {}
    for M in MASSES:
        a = D[("minimal", M)]["red"]["nu"].to_numpy(float); b = D[("legacy", M)]["red"]["nu"].to_numpy(float)
        sa, sb = a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b))
        d = a.mean() - b.mean(); s = math.hypot(sa, sb); z = d / s
        lo, hi = 100 * (d - Z95 * s) / b.mean(), 100 * (d + Z95 * s) / b.mean()
        out = rule_300(z, lo, hi) if M == 300 else ("FAIL" if not math.isfinite(z) or abs(z) >= Z_1500 else "information")   # no z = unresolved
        rows[M] = dict(na=len(a), nb=len(b), ma=a.mean(), sa=sa, mb=b.mean(), sb=sb, d=d, rel=100 * d / b.mean(), lo=lo, hi=hi, z=z, out=out)
    return rows


def print_rule(rows):
    print("| mass | n minimal | n legacy | minimal mean nu | SE | legacy mean nu | SE | difference | relative [%] | 95 % interval [%] | "
          "z | outcome by the rule |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    for M, r in rows.items():
        print(f"| {M} | {r['na']} | {r['nb']} | {r['ma']:.7f} | {r['sa']:.2e} | {r['mb']:.7f} | {r['sb']:.2e} | {r['d']:+.3e} | {r['rel']:+.3f} | "
              f"[{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} | **{r['out']}** |")


def design():
    print("# Test T-prime -- design numbers (gate v3, 261012 sec. 4.4.13), before any of its data\n")
    pz = norm.sf(Z_CONF); p2 = 2 * norm.sf(Z_NOBIAS); p15 = 2 * norm.sf(Z_1500)
    D = load("T"); sd = {(p, M): D[(p, M)]["red"]["nu"].std(ddof=1) for p in POL for M in MASSES}
    nul = {M: D[("legacy", M)]["red"]["nu"].mean() for M in MASSES}
    se = {M: math.sqrt(sd[("minimal", M)] ** 2 + sd[("legacy", M)] ** 2) / math.sqrt(NSEED) / nul[M] * 100 for M in MASSES}   # % of nu
    zci = SIZE_300 / se[300] - Z95
    print("### Thresholds and their false-fail rates under the null (no shift at either mass)\n")
    print("| threshold | role | probability under the null |\n|---|---|---|")
    print(f"| M = 300: z >= {Z_CONF:g} | BIAS CONFIRMED (false confirmation) | {pz:.5f} |")
    print(f"| M = 300: abs(z) >= {Z_NOBIAS:g} or the 95 % interval contains +{SIZE_300} % | FAIL or CONFIRMED (not NO BIAS) | "
          f"{p2:.5f} (the interval condition binds only above z = {zci:.2f}, outside abs(z) < 2 at the design SE) |")
    print(f"| M = 1500: abs(z) >= {Z_1500:g} | FAIL | {p15:.5f} |")
    tot = 1 - (1 - p2) * (1 - p15)
    print(f"| Test T-prime as a whole | not PASS (false fail), the two masses independent | {tot:.5f} |\n")
    print("### Expected precision with 400 + 400 seeds per mass (per-seed SDs of Test T, argmax estimator)\n")
    print("| M | per-seed SD minimal | per-seed SD legacy | legacy mean nu (Test T) | SE of the relative difference [%] | 95 % half-width [%] |\n|---|---|---|---|---|---|")
    for M in MASSES:
        print(f"| {M} | {sd[('minimal', M)]:.3e} | {sd[('legacy', M)]:.3e} | {nul[M]:.7f} | {se[M]:.4f} | {Z95 * se[M]:.4f} |")
    print(f"\nexpected z if the Test T size were real: M = 300 +{SIZE_300} % -> {SIZE_300 / se[300]:+.2f}; M = 1500 {SIZE_1500} % -> {SIZE_1500 / se[1500]:+.2f}")
    zs_T = norm.isf(TT.FW / (2 * TT.NNUM)); se_T = se[300] * math.sqrt(NSEED / len(D[("legacy", 300)]["red"]))   # Test T: 100 + 100
    mu_T = SIZE_300 / se_T
    print(f"for comparison, Test T at M = 300 ({len(D[('legacy', 300)]['red'])} + {len(D[('minimal', 300)]['red'])} seeds, same SDs): SE {se_T:.4f} %, "
          f"expected z for +{SIZE_300} % = {mu_T:+.2f}, power P(abs(z) >= z* = {zs_T:.4f}) = {norm.sf(zs_T - mu_T) + norm.cdf(-zs_T - mu_T):.3f}")
    print("\n### Power at M = 300 (n = 400 + 400): probability of each outcome for a true shift\n")
    print("| true shift [%] | expected z | P(BIAS CONFIRMED) | P(NO BIAS) | P(FAIL) |\n|---|---|---|---|---|")
    for sh in SHIFTS:
        mu = sh / se[300]; pc = norm.sf(Z_CONF - mu); pn = max(0.0, norm.cdf(min(Z_NOBIAS, zci) - mu) - norm.cdf(-Z_NOBIAS - mu))
        print(f"| {sh:+.3f} | {mu:+.2f} | {pc:.4f} | {pn:.4f} | {1 - pc - pn:.4f} |")
    mu = SIZE_1500 / se[1500]
    print(f"\nM = 1500: P(FAIL) = P(abs(z) >= {Z_1500:g}) = {p15:.4f} with no shift; {norm.sf(Z_1500 - mu) + norm.cdf(-Z_1500 - mu):.4f} if {SIZE_1500} % were real")
    print("\n### Fresh seeds\n")
    S = G.seeds(); vals = list(S.values()); used, tset = set(), set()
    for f in glob.glob(os.path.join(HS, "cluster", "confinement_20261013", "tasks_*.txt")):
        for l in open(f):
            p = l.split()
            if p[0] == "B": used.add(int(p[4]))
            elif p[0] in ("A", "AF"): used.add(int(p[3]))
    for l in open(TASKS_T): tset.add(int(l.split()[4]))
    pilot = {T.run_seed(20261013, 0, mi, r) for mi in range(9) for r in range(4)}
    print(f"seeds: run_seed({G.BASE}, {G.STREAM}, mass index, r), r = 0..{NSEED - 1}, M = 300 and 1500: {len(vals)} seeds, {len(set(vals))} distinct")
    print(f"overlap with: the campaign task lists (B/A/AF; the replay re-ran these) {len(set(vals) & used)} of {len(used)}; Test T "
          f"{len(set(vals) & tset)} of {len(tset)}; smoke/pilot/E2 seeds run_seed(20261013, 0, m, r < 4) {len(set(vals) & pilot)}; "
          f"A-fixed 9700-9703 {len(set(vals) & {9700, 9701, 9702, 9703})}")
    print(f"SHA-256 of the seed list (lines 'M r seed', M ascending, r ascending): {hashlib.sha256(G.seed_text(S).encode()).hexdigest()}")
    print(f"SHA-256 of the task list {os.path.relpath(TASKS_TP, HS)} ({sum(1 for _ in open(TASKS_TP))} lines): "
          f"{hashlib.sha256(open(TASKS_TP, 'rb').read()).hexdigest()}")
    print("\n### Cost from the measured KOA times of Test T (same binary, same cell; seconds per trajectory)\n")
    print("| M | minimal [s] | legacy [s] | trajectories | core-h | array task wall time on 8 cores [min] |\n|---|---|---|---|---|---|")
    tot = 0.0; walls = {}
    for M in MASSES:
        t = {p: [int(v) for v in re.findall(r"##RUN .*?\((\d+) s\)", D[(p, M)]["log"])] for p in POL}
        ch = NSEED * (np.mean(t["minimal"]) + np.mean(t["legacy"])) / 3600; tot += ch; walls[M] = ch * 60 / 8
        print(f"| {M} | {np.mean(t['minimal']):.1f} | {np.mean(t['legacy']):.1f} | {2 * NSEED} | {ch:.2f} | {walls[M]:.1f} |")
    tl = math.ceil(2 * max(walls.values()) / 10) * 10
    print(f"\ntotal {tot:.1f} core-h for {2 * NSEED * len(MASSES)} trajectories; --time = 2 x the longest task ({max(walls.values()):.1f} min, "
          f"M = {max(walls, key=walls.get)}) rounded up to 10 min = {tl // 60}:{tl % 60:02d}:00; two array tasks, %2, 8 cores each = 16 cores")


def inventory(D):
    want = {}
    for l in open(TASKS_TP):
        f = l.split(); want.setdefault((f[10], int(f[2])), set()).add((int(f[3]), int(f[4])))
    print("| policy | M | trajectories matched (expected) | n finite > 0 | log sections | health lines (failed runs incl.) | policy line wrong | "
          "max contact gap [px] | contact lines missing | nodes | ok |\n|---|---|---|---|---|---|---|---|---|---|---|")
    ok_all = True; builds = set(); nodes = {p: Counter() for p in POL}
    for pol in POL:
        for M in MASSES:
            e = D[(pol, M)]; red = e["red"]; have = {(int(a), int(b)) for a, b in zip(red["run"], red["seed"])}
            exp = want.get((pol, M), set()); matched = len(exp & have) if len(have) == len(red) else 0
            nfin = int((np.isfinite(red["n"].astype(float)) & (red["n"].astype(float) > 0)).sum())
            secs = e["log"].split("##RUN")[1:]
            nh = sum(s.count("[EDMD-HEALTH]") for s in secs) + sum(f.count("[EDMD-HEALTH]") for f in e["failed"])
            want_line = "minimal (default)" if pol == "minimal" else "legacy (full reschedule"
            bad = sum(f"[EDMD-RESCHED] divider events: {want_line}" not in s for s in secs)
            cg, nmiss = 0.0, 0
            for s in secs:
                m = TT.CONTACT.search(s)
                if not m: nmiss += 1; continue
                cg = max(cg, max(float(v) for v in m.groups()[1:]))
                n = TT.NODE.search(s.split("\n", 1)[0])
                if n: nodes[pol][n.group(3)] += 1
            if e["build"]: builds.add(e["build"])
            ok = matched == len(exp) == NSEED and nfin == len(red) and len(secs) == len(red) and nh == 0 and bad == 0 and nmiss == 0 and cg <= 1e-6
            ok_all &= ok
            print(f"| {pol} | {M} | {matched} ({len(exp)}) | {nfin} | {len(secs)} | {nh} | {bad} | {cg:.2e} | {nmiss} | {len(nodes[pol])} | {'yes' if ok else '**NO**'} |")
    ok_all &= builds == {"00ALLINONE  git 7b08827  target koa"}
    print(f"\nnodes: minimal {dict(nodes['minimal'])}; legacy {dict(nodes['legacy'])}; shared: {sorted(set(nodes['minimal']) & set(nodes['legacy']))}")
    print(f"builds (.build_git): {sorted(builds)} (required: the Test T binary, 00ALLINONE  git 7b08827  target koa); inventory {'clean' if ok_all else '**NOT CLEAN**'}")
    return ok_all


def asan_status():
    reps = sorted(glob.glob(os.path.join(LOC, "asan_261007_*", "report.txt")))
    if not reps: return "MISSING", []
    last = [open(p).read().strip().split("\n")[-1] for p in reps]
    clean = all("): CLEAN" in l for l in last)
    return ("CLEAN" if clean else "NOT CLEAN"), list(zip([os.path.relpath(p, HS) for p in reps], last))


def combined_record(Dp):
    """Test T's eleven numbers (ten from T) and T-prime's two masses: the 95 % intervals that are the stated bias bound on a PASS."""
    with contextlib.redirect_stdout(io.StringIO()):
        c = {x["cid"]: x for x in R.cells()}[CID]; co = dict(c); R.method_B(co)
    Ms, Ns, x = c["Ms"], c["Ns"], np.asarray(co["x"], float)
    nu = {p: {M: pd.read_csv(os.path.join(LOC, TT.REL_T, p, CID, f"m_{M}", "red_nu.csv"))["nu"].to_numpy(float) for M in Ms} for p in POL}
    print("| number | source | relative difference [%] | 95 % interval [%] | z |\n|---|---|---|---|---|")
    for name, a, sa, b, sb in TT.eleven(nu["minimal"], nu["legacy"], Ms, x, Ns):
        d = a - b; s = math.hypot(sa, sb); src = "Test T" + (" (superseded by T-prime)" if name == "nu M=300" else "")
        print(f"| {name} | {src} | {100 * d / b:+.3f} | [{100 * (d - Z95 * s) / b:+.3f}, {100 * (d + Z95 * s) / b:+.3f}] | {d / s:+.2f} |")
    for M, r in compare(Dp).items():
        print(f"| nu M={M} | Test T-prime | {r['rel']:+.3f} | [{r['lo']:+.3f}, {r['hi']:+.3f}] | {r['z']:+.2f} |")


def information(D):
    import warnings; warnings.filterwarnings("ignore")
    import resched_testT_followup_261007 as F
    print("| M | estimator | minimal mean | legacy mean | relative [%] | z | same-seed correlation |\n|---|---|---|---|---|---|---|")
    for M in MASSES:
        ref = 0.5 * (D[("minimal", M)]["red"]["nu"].mean() + D[("legacy", M)]["red"]["nu"].mean()); v = {}
        for pol in POL:
            dd = D[(pol, M)]["dir"]; red = D[(pol, M)]["red"]; z = np.load(os.path.join(dd, "acf_runs.npz")); dt = float(red["dt"].iloc[0]); w = []
            for r in red["run"]:
                try: w.append(F.fit_acf(z[f"run{r}"].astype(float), dt, ref)[4] / (2 * math.pi))
                except Exception: w.append(np.nan)
            v[pol] = np.array(w)
        for lab, a, b in (("argmax (registered)", D[("minimal", M)]["red"]["nu"].to_numpy(float), D[("legacy", M)]["red"]["nu"].to_numpy(float)),
                          ("nu_d (sec. 4.4.13 item 4a)", v["minimal"], v["legacy"])):
            ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
            s = math.hypot(a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b)))
            print(f"| {M} | {lab} | {a.mean():.7f} | {b.mean():.7f} | {100 * (a.mean() - b.mean()) / b.mean():+.3f} | {(a.mean() - b.mean()) / s:+.2f} | "
                  f"{np.corrcoef(a, b)[0, 1]:+.3f} |")


def main():
    if DESIGN: return design()
    which = "T" if DRY else "Tprime"
    print(f"# Test T-prime (gate v3, 261012 sec. 4.4.13){'  [DRY RUN on Test T data, n = 100 + 100: a test of this script, no verdict]' if DRY else ''}\n")
    print(f"RULE: M = 300 BIAS CONFIRMED if z >= {Z_CONF:g}; NO BIAS if abs(z) < {Z_NOBIAS:g} and the 95 % interval excludes +{SIZE_300} %; "
          f"else FAIL. M = 1500: FAIL if abs(z) >= {Z_1500:g}, else information.\n")
    D = load(which)
    if not DRY:
        print("## Inventory\n"); inv = inventory(D)
    rows = compare(D)
    print("\n## The registered numbers\n"); print_rule(rows)
    if DRY:
        ok = round(rows[300]["z"], 2) == 3.02 and round(rows[1500]["z"], 2) == -2.56
        print(f"\ndry run: z(M = 300) = {rows[300]['z']:+.2f}, z(M = 1500) = {rows[1500]['z']:+.2f}; Test T printed +3.02 and -2.56: "
              f"{'REPRODUCED' if ok else '**NOT REPRODUCED**'}")
        print("\n## The combined record (format test)\n"); combined_record(D); return
    tp = rows[300]["out"] == "NO BIAS" and rows[1500]["out"] != "FAIL" and inv
    print(f"\nM = 300: {rows[300]['out']}; M = 1500: {rows[1500]['out']} (z = {rows[1500]['z']:+.2f}); inventory {'clean' if inv else 'NOT CLEAN'}")
    print(f"TEST T-PRIME: {'PASS' if tp else 'FAIL'}")
    st, reps = asan_status()
    for p, l in reps: print(f"ASan report {p}: {l}")
    print(f"ASan (item 3): {st}")
    if st == "MISSING": verdict = "PENDING -- the ASan report has not been fetched (fetch_resched2.sh)"
    elif tp and st == "CLEAN": verdict = "ACCEPTED FOR THE FLUID REGIME (combined record below; second chance disclosed; E0, E1, E2 from gate v2, sec. 4.4.12)"
    else: verdict = "NOT ACCEPTED -- by the stopping rule the fix is shelved and 279282b stays; no further statistical test"
    print(f"\nGATE V3: {verdict}")
    print("\n## The combined record: Test T for ten numbers, T-prime for nu at M = 300 (and M = 1500 as information) -- the stated bias bound on a PASS\n")
    combined_record(D)
    print("\n## Information (no verdict): the refined frequency nu_d of 261012 sec. 4.4.13 item 4a on the T-prime trajectories\n")
    information(D)


if __name__ == "__main__":
    main()
