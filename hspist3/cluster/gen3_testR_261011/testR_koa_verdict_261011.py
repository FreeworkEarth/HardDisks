#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28; 261012 sec. 4.7.57): REGISTRATION of TEST R ON KOA -- its verdict. Committed with the
task list before chunk 1 runs; not edited afterwards. ONE PRINTING: the analysis of the fetched data writes its output file once
(experiments_gen3_testR_koa_261011/testR_verdict_output.txt) and refuses to run again when that file exists.
THE DATA (fetch_testR_koa.sh): experiments_gen3_testR_koa_261011/merged/.merge_complete (the merge record) and per engine and mass
merged/testR/<gen2|gen3>/epi8_H_H10_L10/m_<M>/: red_nu.csv (reduce_B.py, registered), nud_runs.csv, runs_health.csv and
failed_runs.csv (reduce_testR_koa.py on KOA), acf_runs.npz (information rows only), .build_git.
THE RULE (decision 28, verbatim in substance):
  estimators per trajectory: nu_arg = reduce_B.py's argmax (red_nu.csv's nu); nu_d = fit_acf as testG_koa_posthoc_261011.py at 60909289
    computes it, started at the mass's pooled argmax mean of its own engine (nud_runs.csv, fit_ok = 1; a fit that failed is left out
    of the nu_d mean and counted);
  statistic per mass and estimator: d = mean(gen3) / mean(gen2) - 1; SE(d) from each arm's seed scatter, unpaired:
    SE(d)^2 = (SE3 / m2)^2 + (m3 SE2 / m2^2)^2 with SE = SD / sqrt(n); z = d / SE(d);
  per mass: DIFFERENCE CONFIRMED if |z_arg| >= 3.5 or |z_d| >= 3.5; NO DIFFERENCE if |z_arg| < 3.5 and |z_d| < 3.5 and the 95 % interval
    of d from nu_d, d_d +- 1.96 SE(d_d), lies within +-0.15 % (inclusive); NOT DECIDED otherwise;
  overall: TEST R: ENGINES AGREE if NO DIFFERENCE at all three masses; TEST R: ENGINES DIFFER if DIFFERENCE CONFIRMED at any mass;
    TEST R: NOT DECIDED otherwise.
  health (Test G's rule, testG_worker_koa.sh's): gen2 -- no line of the HEALTH pattern, the minimal-policy line, no [EDMD3 line; gen3 --
    exactly one build line, exactly one run record with clean=1 and engine=gen3, no other HEALTH line; both -- the contact-audit line,
    its largest |contact distance| <= 1e-6 px (Test G's inventory). EXCLUDED and counted, never rerun: a run the worker refused
    (failed_runs.csv), a run in the tables that fails the rule, a registered run absent from both (listed apart). More than 1 %
    excluded in any (engine, mass) arm (more than 36 of 3600): TEST R: NOT DECIDED (this script's reading of decision 28: the overall
    verdict, whatever the masses print).
  INVALID, no verdict (Test G's precedent; a pipeline fault, not a data outcome): no merge record; a build line other than
    "00ALLINONE  git a643bc7  target koa"; a row whose (run, seed) is not a registered task of its arm, or a run twice; nu_arg of
    nud_runs.csv not equal to red_nu.csv's (to 1e-14 relative: the two files' decimal texts may read back 1 ulp apart); the task list
    not the registered one.
INFORMATION (no verdict): tau_r, Q, nu_0 and delta_M per engine as in Report G2 (testG_koa_posthoc_261011.py's functions; jackknife
groups of 10; start value the arm's own argmax mean; delta_M against c_ref = c_300, the one mass >= 300); the argmax histograms; d from
Test G and Test R pooled (M = 50 and 300; Test G's nu_d recomputed by the same rule from its fetched ACFs); the fit-failure counts; a
Mac recomputation of nu_d for the first 20 trajectories of every arm against KOA's table.
--design (before any data): every threshold; the per-mass and overall false-alarm rates; the probability of each overall outcome for
true differences of 0, 0.1, 0.2 and 0.3 % at M = 50 (none at 100 and 300), by Monte Carlo (default_rng(20261317), 2e6 draws per
scenario), with the per-trajectory scatter of Test G on KOA (SD of nu_arg and nu_d per engine, their correlation; M = 100
interpolated in log M between 50 and 300) and 3600 trajectories per arm.
--selftest: the formulas, the rule and the inventory on synthetic data with known answers.
usage (from hspist3/):
  python3 cluster/gen3_testR_261011/testR_koa_verdict_261011.py --design
  python3 cluster/gen3_testR_261011/testR_koa_verdict_261011.py --selftest
  python3 cluster/gen3_testR_261011/testR_koa_verdict_261011.py              (the one printing, after the fetch)
"""
import hashlib, io, math, os, re, sys, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS); sys.path.insert(0, os.path.join(HS, "cluster", "gen3_koa_261009"))
warnings.filterwarnings("ignore")
MASSES, ENG, CID, NSEED = (50, 100, 300), ("gen2", "gen3"), "epi8_H_H10_L10", 3600
Z_CONF, BAND, Z95, EXCL = 3.5, 0.0015, 1.959963984540054, 0.01
BUILD = "00ALLINONE  git a643bc7  target koa"
TASKS = os.environ.get("HD_TESTR_TASKS", os.path.join(HERE, "tasks_testR_koa_261011.txt"))                 # env: rehearsal hook only
TASKS_SHA = os.environ.get("HD_TESTR_TASKS_SHA", os.path.join(HERE, "tasks_testR_koa_261011.sha256"))      # env: rehearsal hook only
LOC = os.environ.get("HD_TESTR_LOC", os.path.join(HS, "experiments_gen3_testR_koa_261011"))     # env: test hook only
OUT = os.path.join(LOC, "testR_verdict_output.txt")
TESTG = os.path.join(HS, "experiments_gen3_gate_koa_261009", "testG")


# ---------------------------------------------------------------- the rule
def dstat(a3, a2):
    """d = m3/m2 - 1, SE(d) (delta method, unpaired), z"""
    a3, a2 = np.asarray(a3, float), np.asarray(a2, float)
    m3, m2 = a3.mean(), a2.mean(); s3, s2 = a3.std(ddof=1) / math.sqrt(len(a3)), a2.std(ddof=1) / math.sqrt(len(a2))
    d = m3 / m2 - 1; se = math.sqrt((s3 / m2) ** 2 + (m3 * s2 / m2 ** 2) ** 2)
    return dict(n3=len(a3), n2=len(a2), m3=m3, m2=m2, d=d, se=se, z=d / se, lo=d - Z95 * se, hi=d + Z95 * se)


def mass_verdict(za, zd, lo_d, hi_d):
    if abs(za) >= Z_CONF or abs(zd) >= Z_CONF: return "DIFFERENCE CONFIRMED"
    if -BAND <= lo_d and hi_d <= BAND: return "NO DIFFERENCE"
    return "NOT DECIDED"


def overall(per_mass, excl_ok):
    if not excl_ok: return "TEST R: NOT DECIDED"
    if all(v == "NO DIFFERENCE" for v in per_mass.values()): return "TEST R: ENGINES AGREE"
    if any(v == "DIFFERENCE CONFIRMED" for v in per_mass.values()): return "TEST R: ENGINES DIFFER"
    return "TEST R: NOT DECIDED"


# ---------------------------------------------------------------- inventory and health
def healthy(e, h):
    """Test G's health rule on runs_health.csv rows (a boolean Series)"""
    c = (h["contact_line"] == 1) & (h["contact_max_px"] <= 1e-6) & (h["health_lines"] == 0)
    if e == "gen2": return c & (h["gen2_policy_line"] == 1) & (h["edmd3_lines"] == 0)
    return c & (h["built_lines"] == 1) & (h["gen3_records"] == 1) & (h["gen3_clean"] == 1) & (h["engine_gen3_line"] == 1)


def registered():
    raw = open(TASKS, "rb").read()
    ok = hashlib.sha256(raw).hexdigest() == open(TASKS_SHA).read().strip()
    want = {}
    for l in raw.decode().splitlines():
        f = l.split(); want.setdefault((f[10], int(f[2])), {})[int(f[3])] = int(f[4])
    return ok, want


def inventory(cells, want, merge_record):
    """per arm: the included runs and the counts; the INVALID reasons"""
    inv, invalid = {}, []
    if not merge_record: invalid.append("no merge record (merged/.merge_complete)")
    for (e, M), c in cells.items():
        W = want.get((e, M), {})
        red, nud, h, fl = c["red"], c["nud"], c["health"], c["failed"]
        if c["build"] != BUILD: invalid.append(f"{e} M={M}: build line {c['build']!r}")
        for name, tab in (("red_nu.csv", red), ("nud_runs.csv", nud), ("runs_health.csv", h)):
            if tab["run"].duplicated().any(): invalid.append(f"{e} M={M}: a run twice in {name}")
            bad = [(int(r), int(s)) for r, s in zip(tab["run"], tab["seed"]) if W.get(int(r)) != int(s)]
            if bad: invalid.append(f"{e} M={M}: {len(bad)} rows of {name} not registered, e.g. {bad[:3]}")
        mm = red.merge(nud, on=["run", "seed"], how="outer", indicator=True)
        # ##CHRIS 2026-10-11 (found by the dry test on a synthetic stand-in, before any Test R data): equal to 1e-14 relative, not bit for
        # bit -- red_nu.csv holds the shortest round-trip decimal, nud_runs.csv 17 digits, and pandas' default CSV parser is not
        # always correctly rounded, so the two texts can read back 1 ulp apart (9.7e-17 seen); the check is that the same values came through
        if (mm["_merge"] != "both").any() or not np.allclose(mm["nu"].to_numpy(float), mm["nu_arg"].to_numpy(float), rtol=1e-14, atol=0.0):
            invalid.append(f"{e} M={M}: nud_runs.csv and red_nu.csv differ")
        hk = h.assign(ok=healthy(e, h))
        good = set(hk.loc[hk["ok"], "run"].astype(int)) & set(red["run"].astype(int))
        failed = set(fl["run"].astype(int)) if len(fl) else set()
        unhealthy = set(red["run"].astype(int)) - good
        absent = set(W) - set(red["run"].astype(int)) - failed
        excl = len(set(W) - good)
        inv[(e, M)] = dict(good=good, n_reg=len(W), failed=len(failed), unhealthy=len(unhealthy), absent=len(absent), excluded=excl,
                           frac=excl / max(len(W), 1), fit_fail=int((nud["fit_ok"] == 0).sum()),
                           nodes=h["node"].nunique() if len(h) else 0)
    return inv, invalid


def load_cell(e, M):
    d = os.path.join(LOC, "merged", "testR", e, CID, f"m_{M}")
    rd = lambda n, **k: pd.read_csv(os.path.join(d, n), **k) if os.path.exists(os.path.join(d, n)) else None
    b = os.path.join(d, ".build_git")
    return dict(dir=d, red=rd("red_nu.csv"), nud=rd("nud_runs.csv"), health=rd("runs_health.csv"),
                failed=(rd("failed_runs.csv") if rd("failed_runs.csv") is not None else pd.DataFrame(columns=["folder", "run"])),
                build=open(b).read().strip() if os.path.exists(b) else None)


# ---------------------------------------------------------------- the design (no data)
def scatter_testG():
    """Test G on KOA: per engine and mass (50, 300) the relative SD of nu_arg and nu_d per trajectory and their correlation; nu_d by
    the registered rule (fit_acf started at the arm's own argmax mean)"""
    import resched_testT_followup_261007 as FU
    out = {}
    for M in (50, 300):
        for e in ENG:
            d = os.path.join(TESTG, e, CID, f"m_{M}")
            r = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run"); z = np.load(os.path.join(d, "acf_runs.npz"))
            nl = min(len(z[f"run{k}"]) for k in r["run"]); ref = float(r["nu"].mean()); dt = float(r["dt"].iloc[0])
            nd = np.array([FU.fit_acf(np.asarray(z[f"run{k}"][:nl], float), dt, ref)[4] / (2 * math.pi) for k in r["run"]])
            out[(e, M)] = dict(sa=r["nu"].std(ddof=1) / ref, sd=nd.std(ddof=1) / ref, rho=float(np.corrcoef(r["nu"], nd)[0, 1]), n=len(r))
    for e in ENG:   # M = 100: in log M between 50 and 300
        w = math.log(100 / 50) / math.log(300 / 50); a, b = out[(e, 50)], out[(e, 300)]
        out[(e, 100)] = dict(sa=math.exp((1 - w) * math.log(a["sa"]) + w * math.log(b["sa"])),
                             sd=math.exp((1 - w) * math.log(a["sd"]) + w * math.log(b["sd"])), rho=(1 - w) * a["rho"] + w * b["rho"], n=0)
    return out


def se_rho(S, M, n=NSEED):
    """SE(d) of both estimators and the correlation of d_arg with d_d, for n trajectories per arm"""
    a2, a3 = S[("gen2", M)], S[("gen3", M)]
    va, vd = (a2["sa"] ** 2 + a3["sa"] ** 2) / n, (a2["sd"] ** 2 + a3["sd"] ** 2) / n
    cov = (a2["rho"] * a2["sa"] * a2["sd"] + a3["rho"] * a3["sa"] * a3["sd"]) / n
    return math.sqrt(va), math.sqrt(vd), cov / math.sqrt(va * vd)


def outcome_probs(S, true, ndraw=2_000_000, seed=20261317):
    """Monte Carlo: per mass d_arg, d_d ~ bivariate normal (mean = the true difference, SEs and correlation from S); the rule"""
    rng = np.random.default_rng(seed); res = {}
    per = {}
    for M in MASSES:
        sa, sd, rho = se_rho(S, M); t = true.get(M, 0.0)
        g1 = rng.standard_normal(ndraw); g2 = rho * g1 + math.sqrt(1 - rho ** 2) * rng.standard_normal(ndraw)
        da, dd = t + sa * g1, t + sd * g2
        conf = (np.abs(da / sa) >= Z_CONF) | (np.abs(dd / sd) >= Z_CONF)
        nodiff = ~conf & (dd - Z95 * sd >= -BAND) & (dd + Z95 * sd <= BAND)
        per[M] = (conf, nodiff)
    agree = np.logical_and.reduce([per[M][1] for M in MASSES]); differ = np.logical_or.reduce([per[M][0] for M in MASSES])
    res["AGREE"], res["DIFFER"] = agree.mean(), differ.mean(); res["NOT DECIDED"] = 1 - res["AGREE"] - res["DIFFER"]
    res["per_mass"] = {M: (per[M][0].mean(), per[M][1].mean()) for M in MASSES}
    return res


def design(S=None):
    from scipy.stats import norm
    if S is None: S = scatter_testG()
    print("# TEST R ON KOA -- the design of the verdict (decision 28; 261012 sec. 4.7.57), printed before any Test R data\n")
    print("## Every threshold\n")
    print(f"- per mass and estimator: d = mean(gen3)/mean(gen2) - 1, SE(d)^2 = (SE3/m2)^2 + (m3 SE2/m2^2)^2, SE = SD/sqrt(n), z = d/SE(d)")
    print(f"- DIFFERENCE CONFIRMED: |z_arg| >= {Z_CONF} or |z_d| >= {Z_CONF}")
    print(f"- NO DIFFERENCE: |z_arg| < {Z_CONF} and |z_d| < {Z_CONF} and d_d +- {Z95:.4f} SE(d_d) within +-{100 * BAND:.2f} % (inclusive)")
    print("- NOT DECIDED: otherwise")
    print("- overall: ENGINES AGREE = NO DIFFERENCE at M = 50, 100 and 300; ENGINES DIFFER = DIFFERENCE CONFIRMED at any mass; NOT DECIDED otherwise")
    print(f"- health: Test G's rule (testG_worker_koa.sh); contact <= 1e-6 px; excluded runs counted, never rerun; more than {100 * EXCL:.0f} % "
          f"excluded in any arm (more than {int(EXCL * NSEED)} of {NSEED}) -> TEST R: NOT DECIDED (overall)")
    print(f"- INVALID (no verdict): no merge record, another build line than '{BUILD}', unregistered or duplicated rows, nu_arg of "
          "nud_runs.csv unequal to red_nu.csv, the task list not registered")
    print(f"- nu_d fit failures: left out of the nu_d mean, counted (information)\n")
    p1 = 2 * norm.sf(Z_CONF)
    print(f"one |z| >= {Z_CONF} by chance: {p1:.2e}\n")
    print("## The per-trajectory scatter of Test G on KOA (400 per arm; nu_d by the registered rule), and the SEs at 3600 per arm\n")
    print("| M | engine | SD nu_arg [% of nu] | SD nu_d [% of nu] | corr(nu_arg, nu_d) | source |\n|---|---|---|---|---|---|")
    for M in MASSES:
        for e in ENG:
            s = S[(e, M)]
            print(f"| {M} | {e} | {100 * s['sa']:.3f} | {100 * s['sd']:.3f} | {s['rho']:+.3f} | "
                  f"{'Test G on KOA, ' + str(s['n']) + ' trajectories' if s['n'] else 'interpolated in log M between 50 and 300'} |")
    print("\n| M | SE(d_arg) [%] | SE(d_d) [%] | corr(d_arg, d_d) | d_d must lie within [%] for NO DIFFERENCE |\n|---|---|---|---|---|")
    for M in MASSES:
        sa, sd, rho = se_rho(S, M)
        print(f"| {M} | {100 * sa:.4f} | {100 * sd:.4f} | {rho:+.3f} | +-{100 * (BAND - Z95 * sd):.4f} |")
    print("\n## Probabilities (Monte Carlo, 2e6 draws per scenario, default_rng(20261317); MC SE <= 0.0004)\n")
    r0 = outcome_probs(S, {})
    print("per mass under no true difference: " + "; ".join(f"M = {M}: P(DIFFERENCE CONFIRMED) {c:.2e}, P(NO DIFFERENCE) {nd:.4f}"
                                                          for M, (c, nd) in r0["per_mass"].items()))
    print(f"\nTHE FALSE-ALARM RATE of the overall verdict (no true difference at any mass): P(ENGINES DIFFER) = {r0['DIFFER']:.2e}; "
          f"P(not ENGINES AGREE) = {1 - r0['AGREE']:.4f}\n")
    print("| true difference at M = 50 (none at 100, 300) | P(ENGINES AGREE) | P(ENGINES DIFFER) | P(NOT DECIDED) |\n|---|---|---|---|")
    for t in (0.0, 0.001, 0.002, 0.003):
        r = outcome_probs(S, {50: t}) if t else r0
        print(f"| {100 * t:.1f} % | {r['AGREE']:.4f} | {r['DIFFER']:.4f} | {r['NOT DECIDED']:.4f} |")
    print("\n(Test G on KOA measured gen3 - gen2 at M = 50: +0.635 % by nu_arg and +0.244 % by nu_d; Report G2.)")
    return S


# ---------------------------------------------------------------- the self-test
def selftest():
    ok = {}; rng = np.random.default_rng(1)
    # 1. SE(d) of the ratio of two means against the spread of d over many synthetic replications
    ds = [dstat(1.002 * (1 + 0.02 * rng.standard_normal(3600)), 1 + 0.025 * rng.standard_normal(3600))["d"] for _ in range(2000)]
    r = dstat(1.002 * (1 + 0.02 * rng.standard_normal(3600)), 1 + 0.025 * rng.standard_normal(3600))
    ok["SE(d) by the delta method = the spread of d over 2000 replications (5 %); mean d = 0.2 %"] = \
        abs(r["se"] / np.std(ds, ddof=1) - 1) < 0.05 and abs(np.mean(ds) - 0.002) < 5e-5
    # 2. the per-mass rule at its edges
    ok["rule: |z_arg| = 3.5 -> DIFFERENCE CONFIRMED; 3.49 with the interval inside -> NO DIFFERENCE"] = \
        mass_verdict(3.5, 0, -1e-4, 1e-4) == "DIFFERENCE CONFIRMED" and mass_verdict(3.49, 0, -1e-4, 1e-4) == "NO DIFFERENCE"
    ok["rule: z_d = -3.6 -> DIFFERENCE CONFIRMED; the interval [-0.15, +0.15] % exactly -> NO DIFFERENCE; +0.1501 % -> NOT DECIDED"] = \
        mass_verdict(0, -3.6, 0, 0) == "DIFFERENCE CONFIRMED" and mass_verdict(1, 1, -0.0015, 0.0015) == "NO DIFFERENCE" and \
        mass_verdict(1, 1, -0.001, 0.001501) == "NOT DECIDED"
    ok["overall: AGREE / DIFFER / NOT DECIDED, and NOT DECIDED when more than 1 % is excluded"] = \
        overall({50: "NO DIFFERENCE", 100: "NO DIFFERENCE", 300: "NO DIFFERENCE"}, True) == "TEST R: ENGINES AGREE" and \
        overall({50: "DIFFERENCE CONFIRMED", 100: "NOT DECIDED", 300: "NO DIFFERENCE"}, True) == "TEST R: ENGINES DIFFER" and \
        overall({50: "NOT DECIDED", 100: "NO DIFFERENCE", 300: "NO DIFFERENCE"}, True) == "TEST R: NOT DECIDED" and \
        overall({50: "NO DIFFERENCE", 100: "NO DIFFERENCE", 300: "NO DIFFERENCE"}, False) == "TEST R: NOT DECIDED"
    # 3. the inventory on synthetic tables
    W = {("gen2", 50): {0: 11, 1: 12, 2: 13}, ("gen3", 50): {0: 11, 1: 12, 2: 13}}
    def cell(e, runs, seeds, bad=(), build=BUILD, failed=()):
        red = pd.DataFrame(dict(run=runs, seed=seeds, nu=[0.07] * len(runs), nu_pred=0.068, dt=0.43, n=6752))
        nud = pd.DataFrame(dict(run=runs, seed=seeds, nu_arg=[0.07] * len(runs), nu_d=[0.071] * len(runs), fit_ok=1))
        h = pd.DataFrame(dict(run=runs, seed=seeds, node="cn", seconds=3, built_lines=int(e == "gen3"), gen3_records=int(e == "gen3"),
                              gen3_clean=int(e == "gen3"), engine_gen3_line=int(e == "gen3"), gen2_policy_line=int(e == "gen2"),
                              edmd3_lines=(2 if e == "gen3" else 0), health_lines=[1 if r in bad else 0 for r in runs], contact_line=1,
                              contact_max_px=3e-11))
        return dict(red=red, nud=nud, health=h, failed=pd.DataFrame(dict(folder=[f".failed_run{r}_x_chunk1" for r in failed], run=list(failed))),
                    build=build)
    cells = {("gen2", 50): cell("gen2", [0, 1], [11, 12], failed=(2,)), ("gen3", 50): cell("gen3", [0, 1, 2], [11, 12, 13], bad=(1,))}
    inv, invalid = inventory(cells, W, True)
    ok["inventory: a refused run and an unhealthy run are excluded and counted; no INVALID"] = not invalid and \
        inv[("gen2", 50)]["excluded"] == 1 and inv[("gen2", 50)]["failed"] == 1 and inv[("gen3", 50)]["excluded"] == 1 and \
        inv[("gen3", 50)]["unhealthy"] == 1 and inv[("gen3", 50)]["good"] == {0, 2}
    cells2 = {("gen2", 50): cell("gen2", [0, 1, 5], [11, 12, 99]), ("gen3", 50): cell("gen3", [0, 1, 2], [11, 12, 13], build="00ALLINONE  git 0000000  target koa")}
    inv2, invalid2 = inventory(cells2, W, False)
    c1u = cell("gen2", [0, 1, 2], [11, 12, 13]); c1u["nud"]["nu_arg"] = np.nextafter(c1u["red"]["nu"].to_numpy(float), 1.0)
    c1e = cell("gen2", [0, 1, 2], [11, 12, 13]); c1e["nud"]["nu_arg"] = c1e["red"]["nu"].to_numpy(float) * (1 + 1e-12)
    ok["inventory: nu_arg 1 ulp from red_nu.csv's is the same value; 1e-12 relative is INVALID"] = \
        not inventory({("gen2", 50): c1u}, {("gen2", 50): W[("gen2", 50)]}, True)[1] and \
        any("differ" in x for x in inventory({("gen2", 50): c1e}, {("gen2", 50): W[("gen2", 50)]}, True)[1])
    ok["inventory: an unregistered row, another build line and a missing merge record -> INVALID (all three named)"] = \
        len(invalid2) >= 3 and any("merge record" in x for x in invalid2) and any("build line" in x for x in invalid2) and \
        any("not registered" in x for x in invalid2)
    # 4. the design's Monte Carlo against closed forms
    from scipy.stats import norm
    S0 = {(e, M): dict(sa=0.026, sd=0.015, rho=0.0, n=400) for e in ENG for M in MASSES}
    r0 = outcome_probs(S0, {}, ndraw=400_000)
    p1 = 2 * norm.sf(Z_CONF); exact = 1 - (1 - p1) ** 6
    sd = se_rho(S0, 50)[1]; w = BAND - Z95 * sd; pin = (2 * norm.cdf(w / sd) - 1)
    ok["design MC: P(ENGINES DIFFER) with independent estimators = 1 - (1 - 2 Phi(-3.5))^6 (within MC error)"] = abs(r0["DIFFER"] - exact) < 4 * math.sqrt(exact / 400_000) + 1e-6
    ok["design MC: P(NO DIFFERENCE) at one mass = P(|d_d| <= 0.15 % - 1.96 SE) (within 0.003)"] = abs(r0["per_mass"][50][1] - pin) < 0.003
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else '**FAIL**'} |")
    print(f"\nSELF-TEST: {'PASS' if all(ok.values()) else 'FAIL'} ({sum(ok.values())} of {len(ok)})")
    return all(ok.values())


# ---------------------------------------------------------------- the one printing
def information(cells, inv, S):
    import resched_testT_followup_261007 as FU
    import testG_koa_posthoc_261011 as PH
    if min(len(v["good"]) for v in inv.values()) < 20:
        print("\n(information rows need at least 20 trajectories per arm: skipped)"); return
    print("\n## Information (no verdict): tau_r, Q, nu_0 and delta_M per engine (Report G2's functions; jackknife groups of 10; start: "
          "the arm's own argmax mean; delta_M against c_300)\n")
    print("| engine | M | n | mean nu_arg | mean nu_d | tau_r | Q | nu_0 |\n|---|---|---|---|---|---|---|---|")
    import tests_20260913 as T
    St = {}
    for (e, M), c in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        g = sorted(inv[(e, M)]["good"]); red = c["red"].set_index("run").loc[g]; nud = c["nud"].set_index("run").loc[g]
        z = np.load(os.path.join(c["dir"], "acf_runs.npz")); nl = min(len(z[f"run{r}"]) for r in g)
        acf = np.array([z[f"run{r}"][:nl] for r in g], dtype=float)
        ndv = nud["nu_d"].to_numpy(float); ndv[nud["fit_ok"].to_numpy() == 0] = np.nan
        f, s, _ = PH.block_stats(red["nu"].to_numpy(float), ndv, acf, float(red["dt"].iloc[0]), float(red["nu"].mean()))
        St[(e, M)] = (f, s)
        print(f"| {e} | {M} | {len(g)} | {f[0]:.7f} +- {s[0]:.7f} | {f[1]:.7f} +- {s[1]:.7f} | {f[2]:.1f} +- {s[2]:.1f} | {f[3]:.2f} +- "
              f"{s[3]:.2f} | {f[4]:.7f} +- {s[4]:.7f} |")
    x = {M: T.x_of(M, 10.0) for M in MASSES}
    print("\n| estimator | engine | delta_50 [%] | delta_100 [%] |\n|---|---|---|---|")
    for k, lab in ((0, "nu_arg"), (1, "nu_d"), (4, "nu_0")):
        for e in ENG:
            c = {M: St[(e, M)][0][k] / x[M] for M in MASSES}; s = {M: St[(e, M)][1][k] / x[M] for M in MASSES}
            d, _ = PH.deltas(c, s, (300,))
            print(f"| {lab} | {e} | {100 * d[50][0]:+.3f} +- {100 * d[50][1]:.3f} | {100 * d[100][0]:+.3f} +- {100 * d[100][1]:.3f} |")
    print("\n## Information: the argmax histograms (bin: gen2/gen3)\n")
    for M in MASSES:
        r2 = cells[("gen2", M)]["red"].set_index("run").loc[sorted(inv[("gen2", M)]["good"])]
        r3 = cells[("gen3", M)]["red"].set_index("run").loc[sorted(inv[("gen3", M)]["good"])]
        df = 1 / (float(r2["n"].iloc[0]) * float(r2["dt"].iloc[0])); h2, h3 = PH.hist_bins(r2["nu"], df), PH.hist_bins(r3["nu"], df)
        print(f"M = {M} (df = {100 * df / r2['nu'].mean():.3f} % of nu): " + ", ".join(f"{b}: {h2.get(b, 0)}/{h3.get(b, 0)}" for b in sorted(set(h2) | set(h3))))
    print("\n## Information: d from Test G and Test R pooled (Test G's nu_d by the same rule from its fetched ACFs)\n")
    print("| M | estimator | Test R d [%] | Test G d [%] | pooled d [%] | pooled z |\n|---|---|---|---|---|---|")
    for M in (50, 300):
        G = {}
        for e in ENG:
            d = os.path.join(TESTG, e, CID, f"m_{M}"); r = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run")
            z = np.load(os.path.join(d, "acf_runs.npz")); nl = min(len(z[f"run{k}"]) for k in r["run"]); ref = float(r["nu"].mean())
            G[e] = (r["nu"].to_numpy(float), np.array([FU.fit_acf(np.asarray(z[f"run{k}"][:nl], float), float(r["dt"].iloc[0]), ref)[4] / (2 * math.pi) for k in r["run"]]))
        for i, lab in ((0, "nu_arg"), (1, "nu_d")):
            R = {e: (cells[(e, M)]["nud"].set_index("run").loc[sorted(inv[(e, M)]["good"])]["nu_arg" if i == 0 else "nu_d"]).to_numpy(float) for e in ENG}
            R = {e: v[np.isfinite(v)] for e, v in R.items()}
            a, b, p = dstat(R["gen3"], R["gen2"]), dstat(G["gen3"][i], G["gen2"][i]), dstat(np.r_[R["gen3"], G["gen3"][i]], np.r_[R["gen2"], G["gen2"][i]])
            print(f"| {M} | {lab} | {100 * a['d']:+.3f} +- {100 * a['se']:.3f} | {100 * b['d']:+.3f} +- {100 * b['se']:.3f} | {100 * p['d']:+.3f} +- {100 * p['se']:.3f} | {p['z']:+.2f} |")
    print("\n## Information: KOA's nu_d against a Mac recomputation (the first 20 trajectories of every arm)\n")
    for (e, M), c in sorted(cells.items()):
        nud = c["nud"].sort_values("run").head(20); z = np.load(os.path.join(c["dir"], "acf_runs.npz"))
        nl = min(len(z[f"run{r}"]) for r in c["red"]["run"]); ref = float(nud["nu_ref"].iloc[0]); dt = float(nud["dt"].iloc[0])
        mac = np.array([FU.fit_acf(np.asarray(z[f"run{r}"][:nl], float), dt, ref)[4] / (2 * math.pi) for r in nud["run"]])
        print(f"{e} M = {M}: max |KOA / Mac - 1| = {np.nanmax(np.abs(nud['nu_d'].to_numpy(float) / mac - 1)):.2e}")


def analysis():
    if os.path.exists(OUT): sys.exit(f"STOP: {OUT} exists: the verdict has been printed (one printing)")
    buf = io.StringIO(); real = sys.stdout

    class Tee:
        def write(self, s): real.write(s); buf.write(s)
        def flush(self): real.flush()
    sys.stdout = Tee()
    try:
        print("# TEST R ON KOA -- the verdict (decision 28; 261012 sec. 4.7.57); the one printing\n")
        reg_ok, want = registered()
        mr = os.path.join(LOC, "merged", ".merge_complete")
        print(f"the task list: {'the registered one' if reg_ok else '**NOT the registered one**'}; the merge record: "
              f"{open(mr).read().strip() if os.path.exists(mr) else '**MISSING**'}\n")
        cells = {(e, M): load_cell(e, M) for e in ENG for M in MASSES}
        miss = [f"{e} M={M}" for (e, M), c in cells.items() if any(c[k] is None for k in ("red", "nud", "health"))]
        if miss or not reg_ok:
            print(f"TEST R: INVALID -- {'missing tables: ' + ', '.join(miss) if miss else 'the task list is not the registered one'}; no verdict"); return
        inv, invalid = inventory(cells, want, os.path.exists(mr))
        print("## Inventory and health\n")
        print("| engine | M | registered | included | refused (failed) | unhealthy in the tables | absent | excluded [%] | nu_d fit failures | nodes |\n"
              "|---|---|---|---|---|---|---|---|---|---|")
        for (e, M), v in sorted(inv.items()):
            print(f"| {e} | {M} | {v['n_reg']} | {len(v['good'])} | {v['failed']} | {v['unhealthy']} | {v['absent']} | {100 * v['frac']:.2f} | "
                  f"{v['fit_fail']} | {v['nodes']} |")
        if invalid:
            print("\nINVALID:\n" + "\n".join(f"  - {x}" for x in invalid)); print("\nTEST R: INVALID -- no verdict"); return
        excl_ok = all(v["frac"] <= EXCL for v in inv.values())
        print(f"\nexcluded at most {100 * EXCL:.0f} % in every arm: {'yes' if excl_ok else '**NO** -> TEST R: NOT DECIDED'}\n")
        print("## The registered numbers\n")
        print("| M | estimator | n gen3 | n gen2 | mean gen3 | mean gen2 | d [%] | SE(d) [%] | 95 % interval [%] | z |\n|---|---|---|---|---|---|---|---|---|---|")
        per = {}
        for M in MASSES:
            R = {}
            for lab, col in (("nu_arg", "nu_arg"), ("nu_d", "nu_d")):
                v = {e: cells[(e, M)]["nud"].set_index("run").loc[sorted(inv[(e, M)]["good"])] for e in ENG}
                a = {e: (v[e][col][v[e]["fit_ok"] == 1] if col == "nu_d" else v[e][col]).to_numpy(float) for e in ENG}
                R[lab] = s = dstat(a["gen3"], a["gen2"])
                print(f"| {M} | {lab} | {s['n3']} | {s['n2']} | {s['m3']:.7f} | {s['m2']:.7f} | {100 * s['d']:+.4f} | {100 * s['se']:.4f} | "
                      f"[{100 * s['lo']:+.4f}, {100 * s['hi']:+.4f}] | {s['z']:+.2f} |")
            per[M] = mass_verdict(R["nu_arg"]["z"], R["nu_d"]["z"], R["nu_d"]["lo"], R["nu_d"]["hi"])
        print("\n" + "; ".join(f"M = {M}: {v}" for M, v in per.items()))
        print(overall(per, excl_ok))
        information(cells, inv, None)
    except Exception:
        import traceback; traceback.print_exc(file=sys.stdout); print("STOP: the script failed (the traceback above is part of the printing)")
    finally:
        sys.stdout = real
        os.makedirs(LOC, exist_ok=True)
        with open(OUT, "w") as fh: fh.write(buf.getvalue())


if __name__ == "__main__":
    if "--selftest" in sys.argv: sys.exit(0 if selftest() else 1)
    if "--design" in sys.argv: design(); sys.exit(0)
    analysis()
