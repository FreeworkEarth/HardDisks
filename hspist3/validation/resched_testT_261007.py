#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10 E3, gate version 2): the analysis of TEST T, written and committed before any of its data.
Test T: cell epi8_H_H10_L10, method B, the campaign protocol, 100 fresh seeds per mass, the same seeds on both policies (minimal
and --legacy-resched) of ONE binary (cluster/resched_gate_261005/tasks_T_epi8_H_H10_L10.txt; seed list SHA-256 printed by
validation/resched_testT_design_261007.py). Read from hspist3/experiments_resched_gate2_261007/ (cluster/resched_gate_261005/
fetch_resched2.sh).
Eleven numbers, minimal minus legacy, registered estimators (paper1_populate_cs_err_20261002.slope_with_errors for c_s with its
unscaled error; the inverse-variance mean of (M + 2N_s/3) omega^2/2 over alpha >= 5 for k_S^dyn, as
paper1_confinement_results_261004.identity; the mean nu of each mass with SD/sqrt(n)), z = difference / sqrt(SE_min^2 + SE_leg^2).
RULE (plan author 2026-10-07, fixed before the data; no extension, no "in between"): PASS if all eleven |z| < z* (two-sided
Bonferroni, family-wise false-fail 5 %, n = 11) AND the permutation p of the nine-mass chi2 >= 0.01 (1e5 relabelings within each
mass pool, numpy default_rng(20261008)); otherwise FAIL.
Also printed: the inventory (every trajectory of the task list present, by run and seed; health lines; the contact audit of every
trajectory; the policy line of every trajectory; the node of every trajectory; one build); the 95 % interval of each relative
difference (the bound on any bias of the build); the z each hypothesis of the replay would give at its observed size; for
information (i) the eleven numbers pooled with the 25 campaign (legacy, 279282b) and 25 replay (minimal, 73fc07f) seeds,
(ii) the null calibration: the FIRST 100 seeds of every mass (seed-list order) in four blocks of 25 per policy, the nine-mass chi2
of the six block pairs within each policy (fixed 2026-10-07 before the data: M = 50 has 400 seeds, the others 100); (iii) the
plain-fluid baseline of the second decision of 2026-10-07, item 3: the per-mass implied sound speed c_s,M = nu_M / x_M with SE
for both policies, and the chi2 of the single-c_s model at the registered unweighted slope and at the weighted slope (8 dof).
usage (from hspist3/):  python3 validation/resched_testT_261007.py      [--dry-run-replay: minimal := the 73fc07f replay,
                        legacy := the 279282b campaign, 25 seeds each, no task-list check -- a test of this script only]
"""
import contextlib, glob, io, math, os, re, sys
from collections import Counter
import numpy as np
from scipy.stats import norm, chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import tests_20260913 as T
from paper1_populate_cs_err_20261002 import slope_with_errors

DRY = "--dry-run-replay" in sys.argv
CID = "epi8_H_H10_L10"
LOC = os.environ.get("HD_RESCHED2_LOC", os.path.join(HS, "experiments_resched_gate2_261007"))   # env: test hook only
REL_T = "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/resched_testT_261007"
TASKS = os.path.join(HS, "cluster", "resched_gate_261005", f"tasks_T_{CID}.txt")
NPERM, RNG_SEED, FW, NNUM, P_CHI = 100000, 20261008, 0.05, 11, 0.01
HYP = {"k_S^dyn": 0.0040, "nu M=500": 0.0062, "nu M=1500": 0.0057, "nu M=50": -0.0105}   # replay's observed sizes (sec. 4.4.9)
CONTACT = re.compile(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), "
                     r"outer walls (\S+), divider (\S+), pistons (\S+)\n")
NODE = re.compile(r"run (\d+) seed (\d+) node (\S+) \(")


def cell_dir(pol):
    if DRY:
        return os.path.join(HS, "experiments_resched_gate_261005" if pol == "minimal" else "", R.REL_B, CID).replace("//", "/")
    return os.path.join(LOC, REL_T, pol, CID)


def load(pol, Ms):
    import pandas as pd
    out = {}
    for M in Ms:
        d = os.path.join(cell_dir(pol), f"m_{M}")
        out[M] = dict(red=pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True), log=open(os.path.join(d, "run.log"), errors="ignore").read(),
                      failed=[open(f, errors="ignore").read() for f in glob.glob(os.path.join(d, ".failed_run*", "stdout.log"))],
                      build=(open(os.path.join(d, ".build_git")).read().strip() if os.path.exists(os.path.join(d, ".build_git")) else None))
    return out


def estimators(means, ses, Ms, x, Ns):
    """Registered estimators, vectorised: c_s (unweighted through-origin slope) and its unscaled error; k_S^dyn and its SE."""
    cs = (means * x).sum(-1) / (x * x).sum(); cs_err = np.sqrt((x * x * ses * ses).sum(-1)) / (x * x).sum()
    heavy = np.array([any(abs(M / (2.0 * Ns) - a) < 1e-9 for a in R.HEAVY) for M in Ms])
    Mh = np.array([M + 2.0 * Ns / 3.0 for M in Ms]); om = 2 * math.pi * means; k = Mh * om * om / 2.0
    s = k * 2 * ses / means; w = np.where(heavy, 1.0 / s ** 2, 0.0)
    return cs, cs_err, (w * k).sum(-1) / w.sum(-1), 1.0 / np.sqrt(w.sum(-1))


def eleven(nu_min, nu_leg, Ms, x, Ns):
    rows = []
    mm = np.array([nu_min[M].mean() for M in Ms]); sm = np.array([nu_min[M].std(ddof=1) / math.sqrt(len(nu_min[M])) for M in Ms])
    ml = np.array([nu_leg[M].mean() for M in Ms]); sl = np.array([nu_leg[M].std(ddof=1) / math.sqrt(len(nu_leg[M])) for M in Ms])
    cm, cme, km, kme = estimators(mm, sm, Ms, x, Ns); cl, cle, kl, kle = estimators(ml, sl, Ms, x, Ns)
    rows.append(("k_S^dyn", km, kme, kl, kle)); rows.append(("c_s", cm, cme, cl, cle))
    for j, M in enumerate(Ms): rows.append((f"nu M={M}", mm[j], sm[j], ml[j], sl[j]))
    return rows


def main():
    zs = norm.isf(FW / (2 * NNUM))
    with contextlib.redirect_stdout(io.StringIO()):
        c = {x["cid"]: x for x in R.cells()}[CID]; co = dict(c); R.method_B(co)
    Ms, Ns, x = c["Ms"], c["Ns"], np.asarray(co["x"], float)
    from paper1_populate_cs_err_20261002 import slope_with_errors as swe   # check the vectorised c_s against the registered one
    chk = swe(x, np.array([r["nu"] for r in co["B"]]), np.array([r["se"] for r in co["B"]]))
    cs_v = estimators(np.array([r["nu"] for r in co["B"]]), np.array([r["se"] for r in co["B"]]), Ms, x, Ns)
    if abs(cs_v[0] - chk[0]) > 1e-12 or abs(cs_v[1] - chk[1]) > 1e-12: sys.exit("STOP: vectorised c_s differs from slope_with_errors")
    print(f"# Test T (gate version 2, 261012 sec. 4.4.10 E3){'  [DRY RUN on the replay: a test of this script, no verdict]' if DRY else ''}\n")
    print(f"z* = {zs:.4f} (two-sided Bonferroni, family-wise false-fail {FW:.0%}, n = {NNUM}); permutation chi2 criterion p >= {P_CHI}")
    D = {pol: load(pol, Ms) for pol in ("minimal", "legacy")}
    # ---------------------------------------------------------------- inventory
    print("\n## Inventory\n")
    want = {}
    if not DRY:
        for l in open(TASKS):
            f = l.split(); want.setdefault((f[10], int(f[2])), set()).add((int(f[3]), int(f[4])))
    print("| policy | M | trajectories matched (expected) | n finite > 0 | log sections | health lines (failed runs incl.) | policy line wrong | "
          "max contact gap [px] | contact lines missing | nodes so far | ok |\n|---|---|---|---|---|---|---|---|---|---|---|")
    inv_ok = True; builds = set(); nodes_by = {"minimal": Counter(), "legacy": Counter()}
    for pol in ("minimal", "legacy"):
        for M in Ms:
            e = D[pol][M]; red = e["red"]; have = {(int(a), int(b)) for a, b in zip(red["run"], red["seed"])}
            exp = want.get((pol, M), have) if not DRY else have
            matched = len(exp & have) if len(have) == len(red) else 0
            nfin = int((np.isfinite(red["n"].astype(float)) & (red["n"].astype(float) > 0)).sum())
            secs = e["log"].split("##RUN")[1:]; nh = sum(s.count("[EDMD-HEALTH]") for s in secs) + sum(f.count("[EDMD-HEALTH]") for f in e["failed"])
            want_line = "minimal (default)" if pol == "minimal" else "legacy (full reschedule"
            bad_pol = 0 if DRY else sum(f"[EDMD-RESCHED] divider events: {want_line}" not in s for s in secs)
            cg, nmiss = 0.0, 0
            for s in secs:
                m = CONTACT.search(s)
                if not m: nmiss += 1; continue
                cg = max(cg, max(float(v) for v in m.groups()[1:]))
                n = NODE.search(s.split("\n", 1)[0])
                if n: nodes_by[pol][n.group(3)] += 1
            if e["build"]: builds.add(e["build"])
            ok = matched == len(exp) and nfin == len(red) and len(secs) == len(red) and nh == 0 and bad_pol == 0 and (DRY or (nmiss == 0 and cg <= 1e-6))
            inv_ok &= ok
            print(f"| {pol} | {M} | {matched} ({len(exp)}) | {nfin} | {len(secs)} | {nh} | {bad_pol} | {cg:.2e} | {nmiss} | "
                  f"{len(nodes_by[pol])} | {'yes' if ok else '**NO**'} |")
    shared = set(nodes_by["minimal"]) & set(nodes_by["legacy"])
    print(f"\nnodes: minimal {dict(nodes_by['minimal'])}; legacy {dict(nodes_by['legacy'])}; shared by both policies: {sorted(shared)}")
    print(f"builds (.build_git): {sorted(builds)}; inventory {'clean' if inv_ok and len(builds) <= 1 else '**NOT CLEAN**'}")
    nu = {pol: {M: D[pol][M]["red"]["nu"].to_numpy(float) for M in Ms} for pol in D}
    # ---------------------------------------------------------------- the eleven numbers
    print("\n## The eleven numbers, minimal minus legacy\n")
    print("| number | minimal | SE | legacy | SE | difference | relative [%] | 95 % interval of the relative difference [%] | z | abs(z) < z* |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    rows = eleven(nu["minimal"], nu["legacy"], Ms, x, Ns); zmap = {}
    for name, a, sa, b, sb in rows:
        d = a - b; s = math.hypot(sa, sb); z = d / s if s > 0 and math.isfinite(s) else float("nan"); zmap[name] = (z, s / abs(b))
        print(f"| {name} | {a:.6g} | {sa:.3g} | {b:.6g} | {sb:.3g} | {d:+.3g} | {100 * d / b:+.3f} | [{100 * (d - 1.96 * s) / b:+.3f}, "
              f"{100 * (d + 1.96 * s) / b:+.3f}] | {z:+.2f} | {'yes' if math.isfinite(z) and abs(z) < zs else '**NO**'} |")
    # ---------------------------------------------------------------- permutation chi2
    rng = np.random.default_rng(RNG_SEED)
    zobs = np.array([zmap[f"nu M={M}"][0] for M in Ms]); chi_obs = float((zobs ** 2).sum())
    pchi = np.zeros(NPERM)
    for M in Ms:
        a, b = nu["minimal"][M], nu["legacy"][M]; pool = np.r_[a, b]; na = len(a)
        v = pool[np.argsort(rng.random((NPERM, len(pool))), axis=1)]; pa, pb = v[:, :na], v[:, na:]
        pz = (pa.mean(1) - pb.mean(1)) / np.hypot(pa.std(1, ddof=1) / math.sqrt(na), pb.std(1, ddof=1) / math.sqrt(len(pool) - na))
        pchi += pz * pz
    p_perm = float((pchi >= chi_obs).mean())
    print(f"\nnine-mass chi2 = {chi_obs:.2f} (nominal p {CHI2.sf(chi_obs, 9):.4f}); permutation p = {p_perm:.4f} "
          f"({NPERM} relabelings within each mass pool, default_rng({RNG_SEED}))")
    allz = all(math.isfinite(zmap[r[0]][0]) and abs(zmap[r[0]][0]) < zs for r in rows)
    verdict = allz and p_perm >= P_CHI
    print(f"\nTEST T: {'(dry run, no verdict) ' if DRY else ''}{'PASS' if verdict else 'FAIL'} -- all eleven abs(z) < {zs:.4f}: "
          f"{'yes' if allz else 'NO'}; permutation p of the nine-mass chi2 >= {P_CHI}: {'yes' if p_perm >= P_CHI else 'NO'}; "
          f"inventory {'clean' if inv_ok else 'NOT CLEAN (verdict void)'}")
    print("\n### The replay's hypotheses: the z each would give here at its observed size, and the observed z\n")
    print("| hypothesis | observed size | expected z if real | observed z |\n|---|---|---|---|")
    for name, size in HYP.items():
        z, rel_sd = zmap[name]; print(f"| {name} | {100 * size:+.2f} % | {size / rel_sd:+.2f} | {z:+.2f} |")
    # ---------------------------------------------------------------- information
    if not DRY:
        print("\n## Information (i): pooled with the 25 campaign (legacy, 279282b) and 25 replay (minimal, 73fc07f) seeds\n")
        import pandas as pd
        RB = os.path.join(R.REL_B, CID)
        old = {M: pd.read_csv(os.path.join(HS, RB, f"m_{M}", "red_nu.csv"))["nu"].to_numpy(float) for M in Ms}
        rep = {M: pd.read_csv(os.path.join(HS, "experiments_resched_gate_261005", RB, f"m_{M}", "red_nu.csv"))["nu"].to_numpy(float) for M in Ms}
        print("| number | difference | relative [%] | z |\n|---|---|---|---|")
        for name, a, sa, b, sb in eleven({M: np.r_[nu["minimal"][M], rep[M]] for M in Ms}, {M: np.r_[nu["legacy"][M], old[M]] for M in Ms}, Ms, x, Ns):
            print(f"| {name} | {a - b:+.3g} | {100 * (a - b) / b:+.3f} | {(a - b) / math.hypot(sa, sb):+.2f} |")
    print("\n## Information (ii): null calibration -- nine-mass chi2 between blocks of 25 of the first 100 seeds of every mass, per policy\n")
    print("| policy | block pair | chi2 (9 dof) | nominal p |\n|---|---|---|---|")
    for pol in ("minimal", "legacy"):
        nb = min(len(nu[pol][M]) for M in Ms) // 25; nb = min(nb, 4)   # the first 100 seeds of every mass: four blocks
        for i in range(nb):
            for j in range(i + 1, nb):
                z = []
                for M in Ms:
                    a, b = nu[pol][M][25 * i:25 * i + 25], nu[pol][M][25 * j:25 * j + 25]
                    z.append((a.mean() - b.mean()) / math.hypot(a.std(ddof=1) / 5, b.std(ddof=1) / 5))
                c2 = float(np.sum(np.square(z))); print(f"| {pol} | {i + 1}-{j + 1} | {c2:.2f} | {CHI2.sf(c2, 9):.3f} |")
        if nb < 2: print(f"| {pol} | (only {min(len(nu[pol][M]) for M in Ms)} seeds in some mass: no block pair) | | |")
    print("\n## Information (iii): plain-fluid baseline -- per-mass implied sound speed and the single-c_s fit, per policy\n")
    print("| policy | M | c_s,M = nu/x | SE |\n|---|---|---|---|")
    for pol in ("minimal", "legacy"):
        y = np.array([nu[pol][M].mean() for M in Ms]); sy = np.array([nu[pol][M].std(ddof=1) / math.sqrt(len(nu[pol][M])) for M in Ms])
        for M, xm, ym, sm in zip(Ms, x, y, sy): print(f"| {pol} | {M} | {ym / xm:.5f} | {sm / xm:.5f} |")
    print("\n| policy | unweighted c_s (registered) | chi2 at it (8 dof) | weighted c_s +- SE | chi2 at it (8 dof) | p (weighted) |\n|---|---|---|---|---|---|")
    for pol in ("minimal", "legacy"):
        y = np.array([nu[pol][M].mean() for M in Ms]); sy = np.array([nu[pol][M].std(ddof=1) / math.sqrt(len(nu[pol][M])) for M in Ms])
        su = float((x * y).sum() / (x * x).sum()); w = 1.0 / sy ** 2; sw = float((w * x * y).sum() / (w * x * x).sum())
        cu = float((((y - su * x) / sy) ** 2).sum()); cw = float((((y - sw * x) / sy) ** 2).sum())
        print(f"| {pol} | {su:.5f} | {cu:.1f} | {sw:.5f} +- {1 / math.sqrt(float((w * x * x).sum())):.5f} | {cw:.1f} | {CHI2.sf(cw, 8):.3f} |")


if __name__ == "__main__":
    main()
