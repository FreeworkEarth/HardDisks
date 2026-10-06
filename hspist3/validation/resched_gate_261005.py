#!/usr/bin/env python3
"""##CHRIS 2026-10-05: the engine gate of the minimal divider rescheduling (branch engine-divider-resched, 261012 sec. 4.4),
Mac side, after hspist3/cluster/resched_gate_261005/fetch_resched.sh. Prints G-E2 (the comparison of the KOA G-E2 job, the same
code as on KOA), G-E3 (the three replayed cells clean), G-E4 (the nine numbers against the 279282b campaign) and G-E5 (profile).

Replayed cells (the 279282b campaign's own task lines: same masses, positions, seeds, protocol), read from the separate tree
hspist3/experiments_resched_gate_261005/ (never mixed with the canonical one; the registered estimators are called with the
explicit flag R.ALLOW_NEW_BUILD = True):
  B  e0p10_H_H10_L39.25  (eta 0.10 anchor)   B  epi8_H_H10_L10  (pi/8 anchor)   AF  epi8_H_H10_L10  (pi/8 anchor, held divider)

Registered before the replay data exist (261012 sec. 4.4):
  G-E3  every trajectory present (B 9 x 25, AF 115); no [EDMD-HEALTH] line (overlap_repairs, wall_overdue, forced_advance,
        clamps, past_events -- the last is the count of events popped with a time before the current time, i.e. negative
        collision times); every log carries "[EDMD-RESCHED] divider events: minimal" and every A-fixed summary the replay
        build; divider ledgers: B -- every trajectory's |dE/E| (E = gas + divider KE, [EDMD-ENERGY]) <= 1e-10 in the hold
        and in the record; AF -- u_wall_max = 0 and W_div = 0 in every seed.
  G-E4  nine numbers, three per cell, each with the registered estimator on the 279282b data (old) and on the replay (new):
          B cells:  c_s (through-origin slope, statistical seed error c_s_err), k_S^dyn (heavy masses, C1),
                    Gamma at alpha = 5 (energy-decay rate 2/tau_r of the divider mode, jackknife error);
          AF cell:  F(L_0) = (F_L + F_R)/2 at x = 0, k_T (five-point stencil), k_static = k_T + F^2/(N_s kT).
        z = (new - old) / sqrt(SE_old^2 + SE_new^2): with the same seeds but a different build the trajectories decorrelate
        within the 200 sigma-time equilibration, so old and new are independent realisations. PASS if all nine |z| < 2 (the
        plan author's rule). Also printed: (new - old)/SE_old and the false-fail probability of the rule under the null.
  The old values are also checked against the recorded CSVs (261004_p1_confinement_cells.csv, 261005_p1_identity_afix_cells.csv).
usage (from hspist3/):  python3 validation/resched_gate_261005.py            [--dry-run-old: "new" := the old data, plumbing test]
"""
import contextlib, glob, io, math, os, re, sys
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS); sys.path.insert(0, os.path.join(HS, "cluster", "resched_gate_261005"))
import paper1_confinement_results_261004 as R
import paper1_confinement_afix_261005 as AF
import tests_20260913 as T

DRY = "--dry-run-old" in sys.argv
LOC = HS if DRY else os.path.join(HS, "experiments_resched_gate_261005")
CELLS = (("B", "e0p10_H_H10_L39.25"), ("B", "epi8_H_H10_L10"), ("AF", "epi8_H_H10_L10"))
E_TOL = 1e-10
EN = re.compile(r"\[EDMD-ENERGY\] (\d) .*?E_tot=(\S+) resched=(\w+)")
OLD_PROFILE = {"held 100 legacy": 2.3, "held 400 legacy": 56.4}   # 279282b, KOA job 14983181 (261012 sec. 4.3) [DATA]


@contextlib.contextmanager
def root(mod, path):
    old = mod.HS; mod.HS = path
    try: yield
    finally: mod.HS = old


def b_values(c):
    """c_s, k_S^dyn, Gamma(alpha = 5) with errors, for the current R.HS."""
    R.method_B(c); R.identity(c); R.damping(c)
    g = [d for d in c["damp"] if abs(d["alpha"] - 5.0) < 1e-9][0]
    return dict(cs=(c["cs"], c["cs_err"]), kS=(c["kS"], c["s_kS"]), G5=(g["Gamma"], g["Gamma"] * g["s_tau_r"] / g["tau_r"]))


def af_values(c):
    a = AF.inventory_and_static(c)
    if a.get("error"): sys.exit(f"AF {c['cid']}: {a['error']}")
    s0 = a["st"][0]; sF = 0.5 * math.hypot(s0["sFL"], s0["sFR"])
    s_st = math.hypot(a["s_kT"], 2 * a["F0"] * sF / (c["Ns"] * a["temp"]))
    return dict(F0=(a["F0"], sF), kT=(a["kT"], a["s_kT"]), static=(a["static"], s_st)), a


def g_e3(CS, builds_af):
    print("### G-E3 -- the replayed cells are clean\n")
    print("| cell | method | trajectories (expected) | health lines | logs without '[EDMD-RESCHED] ... minimal' | build(s) | "
          "max abs(dE/E) hold | max abs(dE/E) record | u_wall_max = 0 and W_div = 0 | PASS |\n|---|---|---|---|---|---|---|---|---|---|")
    allok = True
    for mode, cid in CELLS:
        c = CS[cid]
        if mode == "B":
            d0 = os.path.join(LOC, R.REL_B, cid); n = 0; nh = 0; nores = 0; eh = er = 0.0; bad_e = 0
            for M in c["Ms"]:
                n += len(pd.read_csv(os.path.join(d0, f"m_{M}", "red_nu.csv")))
                log = open(os.path.join(d0, f"m_{M}", "run.log"), errors="ignore").read()
                for sec in log.split("##RUN")[1:]:
                    nh += sec.count("[EDMD-HEALTH]")
                    nores += "[EDMD-RESCHED] divider events: minimal" not in sec
                    e = {int(m.group(1)): float(m.group(2)) for m in EN.finditer(sec)}
                    if set(e) != {0, 1, 2}: bad_e += 1; continue
                    eh = max(eh, abs(e[1] / e[0] - 1)); er = max(er, abs(e[2] / e[1] - 1))
            exp = 25 * len(c["Ms"]); bld = "(B: none recorded)"
            ok = n == exp and nh == 0 and nores == 0 and bad_e == 0 and eh <= E_TOL and er <= E_TOL
            led = f"n/a (energy lines missing in {bad_e})" if bad_e else "n/a"
        else:
            d0 = os.path.join(LOC, AF.REL_AF, cid); n = nh = nores = 0; uw = 0; exp = sum(len(v) for v in c["seeds"].values())
            for lab, _ in R.POS:
                for s in c["seeds"][lab]:
                    f = os.path.join(d0, f"x_{lab}", f"red_{s}.csv")
                    if not os.path.exists(f): continue
                    n += 1; r = pd.read_csv(f).iloc[0]; uw += (r["u_wall_max"] != 0.0) or (r["W_div"] != 0.0)
                    log = open(os.path.join(d0, f"x_{lab}", f"run_{s}.log"), errors="ignore").read()
                    nh += log.count("[EDMD-HEALTH]"); nores += "[EDMD-RESCHED] divider events: minimal" not in log
            bld = ",".join(sorted(builds_af)); eh = er = float("nan"); led = "yes" if uw == 0 else f"**NO** ({uw} seeds)"
            ok = n == exp and nh == 0 and nores == 0 and uw == 0 and (DRY or "279282b" not in bld)
        allok &= ok
        print(f"| {cid} | {mode} | {n} ({exp}) | {nh} | {nores} | {bld} | {eh:.2e} | {er:.2e} | {led} | {'yes' if ok else '**NO**'} |")
    print(f"\nG-E3: {'PASS' if allok else 'FAIL'}" + (" (dry run on the old data: the [EDMD-RESCHED]/energy columns cannot pass)" if DRY else ""))
    return allok


def g_e4(CS):
    rec_c = pd.read_csv(os.path.join(T.PLOTS, "261004_p1_confinement_cells.csv"), dtype=str).set_index("cell")
    rec_a = pd.read_csv(os.path.join(T.PLOTS, "261005_p1_identity_afix_cells.csv"), dtype=str).set_index("cell")
    rows = []; repro = []; builds_af = set()
    for mode, cid in CELLS:
        with contextlib.redirect_stdout(io.StringIO()):
            c_old = dict(CS[cid]); c_new = dict(CS[cid])
            if mode == "B":
                R.method_A(c_old); old = b_values(c_old)
                for k in ("static", "s_kT", "kT", "F2term"): c_new[k] = c_old[k]   # identity() needs the static side; only kS is used
                with root(R, LOC): new = b_values(c_new)
                repro += [(cid, "c_s", old["cs"][0], rec_c.loc[cid, "c_s"]), (cid, "c_s_err", old["cs"][1], rec_c.loc[cid, "c_s_err"]),
                          (cid, "k_S_dyn", old["kS"][0], rec_c.loc[cid, "k_S_dyn"])]
                names = (("cs", "c_s"), ("kS", "k_S^dyn"), ("G5", "Gamma(alpha=5)"))
            else:
                old, _ = af_values(c_old)
                with root(AF, LOC): new, a_new = af_values(c_new)
                builds_af = a_new["builds"]
                repro += [(cid, "k_T_afix", old["kT"][0], rec_a.loc[cid, "k_T_afix"]), (cid, "F_L0", old["F0"][0], rec_a.loc[cid, "F_L0"]),
                          (cid, "static_afix", old["static"][0], rec_a.loc[cid, "static_afix"])]
                names = (("F0", "F(L_0)"), ("kT", "k_T"), ("static", "k_T + F^2/(N_s kT)"))
        for k, lab in names:
            (o, so), (n, sn) = old[k], new[k]; sd = math.hypot(so, sn)
            rows.append((f"{mode} {cid}", lab, o, so, n, sn, n - o, sd, (n - o) / sd if sd > 0 else 0.0, (n - o) / so if so > 0 else 0.0))
    print("### G-E4 -- old (279282b) vs new (replay): the nine numbers\n")
    print("| # | cell | quantity | old | SE old | new | SE new | new - old | sigma_diff | z | abs(z) < 2 | (new - old)/SE_old |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(rows, 1):
        print(f"| {i} | {r[0]} | {r[1]} | {r[2]:.6g} | {r[3]:.3g} | {r[4]:.6g} | {r[5]:.3g} | {r[6]:+.3g} | {r[7]:.3g} | {r[8]:+.2f} | "
              f"{'yes' if abs(r[8]) < 2 else '**NO**'} | {r[9]:+.2f} |")
    zs = np.array([r[8] for r in rows]); chi2 = float((zs ** 2).sum())
    p1 = 2 * norm.cdf(2) - 1; pall = p1 ** len(rows)
    print(f"\nsum z^2 = {chi2:.2f} for {len(rows)} numbers (correlated: c_s with k_S^dyn, k_T with k_static)")
    print(f"false-fail probability of 'all {len(rows)} abs(z) < 2' under the null, if independent: 1 - {p1:.4f}^{len(rows)} = {1 - pall:.3f}; "
          f"with the per-number limit 2.77 (Bonferroni, family-wise 5 %): {1 - (2 * norm.cdf(2.77) - 1) ** len(rows):.3f}")
    ok = bool(np.all(np.abs(zs) < 2))
    print(f"\nG-E4 (as registered, all abs(z) < 2): {'PASS' if ok else 'FAIL'}; for information, all abs(z) < 2.77: "
          f"{'yes' if np.all(np.abs(zs) < 2.77) else 'no'}")
    print("\n### Reproduction: the old values recomputed here against the recorded CSVs (at the CSV's printed precision)\n")
    print("| cell | column | recomputed | recorded | equal at the recorded decimals |\n|---|---|---|---|---|")
    rep_ok = True
    for cid, col, v, r in repro:
        dec = len(r.split(".")[1]) if "." in r else 0; eq = abs(v - float(r)) <= 0.5 * 10 ** -dec * (1 + 1e-9); rep_ok &= eq
        print(f"| {cid} | {col} | {v:.10g} | {r} | {'yes' if eq else '**NO**'} |")
    return ok, rep_ok, builds_af


def g_e2():
    import ge2
    ds = sorted(glob.glob(os.path.join(LOC, "resched_gate_261005", "ge2_*")))
    if not ds: print("### G-E2 -- no fetched G-E2 directory yet\n"); return None
    class A: out = ds[-1]
    return ge2.compare(A) == 0


def g_e5():
    fs = sorted(glob.glob(os.path.join(LOC, "profile_edmd_*", "times.tsv")))
    if not fs: print("### G-E5 -- no fetched profile yet\n"); return
    t = pd.read_csv(fs[-1], sep="\t"); print(f"### G-E5 -- profile ({os.path.dirname(fs[-1])})\n")
    print("| kind | N | minimal [s] | legacy [s] | 279282b job 14983181 [s] | legacy / minimal |\n|---|---|---|---|---|---|")
    for kind in ("held", "free"):
        for N in (100, 400):
            m = t[(t.kind == kind) & (t.N == N) & (t.policy == "minimal")].wall_s.iloc[0]
            l = t[(t.kind == kind) & (t.N == N) & (t.policy == "legacy")].wall_s.iloc[0]
            o = OLD_PROFILE.get(f"{kind} {N} legacy", float("nan"))
            print(f"| {kind} | {N} | {m:.2f} | {l:.2f} | {o:.1f} | {l / m:.2f} |")
    for kind in ("held", "free"):
        for pol in ("minimal", "legacy"):
            a = t[(t.kind == kind) & (t.N == 100) & (t.policy == pol)].wall_s.iloc[0]; b = t[(t.kind == kind) & (t.N == 400) & (t.policy == pol)].wall_s.iloc[0]
            print(f"exponent {kind} {pol}: p = ln({b:.2f}/{a:.2f})/ln 4 = {math.log(b / a) / math.log(4):.2f}")


def main():
    R.ALLOW_NEW_BUILD = True   # explicit: this script compares the two build generations (G-E6)
    with contextlib.redirect_stdout(io.StringIO()):
        CS = {c["cid"]: c for c in R.cells()}
    print(f"## Engine gate (261012 sec. 4.4) -- replay tree {LOC}{'  [DRY RUN: new := old]' if DRY else ''}\n")
    r2 = None if DRY else g_e2()
    ok4, rep_ok, builds_af = g_e4(CS)
    print(); ok3 = g_e3(CS, builds_af)
    print(); g_e5() if not DRY else None
    print(f"\nSUMMARY: G-E2 {'n/a' if r2 is None else ('PASS' if r2 else 'FAIL')}; G-E3 {'PASS' if ok3 else 'FAIL'}; "
          f"G-E4 {'PASS' if ok4 else 'FAIL'}; reproduction of the recorded old values {'yes' if rep_ok else 'NO'}")


if __name__ == "__main__":
    main()
