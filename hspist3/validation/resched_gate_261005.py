#!/usr/bin/env python3
"""##CHRIS 2026-10-05: the engine gate of the minimal divider rescheduling (branch engine-divider-resched, 261012 sec. 4.4),
Mac side, after hspist3/cluster/resched_gate_261005/fetch_resched.sh. Prints G-E2 (the comparison of the KOA G-E2 job, the same
code as on KOA), G-E3 (the three replayed cells clean), G-E4 (the nine numbers against the 279282b campaign), G-E5 (profile)
and the one-build check (every gate output from one clean build).

Replayed cells (the 279282b campaign's own task lines: same masses, positions, seeds, protocol), read from the separate tree
hspist3/experiments_resched_gate_261005/ (never mixed with the canonical one). The registered estimators read it only inside
new_tree(), which also sets R.ALLOW_NEW_BUILD for that read alone, so the 279282b guard stays active on the old side.
  B  e0p10_H_H10_L39.25  (eta 0.10 anchor)   B  epi8_H_H10_L10  (pi/8 anchor)   AF  epi8_H_H10_L10  (pi/8 anchor, held divider)

Registered before the replay data exist (261012 sec. 4.4; amended 2026-10-05 after the branch review, still before any data):
  G-E3  every trajectory of the task list present, matched by run index and seed (B) or seed (AF), with a finite estimator
        window n > 0 (B); the number of '##RUN' log sections equals the number of trajectories (B); no [EDMD-HEALTH] line in
        any log, failed trajectories included (B: .failed_run*/stdout.log, AF: every run_<seed>.log); the counters behind it
        are overlap_repairs, wall_overdue, forced_advance, clamp repairs and past_events (= events popped with a time before
        the current time, a heap-order violation; negative collision times cannot occur: the solvers return t > 1e-12, or
        t = 0 for an overdue contact, counted as overlap_repairs / wall_overdue); every log carries "[EDMD-RESCHED] divider events: minimal"; divider ledgers: B -- every trajectory's
        |dE/E| (E = gas + divider KE, [EDMD-ENERGY]) <= 1e-10 in the hold and in the record; AF -- u_wall_max = 0 and W_div = 0
        in every seed; contact audit -- every trajectory's [EDMD-CONTACT] maxima <= 1e-6 px in all four classes. Any value
        that cannot be evaluated (missing or not finite) FAILS.
  G-E4  nine numbers, three per cell, each with the registered estimator on the 279282b data (old) and on the replay (new):
          B cells:  c_s (through-origin slope, statistical seed error c_s_err), k_S^dyn (heavy masses, C1),
                    Gamma at alpha = 5 (energy-decay rate 2/tau_r of the divider mode, jackknife error);
          AF cell:  F(L_0) = (F_L + F_R)/2 at x = 0, k_T (five-point stencil), k_static = k_T + F^2/(N_s kT).
        z = (new - old) / sqrt(SE_old^2 + SE_new^2): with the same seeds but a different build the trajectories decorrelate
        within the 200 sigma-time equilibration, so old and new are independent realisations. PASS if all nine |z| < 2 (the
        plan author's rule); a z that cannot be computed FAILS. Also printed: (new - old)/SE_old and the false-fail
        probability of the rule under the null.
  G-E5  wall times only from runs that exited 0; a row with a failed run is marked invalid.
  BUILD every gate output names one clean build (no -dirty): B .build_git files, AF .build_git files and summaries, the
        replay root's .build_generation, the G-E2 version.txt files of the new binary, the profile summaries.
  The old values are also checked against the recorded CSVs (261004_p1_confinement_cells.csv, 261005_p1_identity_afix_cells.csv).
usage (from hspist3/):  python3 validation/resched_gate_261005.py            [--dry-run-old: "new" := the old data, plumbing test]
"""
import contextlib, glob, io, math, os, re, subprocess, sys
import numpy as np, pandas as pd
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS); sys.path.insert(0, os.path.join(HS, "cluster", "resched_gate_261005"))
import paper1_confinement_results_261004 as R
import paper1_confinement_afix_261005 as AF
import tests_20260913 as T
import edmd_acc_guard   # ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): the loader provenance guard (full name: no alias can be shadowed)

DRY = "--dry-run-old" in sys.argv
LOC = HS if DRY else os.environ.get("HD_RESCHED_LOC", os.path.join(HS, "experiments_resched_gate_261005"))   # env: test hook only
CELLS = (("B", "e0p10_H_H10_L39.25"), ("B", "epi8_H_H10_L10"), ("AF", "epi8_H_H10_L10"))
E_TOL, C_TOL = 1e-10, 1e-6
EN = re.compile(r"\[EDMD-ENERGY\] (\d) .*?E_tot=(\S+) resched=(\w+)")
CONTACT = re.compile(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), "
                     r"outer walls (\S+), divider (\S+), pistons (\S+)\n")
GIT = re.compile(r"git\s+([0-9a-f]{7,40}(?:-dirty)?)")
OLD_PROFILE = {"held 100 legacy": 2.3, "held 400 legacy": 56.4}   # 279282b, KOA job 14983181 (261012 sec. 4.3) [DATA]
CONF = os.path.join(HS, "cluster", "confinement_20261013")


@contextlib.contextmanager
def new_tree():
    """Read the replay tree with the registered estimators: R.HS and AF.HS point to LOC, and the 279282b guard of R is
    lifted for this read only (explicit flag, G-E6)."""
    old = (R.HS, AF.HS, R.ALLOW_NEW_BUILD); R.HS = LOC; AF.HS = LOC; R.ALLOW_NEW_BUILD = True
    try: yield
    finally: R.HS, AF.HS, R.ALLOW_NEW_BUILD = old


def fin(*xs):
    return all(isinstance(x, (int, float, np.floating)) and math.isfinite(float(x)) for x in xs)


def b_values(c):
    """c_s, k_S^dyn, Gamma(alpha = 5) with errors, for the current R.HS."""
    R.method_B(c); R.identity(c); R.damping(c)
    g = [d for d in c["damp"] if abs(d["alpha"] - 5.0) < 1e-9][0]
    sg = g["Gamma"] * g["s_tau_r"] / g["tau_r"] if g.get("ok") else float("nan")
    return dict(cs=(c["cs"], c["cs_err"]), kS=(c["kS"], c["s_kS"]), G5=(g.get("Gamma", float("nan")), sg))


def af_values(c):
    a = AF.inventory_and_static(c)
    if a.get("error"): sys.exit(f"AF {c['cid']}: {a['error']}")
    s0 = a["st"][0]; sF = 0.5 * math.hypot(s0["sFL"], s0["sFR"])
    s_st = math.hypot(a["s_kT"], 2 * a["F0"] * sF / (c["Ns"] * a["temp"]))
    return dict(F0=(a["F0"], sF), kT=(a["kT"], a["s_kT"]), static=(a["static"], s_st)), a


def tasks(mode, cid):
    """The campaign's task list: B -> {M: {(run, seed)}}, AF -> {position label: [seeds]}."""
    out = {}
    for l in open(os.path.join(CONF, f"tasks_{mode}_{cid}.txt")):
        f = l.split()
        if mode == "B": out.setdefault(int(f[2]), set()).add((int(f[3]), int(f[4])))
        else: out.setdefault(os.path.basename(f[1])[2:], []).append(int(f[3]))
    return out


def contact_ok(text):
    m = CONTACT.findall(text)
    if len(m) != 1: return False, float("nan")
    g = [float(x) for x in m[0][1:]]
    return (all(math.isfinite(x) and x <= C_TOL for x in g)), max(g)


def g_e3(CS):
    print("### G-E3 -- the replayed cells are clean\n")
    print("| cell | method | trajectories matched (expected) | log sections | health lines (failed runs incl.) | logs without the minimal "
          "policy line | max abs(dE/E) hold | max abs(dE/E) record | max contact gap [px] | divider ledger | PASS |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    allok = True
    for mode, cid in CELLS:
        tk = tasks(mode, cid)
        if mode == "B":
            d0 = os.path.join(LOC, R.REL_B, cid); matched = exp = nsec = nh = nores = bad = 0; eh = er = cg = 0.0
            for M, want in tk.items():
                exp += len(want); d = os.path.join(d0, f"m_{M}")
                rn = pd.read_csv(edmd_acc_guard.guard(os.path.join(d, "red_nu.csv")))
                have = {(int(r), int(s)) for r, s in zip(rn["run"], rn["seed"])}
                matched += len(want & have) if len(rn) == len(have) else 0
                bad += int((~np.isfinite(rn["n"].astype(float)) | (rn["n"].astype(float) <= 0)).sum())
                secs = open(edmd_acc_guard.guard(os.path.join(d, "run.log")), errors="ignore").read().split("##RUN")[1:]
                nsec += len(secs); bad += abs(len(secs) - len(rn))
                for sec in secs:
                    nh += sec.count("[EDMD-HEALTH]"); nores += "[EDMD-RESCHED] divider events: minimal" not in sec
                    e = {int(m.group(1)): float(m.group(2)) for m in EN.finditer(sec)}
                    c_ok, c_max = contact_ok(sec)
                    if set(e) != {0, 1, 2} or not c_ok and not DRY: bad += 1
                    if set(e) == {0, 1, 2}:
                        h, r = abs(e[1] / e[0] - 1), abs(e[2] / e[1] - 1)
                        if not fin(h, r): bad += 1; continue
                        eh, er = max(eh, h), max(er, r)
                    if fin(c_max): cg = max(cg, c_max)
                for f in glob.glob(os.path.join(d, ".failed_run*", "stdout.log")):
                    nh += open(edmd_acc_guard.guard(f), errors="ignore").read().count("[EDMD-HEALTH]")
            led = "energy (above)"
            ok = matched == exp and nsec == exp and nh == 0 and nores == 0 and bad == 0 and eh <= E_TOL and er <= E_TOL
        else:
            d0 = os.path.join(LOC, AF.REL_AF, cid); matched = nh = nores = bad = uw = 0; exp = sum(len(v) for v in tk.values())
            cg = 0.0; eh = er = float("nan")
            for lab, seeds in tk.items():
                for lg in glob.glob(os.path.join(d0, f"x_{lab}", "run_*.log")):
                    nh += open(edmd_acc_guard.guard(lg), errors="ignore").read().count("[EDMD-HEALTH]")
                for s in seeds:
                    f = os.path.join(d0, f"x_{lab}", f"red_{s}.csv")
                    if not os.path.exists(f): continue
                    matched += 1; r = pd.read_csv(edmd_acc_guard.guard(f)).iloc[0]; uw += (r["u_wall_max"] != 0.0) or (r["W_div"] != 0.0)
                    log = open(edmd_acc_guard.guard(os.path.join(d0, f"x_{lab}", f"run_{s}.log")), errors="ignore").read()
                    nores += "[EDMD-RESCHED] divider events: minimal" not in log
                    c_ok, c_max = contact_ok(log)
                    if not c_ok and not DRY: bad += 1
                    if fin(c_max): cg = max(cg, c_max)
            nsec = matched; led = "yes" if uw == 0 else f"**NO** ({uw} seeds)"
            ok = matched == exp and nh == 0 and nores == 0 and uw == 0 and bad == 0
        allok &= ok
        print(f"| {cid} | {mode} | {matched} ({exp}) | {nsec} | {nh} | {nores} | {eh:.2e} | {er:.2e} | {cg:.2e} | {led} | "
              f"{'yes' if ok else '**NO**'}{'' if not bad else f' ({bad} unusable)'} |")
    print(f"\nG-E3: {'PASS' if allok else 'FAIL'}" + (" (dry run on the old data: the policy/energy/contact columns cannot pass)" if DRY else ""))
    return allok


def g_e4(CS):
    rec_c = pd.read_csv(os.path.join(T.PLOTS, "261004_p1_confinement_cells.csv"), dtype=str).set_index("cell")
    rec_a = pd.read_csv(os.path.join(T.PLOTS, "261005_p1_identity_afix_cells.csv"), dtype=str).set_index("cell")
    rows = []; repro = []; builds_af = set(); permass = []
    for mode, cid in CELLS:
        with contextlib.redirect_stdout(io.StringIO()):
            c_old = dict(CS[cid]); c_new = dict(CS[cid])
            if mode == "B":
                R.method_A(c_old); old = b_values(c_old)
                for k in ("static", "s_kT", "kT", "F2term"): c_new[k] = c_old[k]   # identity() needs the static side; only kS is used
                with new_tree(): new = b_values(c_new)
                for ro, rn in zip(c_old["B"], c_new["B"]):   # same masses, same order (c["Ms"])
                    permass.append((cid, ro["M"], ro["alpha"], ro["nu"], ro["se"], ro["n"], rn["nu"], rn["se"], rn["n"]))
                repro += [(cid, "c_s", old["cs"][0], rec_c.loc[cid, "c_s"]), (cid, "c_s_err", old["cs"][1], rec_c.loc[cid, "c_s_err"]),
                          (cid, "k_S_dyn", old["kS"][0], rec_c.loc[cid, "k_S_dyn"])]
                names = (("cs", "c_s"), ("kS", "k_S^dyn"), ("G5", "Gamma(alpha=5)"))
            else:
                old, _ = af_values(c_old)
                with new_tree(): new, a_new = af_values(c_new)
                builds_af = a_new["builds"]
                repro += [(cid, "k_T_afix", old["kT"][0], rec_a.loc[cid, "k_T_afix"]), (cid, "F_L0", old["F0"][0], rec_a.loc[cid, "F_L0"]),
                          (cid, "static_afix", old["static"][0], rec_a.loc[cid, "static_afix"])]
                names = (("F0", "F(L_0)"), ("kT", "k_T"), ("static", "k_T + F^2/(N_s kT)"))
        for k, lab in names:
            (o, so), (n, sn) = old[k], new[k]; sd = math.hypot(so, sn) if fin(so, sn) else float("nan")
            z = (n - o) / sd if fin(sd, n, o) and sd > 0 else float("nan")
            rows.append((f"{mode} {cid}", lab, o, so, n, sn, n - o, sd, z, (n - o) / so if fin(so) and so > 0 else float("nan")))
    print("### G-E4 -- old (279282b) vs new (replay): the nine numbers\n")
    print("| # | cell | quantity | old | SE old | new | SE new | new - old | sigma_diff | z | abs(z) < 2 | (new - old)/SE_old |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(rows, 1):
        verdict = "**not computable**" if not fin(r[8]) else ("yes" if abs(r[8]) < 2 else "**NO**")
        print(f"| {i} | {r[0]} | {r[1]} | {r[2]:.6g} | {r[3]:.3g} | {r[4]:.6g} | {r[5]:.3g} | {r[6]:+.3g} | {r[7]:.3g} | {r[8]:+.2f} | "
              f"{verdict} | {r[9]:+.2f} |")
    zs = np.array([r[8] for r in rows], float); good = np.isfinite(zs)
    p1 = 2 * norm.cdf(2) - 1; pall = p1 ** len(rows)
    print(f"\nsum z^2 = {float(np.nansum(zs ** 2)):.2f} for {int(good.sum())} computable of {len(rows)} numbers "
          f"(correlated: c_s with k_S^dyn, k_T with k_static)")
    print(f"false-fail probability of 'all {len(rows)} abs(z) < 2' under the null, if independent: 1 - {p1:.4f}^{len(rows)} = {1 - pall:.3f}; "
          f"with the per-number limit 2.77 (Bonferroni, family-wise 5 %): {1 - (2 * norm.cdf(2.77) - 1) ** len(rows):.3f}")
    ok = bool(good.all() and np.all(np.abs(zs) < 2))
    print(f"\nG-E4 (as registered, all abs(z) < 2): {'PASS' if ok else 'FAIL'}; for information, all abs(z) < 2.77: "
          f"{'yes' if good.all() and np.all(np.abs(zs) < 2.77) else 'no'}")
    print("\n### Per-mass divider frequency, old vs new (information only, added 2026-10-06 before the replay data; the failed KOA smoke")
    print("### test (261012 sec. 4.4.7) had its deficit at alpha = 0.5 and 1, so those rows test a light-mass shift directly)\n")
    print("| cell | M | alpha | old nu (25 seeds) | SE | new nu | SE | (new - old)/old [%] | z |\n|---|---|---|---|---|---|---|---|---|")
    for cid, M, al, no, so, nno, nn, sn, nnn in permass:
        sd = math.hypot(so, sn); z = (nn - no) / sd if fin(sd) and sd > 0 else float("nan")
        print(f"| {cid} | {M} | {al:g} | {no:.6f} | {so:.6f} | {nn:.6f} | {sn:.6f} | {100 * (nn - no) / no:+.2f} | {z:+.2f}"
              f"{' **light**' if al <= 1.0 else ''} |")
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
    if not fs: print("### G-E5 -- no fetched profile yet\n"); return None
    t = pd.read_csv(fs[-1], sep="\t"); print(f"### G-E5 -- profile ({os.path.dirname(fs[-1])})\n")
    wall = {}
    for _, r in t.iterrows():
        wall[(r["kind"], int(r["N"]), r["policy"])] = float(r["wall_s"]) if int(r["exit"]) == 0 else float("nan")
        if int(r["exit"]) != 0: print(f"**INVALID: {r['kind']} N = {int(r['N'])} {r['policy']} exited {int(r['exit'])}; its time is not used**")
    print("| kind | N | minimal [s] | legacy [s] | 279282b job 14983181 [s] | legacy / minimal |\n|---|---|---|---|---|---|")
    for kind in ("held", "free"):
        for N in (100, 400):
            m, l = wall.get((kind, N, "minimal"), float("nan")), wall.get((kind, N, "legacy"), float("nan"))
            print(f"| {kind} | {N} | {m:.2f} | {l:.2f} | {OLD_PROFILE.get(f'{kind} {N} legacy', float('nan')):.1f} | "
                  f"{l / m if fin(l, m) else float('nan'):.2f} |")
    for kind in ("held", "free"):
        for pol in ("minimal", "legacy"):
            a, b = wall.get((kind, 100, pol), float("nan")), wall.get((kind, 400, pol), float("nan"))
            print(f"exponent {kind} {pol}: p = ln({b:.2f}/{a:.2f})/ln 4 = {math.log(b / a) / math.log(4) if fin(a, b) else float('nan'):.2f}")
    return True


def one_build(builds_af):
    """Every gate output from one clean build; prints where each build name came from."""
    found = {}
    def add(src, text):
        m = GIT.search(text or "")
        found.setdefault(m.group(1) if m else f"(unreadable: {str(text).strip()[:40]})", []).append(src)
    for mode, cid in CELLS:
        base = os.path.join(LOC, R.REL_B if mode == "B" else AF.REL_AF, cid)
        for f in sorted(glob.glob(os.path.join(base, "*", ".build_git"))): add(os.path.relpath(f, LOC), open(f).read())
    for b in sorted(builds_af): add("AF summaries", f"git {b}")
    g = os.path.join(LOC, ".build_generation")
    if os.path.exists(g): add(".build_generation", open(g).read())
    for f in sorted(glob.glob(os.path.join(LOC, "resched_gate_261005", "ge2_*", "*_m*", "version.txt"))
                    + glob.glob(os.path.join(LOC, "resched_gate_261005", "ge2_*", "*_legacy", "version.txt"))):
        add(os.path.relpath(f, LOC), open(f).read())
    for f in sorted(glob.glob(os.path.join(LOC, "profile_edmd_*", "*", "summary.csv"))):
        try: add(os.path.relpath(f, LOC), "git " + str(pd.read_csv(f)["build_git"].iloc[-1]))
        except Exception as ex: add(os.path.relpath(f, LOC), f"unreadable {ex}")
    print("### Build check -- one clean build behind every gate output\n\n| build | sources |\n|---|---|")
    for b, srcs in found.items(): print(f"| {b} | {len(srcs)}: {', '.join(srcs[:3])}{' ...' if len(srcs) > 3 else ''} |")
    keys = list(found)
    same = len(keys) >= 1 and all(k.startswith(keys[0][:7]) and not k.endswith("-dirty") and "unreadable" not in k for k in keys)
    try: head = subprocess.run(["git", "-C", HS, "rev-parse", "--short", "engine-divider-resched"], capture_output=True, text=True).stdout.strip()
    except Exception: head = "?"
    print(f"\nlocal branch head engine-divider-resched: {head}; the gate's build is "
          f"{'the same' if same and head and keys[0].startswith(head[:7]) else 'NOT the local branch head -- check which commit KOA built'}")
    print(f"BUILD: {'PASS (one clean build: ' + keys[0] + ')' if same else 'FAIL (more than one build, a -dirty build, or an unreadable record)'}")
    return same


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        CS = {c["cid"]: c for c in R.cells()}
    print(f"## Engine gate (261012 sec. 4.4) -- replay tree {LOC}{'  [DRY RUN: new := old]' if DRY else ''}\n")
    r2 = None if DRY else g_e2()
    ok4, rep_ok, builds_af = g_e4(CS)
    print(); ok3 = g_e3(CS)
    print(); r5 = None if DRY else g_e5()
    print(); okb = None if DRY else one_build(builds_af)
    print(f"\nSUMMARY: G-E2 {'n/a' if r2 is None else ('PASS' if r2 else 'FAIL')}; G-E3 {'PASS' if ok3 else 'FAIL'}; "
          f"G-E4 {'PASS' if ok4 else 'FAIL'}; G-E5 {'n/a' if r5 is None else 'printed'}; build {'n/a' if okb is None else ('PASS' if okb else 'FAIL')}; "
          f"reproduction of the recorded old values {'yes' if rep_ok else 'NO'}")


if __name__ == "__main__":
    main()
