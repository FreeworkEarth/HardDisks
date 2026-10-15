#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.14): every table of stage A (M3), printed from the evidence folder of stageA_evidence.py.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m3_261009/stageA_tables.py --bin-dir <frozen binaries> --out <evidence dir> > stageA_tables_output.txt
"""
import argparse, glob, hashlib, math, os, re, shutil, subprocess, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); WT = os.path.dirname(HS)
MAIN_HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
GATE = os.path.join(MAIN_HS, "cluster", "resched_gate_261005")
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REC = re.compile(r"\[EDMD3-HEALTH\] (.*?): clean=(\S+) (.*)")
CONTACT = re.compile(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), outer walls (\S+), divider (\S+), pistons (\S+)")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rec(log):
    t = open(log, errors="ignore").read(); m = REC.findall(t)
    if not m: return None
    rid, clean, rest = m[-1]
    kv = dict(re.findall(r"(\w+)=(\"[^\"]*\"|\S+)", rest)); kv["clean"] = clean; kv["id"] = rid; kv["n_records"] = len(m)
    kv["n_builds"] = len(re.findall(r"\[EDMD3\] built #", t))
    c = CONTACT.search(t); kv["contact"] = tuple(float(x) for x in c.groups()[1:]) if c else None
    return kv
def events(k): return sum(int(k[f]) for f in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# Stage A (M3) evidence tables, printed by experiments_gen3_m3_261009/stageA_tables.py\n")
    # ------------------------------------------------------------ binaries
    print("## 1. The frozen binaries (rule 5): built from the committed tree by build_clean.sh\n")
    rec_txt = open(os.path.join(HERE, "build_record_f42befb.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec_txt.split("## SHA-256")[1], re.M))
    print("| binary | SHA-256 (frozen copy) | = build record |\n|---|---|---|")
    for n, h in want.items():
        got = sha(os.path.join(B, n)); print(f"| {n} | {got} | {'yes' if got == h else '**NO**'} |")
    v = subprocess.run([os.path.join(B, "00ALLINONE"), "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    print(f"\n`00ALLINONE --version`: {v} ({'no -dirty' if '-dirty' not in v else '**-dirty**'})")
    # ------------------------------------------------------------ e0
    print("\n## 2. Acceptance 1 / rule 4: gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg)\n")
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "e0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "e0", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "e0", f"report_{tag}.txt"), "w").write(r)
        print(f"### {'default build' if tag == 'default' else '--engine=gen2 (every command plus the flag)'}\n")
        keep = False
        for l in r.splitlines():          # the report's two tables, as printed by audit_runs_261007.py report
            if l.startswith("| case |"): keep = True; print()
            if keep and l.startswith("|"): print(l)
            elif keep and not l.startswith("|"): keep = False
        print()
    # ------------------------------------------------------------ harness
    print("## 3. Rule 8 after the engine changes: the M1 and M2 harness outputs, byte for byte\n")
    H = os.path.join(O, "harness")
    print("| output | this build | committed output | SHA-256 (this build) | identical |\n|---|---|---|---|---|")
    for f, ref in (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
                   ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
                   ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt")):
        p, q = os.path.join(H, f), os.path.join(HS, ref)
        ok = os.path.exists(p) and open(p, "rb").read() == open(q, "rb").read()
        print(f"| {f} | {os.path.getsize(p) if os.path.exists(p) else 'missing'} bytes | {ref} | {sha(p) if os.path.exists(p) else '-'} | {'IDENTICAL' if ok else '**DIFFERENT**'} |")
    p, q = os.path.join(H, "body_rule_test_output.txt"), os.path.join(HS, "experiments_gen3_m2_261009/body_rule_test_output.txt")
    a_, b_ = open(p).read(), open(q).read()
    l1 = a_.split("\n", 1)[0]; same1 = l1 == b_.split("\n", 1)[0]
    print(f"| body_rule_test_output.txt | {len(a_)} bytes | experiments_gen3_m2_261009/body_rule_test_output.txt | {sha(p)} | "
          f"random-case line (the first line) {'IDENTICAL' if same1 else '**DIFFERENT**'}; the committed file's second line 'exit=0' was "
          f"appended by the M2 run command; new: the constructed categories of amendment b |")
    print("\n### Amendment b: the constructed categories (gen3_body_rule_test 20000, the part after the random cases)\n")
    print(a_.split("\n", 1)[1].strip())
    print("\n### Amendment a: the band-edge stress cells (gen3_band_edge_test)\n")
    print(re.sub(r"^\[EDMD3-AUDIT\].*\n", "", open(os.path.join(H, "band_edge_output.txt")).read(), flags=re.M).strip())
    na = len(re.findall(r"^\[EDMD3-AUDIT\]", open(os.path.join(H, "band_edge_output.txt")).read(), re.M))
    print(f"\n(audit lines printed by the engine, not repeated here: {na})")
    print("\n### Amendment a: the numbers behind each audit finding of the exact cells (gen3_band_edge_ties)\n")
    print(open(os.path.join(H, "band_edge_ties_output.txt")).read().split("\n", 2)[2].strip())
    print("\n" + open(os.path.join(H, "run.txt")).read().strip())
    # ------------------------------------------------------------ replays
    print("\n## 4. Acceptances 2 and 3: the harness cells replayed through the driver (run_replays.sh)\n")
    t = open(os.path.join(O, "replay", "replay_output.txt")).read()
    rows = [l for l in t.splitlines() if re.match(r"^[abcd] m[12]_", l)]
    print("| variant | replays | MATCH | clean=1 |\n|---|---|---|---|")
    for v in "abcd":
        r = [l for l in rows if l.startswith(v + " ")]
        print(f"| {v} | {len(r)} | {sum(': MATCH;' in l for l in r)} | {sum('clean=1' in l for l in r)} |")
    print("\n" + t.strip().splitlines()[-1])
    # ------------------------------------------------------------ prod
    print("\n## 5. Amendment d: one production-length run through the driver (N = 400, eta 0.70, free divider M = 500, 2e4 sigma-time)\n")
    lg = os.path.join(O, "prod", "run.log"); k = rec(lg); t = open(lg, errors="ignore").read()
    for tag in ("[EDMD-CONTACT]", "[EDMD3-LEDGER]", "[EDMD3-GAP]", "[EDMD3-AUDIT] audited"):
        for l in t.splitlines():
            if l.startswith(tag): print(f"    {l}")
    print(f"\n| clean | t_end [sigma] | origin shifts | events | validator every [steps] | schedule audits | build | run [s] | engine [s] | driver share |\n|---|---|---|---|---|---|---|---|---|---|")
    na = re.search(r"audited states (\d+)", t).group(1)
    print(f"| {k['clean']} | {float(k['t_end']) / 24:.1f} | {k['origin_shifts']} | {events(k)} | {k['validator_every']} | {na} | {k['build']} {k['target']} | "
          f"{k['run_s']} | {k['engine_s']} | {k['driver_share']} |")
    print(f"\nhealth line: {open(lg).read().split('[EDMD3-HEALTH] ')[-1].split(' engine=')[0]}")
    # ------------------------------------------------------------ a3, acc3, acc6
    print("\n## 6. Stage A3 and acceptance 3: the event hash at every output cadence; acceptance 6: same-seed determinism\n")
    print("| loop | run | validator every | trace / log / psi6 cadence | event hash | events | t_end [sigma] | clean | driver share |\n|---|---|---|---|---|---|---|---|---|")
    for tag, cad in (("default", "trace every 600 steps"), ("validator1", "trace every 600 steps"), ("validator1_trace1", "trace every step")):
        k = rec(os.path.join(O, "a3", tag, "run.log"))
        print(f"| energy transfer (A-fixed protocol) | {tag} | {k['validator_every']} | {cad} | {k['hash']} | {events(k)} | {float(k['t_end']) / 24:.4f} | {k['clean']} | {k['driver_share']} |")
    for tag, cad in (("default", "log 60 steps, psi6 0.25"), ("max_cadence", "log every step, psi6 0.05"), ("min_cadence", "log 600 steps, psi6 1000"),
                     ("default_repeat", "log 60 steps, psi6 0.25")):
        k = rec(os.path.join(O, "acc3", tag, "run.log"))
        print(f"| speed of sound (T-prime protocol, M = 300) | {tag} | {k['validator_every']} | {cad} | {k['hash']} | {events(k)} | {float(k['t_end']) / 24:.4f} | {k['clean']} | {k['driver_share']} |")
    a3 = {t_: os.path.join(O, "a3", t_) for t_ in ("default", "validator1", "validator1_trace1")}
    same_ev = all(open(os.path.join(a3["default"], "ev_9700.csv"), "rb").read() == open(os.path.join(a3[t_], "ev_9700.csv"), "rb").read() for t_ in a3)
    tr0 = open(os.path.join(a3["default"], "tr_9700.csv")).read().splitlines(); tr1 = open(os.path.join(a3["validator1"], "tr_9700.csv")).read().splitlines()
    trf = open(os.path.join(a3["validator1_trace1"], "tr_9700.csv")).read().splitlines()
    sub = set(tr0[1:]) <= set(trf[1:]) and tr0[0] == trf[0]
    print(f"\nA3: event logs (ev_9700.csv) of the three runs byte-identical: {'yes' if same_ev else '**NO**'}; trace default = validator1: "
          f"{'yes' if tr0 == tr1 else '**NO**'}; every row of the default trace ({len(tr0) - 1}) is a row of the every-step trace ({len(trf) - 1}): {'yes' if sub else '**NO**'}")
    d0, d1 = os.path.join(O, "acc3", "default"), os.path.join(O, "acc3", "default_repeat")
    fs = sorted(os.path.basename(f) for f in glob.glob(os.path.join(d0, "*.csv")))
    same = [f for f in fs if open(os.path.join(d0, f), "rb").read() == open(os.path.join(d1, f), "rb").read()]
    mask = lambda s: re.sub(r"run_s=\S+ engine_s=\S+ driver_share=\S+", "", re.sub(r"\d+(\.\d+)? ?(s|ms|sec|seconds)\b", "", s))
    l0 = [l for l in open(os.path.join(d0, "run.log"), errors="ignore") if l.startswith("[EDMD3")]
    l1 = [l for l in open(os.path.join(d1, "run.log"), errors="ignore") if l.startswith("[EDMD3")]
    print(f"acceptance 6: the repeat's CSV outputs byte-identical: {len(same)} of {len(fs)} ({', '.join(fs)}); its [EDMD3...] log lines identical "
          f"apart from the wall-time fields (run_s, engine_s, driver_share): {'yes' if [mask(x) for x in l0] == [mask(x) for x in l1] else '**NO**'}")
    # ------------------------------------------------------------ acc4
    print("\n## 7. Acceptance 4: the gate's cases under gen3 (plain runs, contact audit); the validator cadence\n")
    print("| case | records / builds | clean | validator every [steps] | events | max contact gap [px] disk-disk / walls / divider | t_end [sigma] | driver share |\n|---|---|---|---|---|---|---|---|")
    for case in ("free_M50", "free_M500", "free_M1500", "free_M2000", "dense_M50", "dense_M2000", "ctrl_min", "afix"):
        k = rec(os.path.join(O, "acc4", case, "run.log"))
        if not k: print(f"| {case} | **no run record** | | | | | | |"); continue
        c = k["contact"]; cs = f"{c[0]:.2e} / {c[1]:.2e} / {c[2]:.2e}" if c else "no contact line"
        print(f"| {case} | {k['n_records']} / {k['n_builds']} | {k['clean']} | {k['validator_every']} | {events(k)} | {cs} | {float(k['t_end']) / 24:.1f} | {k['driver_share']} |")
    print("\nThe driver's validator (experiment_validation.c) runs every --validator-every steps; under gen3 the default is 60 steps of "
          "1/60 sigma-time = once per sigma-time (G3_VALIDATOR_EVERY), under gen2 every step as before; between validator steps the "
          "Left/Right counts are carried over (the trace documentation states it, stage A2).")
    # ------------------------------------------------------------ acc5
    print("\n## 8. Acceptance 5, information: gen2 and gen3 on one fluid cell (N = 400, eta = pi/8, H = L0 = 20, M = 300; 8 seeds each, "
          "the same seeds; reduce_B.py's argmax estimator)\n")
    r = {e: pd.read_csv(os.path.join(O, "acc5", e, "cell", "m_300", "red_nu.csv")).sort_values("run") for e in ("gen2", "gen3")}
    a_, b_ = r["gen3"]["nu"].to_numpy(float), r["gen2"]["nu"].to_numpy(float)
    s = math.hypot(a_.std(ddof=1) / math.sqrt(len(a_)), b_.std(ddof=1) / math.sqrt(len(b_)))
    print("| engine | n | mean nu | SE | seeds equal |\n|---|---|---|---|---|")
    for e in ("gen2", "gen3"):
        x = r[e]["nu"].to_numpy(float); print(f"| {e} | {len(x)} | {x.mean():.6f} | {x.std(ddof=1) / math.sqrt(len(x)):.2e} | "
                                              f"{'yes' if list(r['gen2']['seed']) == list(r['gen3']['seed']) else 'NO'} |")
    print(f"\ngen3 - gen2 = {a_.mean() - b_.mean():+.2e} ({100 * (a_.mean() - b_.mean()) / b_.mean():+.2f} %), z = {(a_.mean() - b_.mean()) / s:+.2f} "
          "(information; Test G is the test)")
    outs = [("divider trace (dynamic)", glob.glob(os.path.join(O, "acc5", "gen3", "cell", "m_300", "wall_x_positions_*_run*.csv"))),
            ("impulses on the held divider and the walls (static: the event log)", glob.glob(os.path.join(O, "acc4", "afix", "ev_9700.csv"))),
            ("psi6(t)", glob.glob(os.path.join(O, "acc5", "gen3", "cell", "m_300", "r*", "psi6_t_*.csv")))]
    print("\noutputs both methods need, under gen3: " + "; ".join(f"{n}: {len(f)} files" for n, f in outs))
    ev = pd.read_csv(os.path.join(O, "acc4", "afix", "ev_9700.csv"))
    print("event-log rows by kind (acc4 afix, gen3): " + ", ".join(f"{k_} {v_}" for k_, v_ in ev["kind"].value_counts().sort_index().items()))
    # ------------------------------------------------------------ a4
    print("\n## 9. Stage A4: the driver's share of the run time (one run at a time, nothing else of mine running)\n")
    print("| loop | N | eta | record [sigma] | events | run [s] | engine [s] | driver share | events per s (run) |\n|---|---|---|---|---|---|---|---|---|")
    for loop in ("sos", "et"):
        for N in (400, 1600):
            for lab, eta in (("pi8", "pi/8"), ("070", "0.70")):
                lg = os.path.join(O, "a4", f"{loop}_N{N}_{lab}", "run.log")
                k = rec(lg) if os.path.exists(lg) else None
                if not k: print(f"| {loop} | {N} | {eta} | not run | | | | | |"); continue
                rs = float(k["run_s"])
                print(f"| {'speed of sound' if loop == 'sos' else 'energy transfer'} | {N} | {eta} | {float(k['t_end']) / 24:.0f} | {events(k)} | {rs:.1f} | "
                      f"{float(k['engine_s']):.1f} | {float(k['driver_share']):.3f} | {events(k) / rs:.3g} |")


if __name__ == "__main__":
    main()
