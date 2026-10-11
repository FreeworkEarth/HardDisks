#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 19, item 2; 261012 sec. 4.7.50): the acceptance of gen3 in the window (a provisional
DISPLAY feature), printed from d19_evidence.py's folder and d19_build.sh's bin folder. Acceptance (decision 19, verbatim in sec. 4.7.44):
rule 4 identical (default and --engine=gen2, ctrl_min, ctrl_leg); M1/M2 harness outputs identical; headless gen3 unchanged (P4a A seed 1
gives bd25077c5f1c9c9f, and one stage-A replay cell gives its recorded hash); the display does not steer (the same seed and protocol
with the window, SDL dummy video driver, no input, gives the headless event hash); a table for N = 5000 and 8000 at eta 0.31 (frames per
second and sigma-time per wall second, gen3 next to gen2); 0 new warnings.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d19_261010/d19_tables.py --bin-dir <bin dir> --out <evidence dir>
"""
import argparse, filecmp, glob, hashlib, os, re, shutil, subprocess, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REF = (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
       ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
       ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt"))
P4A1 = os.path.join(HS, "experiments_gen3_p4_261010", "data_06e7248", "N400_eta0.740", "A", "seed1")
P4_HASH, RP_HASH = "bd25077c5f1c9c9f", "694fdf127f35d4b4"


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def same(a, b): return os.path.exists(a) and os.path.exists(b) and filecmp.cmp(a, b, shallow=False)


def g3hash(d):
    t = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    h = re.findall(r"^\[EDMD3-HEALTH\] .*? clean=(\S+) .* hash=([0-9a-f]+)", t, re.M)
    return h[-1] if h else (None, None)


def warnings(p):
    return Counter(re.sub(r"^.*?:\d+:\d+: ", "", l.strip()) for l in open(p, errors="ignore") if "warning:" in l)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); ok = {}
    print("# Decision 19 item 2 (gen3 in the window, a provisional DISPLAY feature): the acceptance, printed by "
          "experiments_gen3_d19_261010/d19_tables.py\n")
    rec = open(os.path.join(B, "build_record.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The binaries of the committed tree\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '(not in the record)'} |")
    # ---- rule 4
    print("\n## 2. Rule 4: gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg; the default build and --engine=gen2)\n")
    ok0 = True
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "h0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "h0", tag)],
                           capture_output=True, text=True, cwd=HS).stdout
        open(os.path.join(O, "h0", f"report_{tag}.txt"), "w").write(r)
        rows = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        ok0 &= len(rows) == 2 and all(l.split("|")[-3].strip() == "IDENTICAL" and l.split("|")[-2].strip() == "IDENTICAL" for l in rows)
        print(f"{'default build' if tag == 'default' else '--engine=gen2'}: " + "; ".join(
            f"{l.split('|')[1].strip()}: audit vs plain {l.split('|')[-3].strip()}, plain vs ref {l.split('|')[-2].strip()}" for l in rows))
    ok["rule 4 identical (default and --engine=gen2, ctrl_min and ctrl_leg)"] = ok0
    # ---- M1/M2
    print("\n## 3. The M1 and M2 harness outputs against the committed ones\n\n| output | identical |\n|---|---|")
    ok1 = True
    for f, ref in REF:
        e = same(os.path.join(O, "h1", f), os.path.join(HS, ref)); ok1 &= e
        print(f"| {f} | {'IDENTICAL' if e else '**DIFFERENT or missing**'} |")
    ok["M1/M2 harness outputs identical"] = ok1
    # ---- headless gen3 unchanged
    print("\n## 4. Headless gen3 unchanged\n")
    c, h = g3hash(os.path.join(O, "p4a_A1"))
    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(P4A1, "*")) if os.path.isfile(p) and not p.endswith(("run.log", "command.txt", "exit_code.txt")))
    eq = {f: same(os.path.join(O, "p4a_A1", f), os.path.join(P4A1, f)) for f in files}
    print(f"P4a start A seed 1, P4's command with the new binary: clean={c}, event hash {h} (P4: {P4_HASH}): {'EQUAL' if h == P4_HASH else '**DIFFERENT**'}; "
          f"its files against P4's: {sum(eq.values())} of {len(eq)} byte-identical" + (f" (differ: {[f for f, v in eq.items() if not v]})" if not all(eq.values()) else ""))
    t = open(os.path.join(O, "replay", "run.log"), errors="ignore").read() if os.path.exists(os.path.join(O, "replay", "run.log")) else ""
    rl = re.findall(r"^\[EDMD3-REPLAY\] (.*)$", t, re.M); rc = re.findall(r"^\[EDMD3-HEALTH\] .*?: (clean=[01])", t, re.M)
    rp_ok = bool(rl) and RP_HASH in rl[0] and ": MATCH;" in rl[0] and rc[:1] == ["clean=1"]
    print(f"stage A's replay cell m2_driven_pi8, variant a (readers at every target, one stop per interval): {rc[:1]} {rl[0] if rl else '**no replay line**'}"
          f" -- recorded hash {RP_HASH}: {'EQUAL, MATCH' if rp_ok else '**NOT EQUAL or no MATCH**'}")
    ok["headless gen3 unchanged (P4a A seed 1; the replay cell)"] = h == P4_HASH and all(eq.values()) and rp_ok
    # ---- the display does not steer
    print("\n## 5. The display does not steer (N = 400, eta 0.31, the right piston stepping in; 40 sigma-time; the window: SDL's dummy video "
          "driver, the software renderer, no input)\n")
    print("| engine | window: clean, event hash | headless: clean, event hash | equal | the scene's output files, window = headless |\n|---|---|---|---|---|")
    st_ok = True
    for eng in ("gen3", "gen2"):
        w, hl = os.path.join(O, "steer", f"{eng}_window"), os.path.join(O, "steer", f"{eng}_headless")
        cw, hw = g3hash(w); ch, hh = g3hash(hl)
        outs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(w, "*.csv")) + glob.glob(os.path.join(w, "*.json")))
        eqf = {f: same(os.path.join(w, f), os.path.join(hl, f)) for f in outs}
        if eng == "gen3": st_ok = hw is not None and hw == hh and cw == ch == "1"
        print(f"| {eng} | {cw}, {hw} | {ch}, {hh} | {('yes' if hw == hh else '**NO**') if eng == 'gen3' else '- (gen2: no event hash)'} | "
              f"{sum(eqf.values())} of {len(eqf)}" + (f" (differ: {[f for f, v in eqf.items() if not v]})" if not all(eqf.values()) else "") + " |")
    ok["the display does not steer (gen3: window = headless event hash)"] = st_ok
    # ---- performance
    print("\n## 6. The window at N = 5000 and 8000, eta 0.31 (SDL's dummy video driver, the software renderer; --auto-release, the "
          "window's own pacing; 30 wall seconds each, one run at a time, nothing else of ours running)\n")
    print("| N | time scale | engine | frames per second | sigma-time per wall second | events per wall second | frames | steps | clean (gen3) |\n|---|---|---|---|---|---|---|---|---|")
    pf_ok = True
    for N in (5000, 8000):
        for ts in ("1", "10000"):
            for eng in ("gen2", "gen3"):
                d = os.path.join(O, "perf", f"N{N}_ts{ts}_{eng}")
                t = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
                g = re.findall(r"^\[GUI-STATS\] (.*)$", t, re.M)
                kv = dict(re.findall(r"(\w+)=(\S+)", g[-1])) if g else {}
                cl = g3hash(d)[0] if eng == "gen3" else "-"
                pf_ok &= bool(kv) and (eng == "gen2" or cl == "1")
                print(f"| {N} | {ts} | {eng} | {kv.get('fps', '**no [GUI-STATS]**')} | {kv.get('sigma_per_wall_s', '-')} | {kv.get('events_per_wall_s', '-') if eng == 'gen3' else '-'} | "
                      f"{kv.get('frames', '-')} | {kv.get('steps', '-')} | {cl} |")
    ok["the table N = 5000 and 8000 printed (every run's [GUI-STATS]; gen3 clean)"] = pf_ok
    # ---- warnings
    new, old = warnings(os.path.join(B, "00ALLINONE.Wall.txt")), warnings(os.path.join(B, "00ALLINONE.Wall_parent_f72d764.txt"))
    added = new - old
    print(f"\n## 7. Warnings (build_clean.sh's -Wall -Wextra pass of the driver's sources; compared without file and line)\n\n"
          f"this build {sum(new.values())}, the parent f72d764 {sum(old.values())}; new: {sum(added.values())}" + (f": {sorted(added)}" if added else ""))
    for n in ("gen3_m1", "gen3_m2"):
        print(f"{n}: {sum(1 for l in open(os.path.join(B, f'{n}.warnings.txt'), errors='ignore') if 'warning:' in l)} warnings")
    ok["0 new warnings"] = not added
    print("\n## Acceptance\n\n| item | result |\n|---|---|")
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else '**FAIL**'} |")
    print(f"\nDECISION 19 ITEM 2 ACCEPTANCE: {'PASS' if all(ok.values()) else 'FAIL'}")


if __name__ == "__main__":
    main()
