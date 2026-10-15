#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.18): the tables of stage E2, printed from the evidence folder of stageE_evidence.sh.
  1. the frozen binaries (78ff48d: default and long-double) against the build record
  2. E2a: the default build's M1 and M2 harness outputs against the committed ones (byte for byte)
  3. E2b: the long-double harness outputs against the default build's (byte for byte; on arm64 long double = double)
  4. E2c: the long-double driver against the default driver on two harness-cell replays (event hash, events, MATCH)
  5. rule 4: gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg; default and --engine=gen2)
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m5_261009/stageE_tables.py --bin-dir <frozen binaries> --out <evidence dir> > stageE_tables_output.txt
"""
import argparse, hashlib, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REF = (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
       ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
       ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt"))


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def same(p, q): return os.path.exists(p) and os.path.exists(q) and open(p, "rb").read() == open(q, "rb").read()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# Stage E2 (M5) evidence tables, printed by experiments_gen3_m5_261009/stageE_tables.py\n")
    rec = open(os.path.join(HERE, "build_record_78ff48d.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The frozen binaries\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    for n in ("00ALLINONE", "00ALLINONE_ld", "gen3_m1", "gen3_m1_ld", "gen3_m2", "gen3_m2_ld"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '**NO**'} |")
    print("\n## 2. E2a: the default build on the M1 and M2 harness outputs (the committed ones)\n\n| output | identical |\n|---|---|")
    for f, ref in REF: print(f"| {f} | {'IDENTICAL' if same(os.path.join(O, 'e2a', f), os.path.join(HS, ref)) else '**DIFFERENT or missing**'} |")
    print("\n## 3. E2b: the long-double harnesses against the default build's outputs (this Mac: long double = double)\n\n| output | identical |\n|---|---|")
    for f, _ in REF: print(f"| {f} | {'IDENTICAL' if same(os.path.join(O, 'e2b', f), os.path.join(O, 'e2a', f)) else '**DIFFERENT or missing**'} |")
    print("\n## 4. E2c: the long-double driver against the default driver, two harness cells replayed\n")
    t = open(os.path.join(O, "e2c", "replays.txt")).read()
    print("```\n" + t.strip() + "\n```")
    hs = re.findall(r"^(default|ld) (\S+): \[EDMD3-REPLAY\] \S+: hash (\w+), harness (\w+); events (\d+), harness (\d+): (\S+);", t, re.M)
    byc = {}
    for tag, cell, h, hh, e, he, m in hs: byc.setdefault(cell, {})[tag] = (h, e, m)
    for cell, d in byc.items():
        ok = d.get("default") and d.get("ld") and d["default"][:2] == d["ld"][:2] and d["default"][2] == d["ld"][2] == "MATCH"
        print(f"{cell}: default {d.get('default')}, long double {d.get('ld')}: {'IDENTICAL hash and events, both MATCH' if ok else '**NOT IDENTICAL**'}")
    print("\n## 5. Rule 4: gen2 byte identity with 7b08827 after the Makefile and 00ALLINONE.c changes (ctrl_min, ctrl_leg)\n")
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "e0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "e0", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "e0", f"report_{tag}.txt"), "w").write(r)
        rows = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        print(f"{'default build' if tag == 'default' else '--engine=gen2'}: " + "; ".join(
            f"{l.split('|')[1].strip()}: audit vs plain {l.split('|')[-3].strip()}, plain vs ref {l.split('|')[-2].strip()}" for l in rows))


if __name__ == "__main__":
    main()
