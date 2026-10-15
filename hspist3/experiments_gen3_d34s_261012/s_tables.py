#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 34, the rest; decision 38 item 4; 261012 sec. 4.7.67): the acceptance of the gen3 window
fitted to the screen (a DISPLAY change), printed from s_evidence.py's folder and s_build.sh's bin folder. Decision 34 item 4: rule 4
identical; the M1/M2 harness outputs identical; headless gen3 unchanged (P4a start A seed 1, the replay cell); the display does not steer
(gen3 with the window at two window sizes = the headless event hash; the output files byte-identical); 0 new warnings; screenshots at
N = 400, 5000 and 8000. With decision 38 item 4: no exit on a scene that cannot be built.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d34s_261012/s_tables.py --bin-dir <bin dir> --out <evidence dir>
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
P4_HASH, RP_HASH, ST_HASH_D19 = "bd25077c5f1c9c9f", "694fdf127f35d4b4", "daea5de1378282dc"
RECORDS = ("00_COMMAND.md", "01_PLOT_COMMANDS.md", "00_COMMAND_AND_PLOT_COMMAND.md")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def same(a, b): return os.path.exists(a) and os.path.exists(b) and filecmp.cmp(a, b, shallow=False)
def log(d): return open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
def code(d): return open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"


def health(d):
    h = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) .*?body_findings=(\d+) .* hash=([0-9a-f]+) t_end=(\S+)", log(d), re.M)
    return {"clean": h[-1][0], "body": h[-1][1], "hash": h[-1][2], "t_end": float(h[-1][3])} if h else {}


def warnings(p):
    return Counter(re.sub(r"^.*?:\d+:\d+: ", "", l.strip()) for l in open(p, errors="ignore") if "warning:" in l)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); ok = {}
    print("# Decision 34, the rest (the gen3 window fitted to the screen, a DISPLAY change): the acceptance, printed by "
          "experiments_gen3_d34s_261012/s_tables.py\n")
    rec = open(os.path.join(B, "build_record.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The binaries of the committed tree\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '(not in the record)'} |")
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
    print("\n## 3. The M1 and M2 harness outputs against the committed ones\n\n| output | identical |\n|---|---|")
    ok1 = True
    for f, ref in REF:
        e = same(os.path.join(O, "h1", f), os.path.join(HS, ref)); ok1 &= e
        print(f"| {f} | {'IDENTICAL' if e else '**DIFFERENT or missing**'} |")
    ok["M1/M2 harness outputs identical"] = ok1
    print("\n## 4. Headless gen3 unchanged\n")
    hp = health(os.path.join(O, "p4a_A1"))
    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(P4A1, "*")) if os.path.isfile(p)
                   and not p.endswith(("run.log", "command.txt", "exit_code.txt")) and os.path.basename(p) not in RECORDS)
    eq = {f: same(os.path.join(O, "p4a_A1", f), os.path.join(P4A1, f)) for f in files}
    print(f"P4a start A seed 1, P4's command with the new binary: clean={hp.get('clean')}, event hash {hp.get('hash')} (P4: {P4_HASH}): "
          f"{'EQUAL' if hp.get('hash') == P4_HASH else '**DIFFERENT**'}; its data files against P4's: {sum(eq.values())} of {len(eq)} byte-identical "
          f"({', '.join(files)})" + (f" (differ: {[f for f, v in eq.items() if not v]})" if not all(eq.values()) else ""))
    t = log(os.path.join(O, "replay"))
    rl = re.findall(r"^\[EDMD3-REPLAY\] (.*)$", t, re.M); rc = re.findall(r"^\[EDMD3-HEALTH\] .*?: (clean=[01])", t, re.M)
    rp_ok = bool(rl) and RP_HASH in rl[0] and ": MATCH;" in rl[0] and rc[:1] == ["clean=1"]
    print(f"stage A's replay cell m2_driven_pi8, variant a: {rc[:1]} {rl[0] if rl else '**no replay line**'} -- recorded hash {RP_HASH}: "
          f"{'EQUAL, MATCH' if rp_ok else '**NOT EQUAL or no MATCH**'}")
    ok["headless gen3 unchanged (P4a A seed 1; the replay cell)"] = hp.get("hash") == P4_HASH and all(eq.values()) and rp_ok
    print("\n## 5. The display does not steer (decision 19's ST scene: N = 400, eta 0.31, the right piston stepping in, 40 sigma-time; "
          "gen3; SDL's dummy video driver, the software renderer, no input)\n")
    print("| run | window | clean | event hash | = headless | = decision 19's ST hash | output files = headless |\n|---|---|---|---|---|---|---|")
    hd = os.path.join(O, "steer", "headless"); hh = health(hd)
    st_ok = bool(hh) and hh["clean"] == "1" and hh["hash"] == ST_HASH_D19
    for mode, wdesc in (("headless", "-"), ("window_A", "1024 x 768 points, UI scale 1"), ("window_B", "4112 x 2462 device pixels, UI scale 2 (2056 x 1231 points)")):
        d = os.path.join(O, "steer", mode); h = health(d)
        outs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "*.csv")) + glob.glob(os.path.join(d, "*.json")))
        eqf = {f: same(os.path.join(d, f), os.path.join(hd, f)) for f in outs}
        if mode != "headless": st_ok &= bool(h) and h["clean"] == "1" and h["hash"] == hh.get("hash") and bool(eqf) and all(eqf.values())
        print(f"| {mode} | {wdesc} | {h.get('clean')} | {h.get('hash')} | {'-' if mode == 'headless' else ('yes' if h.get('hash') == hh.get('hash') else '**NO**')} | "
              f"{'yes' if h.get('hash') == ST_HASH_D19 else '**NO**'} | {'-' if mode == 'headless' else f'{sum(eqf.values())} of {len(eqf)}'} |")
    ok["the display does not steer (gen3, two window sizes = the headless event hash, output files byte-identical)"] = st_ok
    print("\n## 6. Screenshots (eta 0.31, decision 19's display pair scaled; the divider free from the first step; the frame after step 60)\n")
    print("| N | window | UI scale | screenshot | gen3 clean | frames per second |\n|---|---|---|---|---|---|")
    ss_ok = True
    for tag, N, win, us in (("N400_screen", 400, "2056 x 1231 points (Chris's display)", 1), ("N5000_screen", 5000, "2056 x 1231 points", 1),
                            ("N8000_screen", 8000, "2056 x 1231 points", 1), ("N5000_retina", 5000, "4112 x 2462 device pixels", 2),
                            ("N5000_small", 5000, "1440 x 900 points (re-fitted)", 1)):
        d = os.path.join(O, "screens", tag); png = os.path.join(O, "screenshots", f"SS_{tag}.png"); h = health(d)
        g = re.findall(r"^\[GUI-STATS\] .*? fps=(\S+)", log(d), re.M)
        ss_ok &= os.path.exists(png) and h.get("clean") == "1" and code(d) == "0"
        print(f"| {N} | {win} | {us} | {os.path.relpath(png, HS) if os.path.exists(png) else '**missing**'} | {h.get('clean')} | {g[-1] if g else '-'} |")
    ok["screenshots at N = 400, 5000 and 8000 (and N = 5000 at UI scale 2 and in a smaller window)"] = ss_ok
    print("\n## 7. No exit on a scene that cannot be built\n")
    da = os.path.join(O, "noexit", "apply_refused"); ta = log(da)
    ra = re.findall(r"^\[GUI\] APPLY refused: .*$", ta, re.M)
    print(f"- (a) the F3 APPLY of a scene the seeder cannot place: `{ra[0] if ra else '**no refusal line**'}`; rebuilt: "
          f"{'**yes**' if 'Scene reset queued' in ta else 'no'}; 'Could only place' in the log: {'**yes**' if 'Could only place' in ta else 'no'}; exit code {code(da)}")
    db = os.path.join(O, "noexit", "gen3_refused"); tb = log(db)
    rb = re.findall(r"^\[GUI-STOP\] .*$", tb, re.M); png = os.path.join(O, "screenshots", "NX_gen3_refused_pause.png")
    print(f"- (b) a scene gen3 cannot build: `{rb[0] if rb else '**no [GUI-STOP] line**'}`; exit code {code(db)}; screenshot "
          f"{os.path.relpath(png, HS) if os.path.exists(png) else '**missing**'}")
    dc = os.path.join(O, "noexit", "gen3_refused_old_bin_8bc9ac9"); tc = log(dc)
    rcx = re.findall(r"^STOP \(--engine=gen3\): .*$", tc, re.M)
    print(f"- (c) information, (b) with the addition's binary bin_8bc9ac9: `{rcx[0] if rcx else 'no STOP line'}`, exit code {code(dc)}")
    ok["no exit: the APPLY refused with a note (scene not rebuilt), the gen3 refusal pauses the window"] = (
        bool(ra) and "Scene reset queued" not in ta and "Could only place" not in ta and code(da) == "0"
        and bool(rb) and code(db) == "0" and os.path.exists(png))
    new, old = warnings(os.path.join(B, "00ALLINONE.Wall.txt")), warnings(os.path.join(B, "00ALLINONE.Wall_parent_f72d764.txt"))
    added = new - old
    print(f"\n## 8. Warnings (build_clean.sh's -Wall -Wextra pass of the driver's sources; compared without file and line)\n\n"
          f"this build {sum(new.values())}, the parent f72d764 {sum(old.values())}; new: {sum(added.values())}" + (f": {sorted(added)}" if added else ""))
    ok["0 new warnings"] = not added
    print("\n## Acceptance\n\n| item | result |\n|---|---|")
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else '**FAIL**'} |")
    print(f"\nDECISION 34 (THE REST) ACCEPTANCE: {'PASS' if all(ok.values()) else 'FAIL'}")


if __name__ == "__main__":
    main()
