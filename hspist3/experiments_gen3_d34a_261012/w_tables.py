#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 34, addition; 261012 sec. 4.7.73): the acceptance of the energy-transfer scene in the gen3
window (a DISPLAY change), printed from w_evidence.py's folder and w_build.sh's bin folder. Acceptance: decision 34's item 4 as decision
19 item 2 -- rule 4 identical; the M1/M2 harness outputs identical; headless gen3 unchanged (P4a A seed 1, the replay cell); the
display does not steer (gen3 with the window at two window sizes = the headless event hash, the output files byte-identical); 0 new
warnings -- and the addition's item E: a spring scene at two values of k (KE + spring energy conserved to rounding over 200 sigma-time
after the release, gen3's ledger line; the same scene headless = the window's event hash); a divider that reaches its limit pauses the
window (screenshot). Screenshots at N = 5000 and 8000 belong to decision 34's drawing scale (the rest of decision 34), not to this part.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d34a_261012/w_tables.py --bin-dir <bin dir> --out <evidence dir>
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
RECORDS = ("00_COMMAND.md", "01_PLOT_COMMANDS.md", "00_COMMAND_AND_PLOT_COMMAND.md")   # provenance: the run's folder, binary, time, command
SPRING_KS = ("11.52", "115.2")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def same(a, b): return os.path.exists(a) and os.path.exists(b) and filecmp.cmp(a, b, shallow=False)
def log(d): return open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
def code(d): return open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"


def health(d):
    """the last run record: clean, event hash, t_end [internal units], body_findings"""
    h = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) .*?body_findings=(\d+) .* hash=([0-9a-f]+) t_end=(\S+)", log(d), re.M)
    return {"clean": h[-1][0], "body": h[-1][1], "hash": h[-1][2], "t_end": float(h[-1][3])} if h else {}


def ledger(d):
    """the last [EDMD3-LEDGER] line's energy part: residual, scale, ratio, E0, E, W"""
    m = re.findall(r"^\[EDMD3-LEDGER\] .*energy: residual (\S+), scale (\S+), ratio (\S+) \(E0 (\S+), E (\S+), W (\S+)\)", log(d), re.M)
    return dict(zip(("residual", "scale", "ratio", "E0", "E", "W"), (float(x) for x in m[-1]))) if m else {}


def data_files(d):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "*")) if os.path.isfile(p)
                  and os.path.basename(p) not in RECORDS + ("run.log", "command.txt", "exit_code.txt", "Roboto-Regular.ttf"))


def warnings(p):
    return Counter(re.sub(r"^.*?:\d+:\d+: ", "", l.strip()) for l in open(p, errors="ignore") if "warning:" in l)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); ok = {}
    print("# Decision 34, addition (the energy-transfer scene in the gen3 window, a DISPLAY change): the acceptance, printed by "
          "experiments_gen3_d34a_261012/w_tables.py\n")
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
    hp = health(os.path.join(O, "p4a_A1"))
    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(P4A1, "*")) if os.path.isfile(p)
                   and not p.endswith(("run.log", "command.txt", "exit_code.txt")) and os.path.basename(p) not in RECORDS)
    eq = {f: same(os.path.join(O, "p4a_A1", f), os.path.join(P4A1, f)) for f in files}
    print(f"P4a start A seed 1, P4's command with the new binary: clean={hp.get('clean')}, event hash {hp.get('hash')} (P4: {P4_HASH}): "
          f"{'EQUAL' if hp.get('hash') == P4_HASH else '**DIFFERENT**'}; its data files against P4's: {sum(eq.values())} of {len(eq)} byte-identical "
          f"({', '.join(files)})" + (f" (differ: {[f for f, v in eq.items() if not v]})" if not all(eq.values()) else "")
          + f"; the command records {', '.join(RECORDS)} name the run's folder, binary and time (provenance; not compared)")
    t = log(os.path.join(O, "replay"))
    rl = re.findall(r"^\[EDMD3-REPLAY\] (.*)$", t, re.M); rc = re.findall(r"^\[EDMD3-HEALTH\] .*?: (clean=[01])", t, re.M)
    rp_ok = bool(rl) and RP_HASH in rl[0] and ": MATCH;" in rl[0] and rc[:1] == ["clean=1"]
    print(f"stage A's replay cell m2_driven_pi8, variant a (readers at every target, one stop per interval): {rc[:1]} {rl[0] if rl else '**no replay line**'}"
          f" -- recorded hash {RP_HASH}: {'EQUAL, MATCH' if rp_ok else '**NOT EQUAL or no MATCH**'}")
    ok["headless gen3 unchanged (P4a A seed 1; the replay cell)"] = hp.get("hash") == P4_HASH and all(eq.values()) and rp_ok
    # ---- the display does not steer, two window sizes
    print("\n## 5. The display does not steer (decision 19's ST scene: N = 400, eta 0.31, the right piston stepping in, 40 sigma-time; "
          "gen3; the window: SDL's dummy video driver, the software renderer, no input)\n")
    print("| run | window | clean | event hash | = headless | = decision 19's ST hash | output files = headless |\n|---|---|---|---|---|---|---|")
    hd = os.path.join(O, "steer", "headless"); hh = health(hd)
    st_ok = bool(hh) and hh["clean"] == "1" and hh["hash"] == ST_HASH_D19
    for mode, wdesc in (("headless", "-"), ("window_A", "1024 x 768 points, UI scale 1"), ("window_B", "4112 x 2462 device pixels, UI scale 2 (2056 x 1231 points)")):
        d = os.path.join(O, "steer", mode); h = health(d)
        outs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "*.csv")) + glob.glob(os.path.join(d, "*.json")))
        eqf = {f: same(os.path.join(d, f), os.path.join(hd, f)) for f in outs}
        if mode != "headless": st_ok &= bool(h) and h["clean"] == "1" and h["hash"] == hh.get("hash") and bool(eqf) and all(eqf.values())
        print(f"| {mode} | {wdesc} | {h.get('clean')} | {h.get('hash')} | {'-' if mode == 'headless' else ('yes' if h.get('hash') == hh.get('hash') else '**NO**')} | "
              f"{'yes' if h.get('hash') == ST_HASH_D19 else '**NO**'} | {'-' if mode == 'headless' else f'{sum(eqf.values())} of {len(eqf)}'}"
              + (f" (differ: {[f for f, v in eqf.items() if not v]})" if mode != "headless" and not all(eqf.values()) else "") + " |")
    ok["the display does not steer (gen3, two window sizes = the headless event hash, output files byte-identical)"] = st_ok
    # ---- item E: the spring scene at two k
    print("\n## 6. Item E: the energy-transfer scene with a spring, 200 sigma-time after the release (N = 400 on the right, r = 0.2159 sigma, "
          "eta 0.116; spring wall 0, L0 = 25.34 sigma, M = 50; free from the first step)\n")
    print("gen3's ledger line: E = KE (disks and divider) + the spring energy, W = the work of driven bodies and API changes; "
          "(E - E0) - W against its rounding scale (u x the sum of the rounded terms); ratio <= 1: conserved to rounding.\n")
    print("| k [kT/sigma^2] | run | clean | event hash | sigma-time (t_end/24) | E0 | E | W | (E - E0) - W | scale | ratio | [ENERGY] first and last: E_spring, E_tot | data files = headless |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    sp_ok = True
    for k in SPRING_KS:
        dh = os.path.join(O, "spring", f"k{k}_headless")
        for mode in ("headless", "window_A"):
            d = os.path.join(O, "spring", f"k{k}_{mode}"); h = health(d); L = ledger(d)
            en = re.findall(r"^\[ENERGY\] step=(\d+) .*?E_spring=(\S+)  E_tot=(\S+)$", log(d), re.M)
            fl = data_files(d); eqf = {f: same(os.path.join(d, f), os.path.join(dh, f)) for f in fl}
            row_ok = bool(h) and h["clean"] == "1" and bool(L) and L["ratio"] <= 1.0 and h["t_end"] / 24.0 >= 200.0
            if mode != "headless": row_ok &= h.get("hash") == health(dh).get("hash") and bool(eqf) and all(eqf.values())
            sp_ok &= row_ok
            print(f"| {k} | {mode} | {h.get('clean')} | {h.get('hash')} | {h.get('t_end', float('nan')) / 24.0:.3f} | {L.get('E0')} | {L.get('E')} | {L.get('W')} | "
                  f"{L.get('residual')} | {L.get('scale')} | {L.get('ratio')} | "
                  + (f"step {en[0][0]}: {en[0][1]}, {en[0][2]}; step {en[-1][0]}: {en[-1][1]}, {en[-1][2]}" if en else "-")
                  + f" | {'-' if mode == 'headless' else f'{sum(eqf.values())} of {len(eqf)}'} |")
    ok["item E: the spring scene at two k -- ledger to rounding over 200 sigma-time, window = headless event hash"] = sp_ok
    # ---- item B through the F3 panel
    print("\n## 7. Item B through the F3 panel: the scene rebuilt by APPLY (--gui-keys=1:apply-gen3), k = 11.52 from the k row, window B\n")
    d = os.path.join(O, "spring_f3", "k11.52_window_B"); t = log(d); h = health(d); L = ledger(d)
    q = re.findall(r"^Scene reset queued: .*$", t, re.M); s = re.findall(r"^.*spring: k = .*$", t, re.M)
    print(f"- the panel's line: `{q[-1] if q else '**none**'}`")
    print(f"- the engine's spring after the rebuild: `{s[-1].strip() if s else '**none**'}`")
    print(f"- the run record: clean={h.get('clean')}, event hash {h.get('hash')}, {h.get('t_end', float('nan')) / 24.0:.3f} sigma-time; the ledger: "
          f"(E - E0) - W = {L.get('residual')}, scale {L.get('scale')}, ratio {L.get('ratio')} (E0 {L.get('E0')}, E {L.get('E')}, W {L.get('W')})")
    png = os.path.join(O, "screenshots", "SF_k11.52_window_B_step3000.png")
    print(f"- the screenshot after step 3000: {os.path.relpath(png, HS) if os.path.exists(png) else '**missing**'}")
    ok["item B: the k row reaches gen3 through APPLY (springk=11.52, the engine's k = 0.02 /px^2, ledger to rounding)"] = (
        bool(q) and "springk=11.52" in q[-1] and bool(s) and "k = 0.02 /px^2" in s[-1] and h.get("clean") == "1" and bool(L) and L["ratio"] <= 1.0
        and os.path.exists(png))
    # ---- item C: the pause at the limit
    print("\n## 8. Item C: a divider that reaches its limit pauses the window (no spring: Chris's case; window B)\n")
    d = os.path.join(O, "limit", "new_window_B"); t = log(d)
    g = re.findall(r"^\[GUI-LIMIT\] step .*$", t, re.M)
    lr = re.findall(r"^\[EDMD3-HEALTH\] scene, a divider at the end of its range: clean=(\S+) .*?body_findings=(\d+) .* hash=([0-9a-f]+)", t, re.M)
    png = os.path.join(O, "screenshots", "LM_window_B_pause.png")
    print(f"- `{g[0] if g else '**no [GUI-LIMIT] line**'}`")
    print(f"- the run record at the pause: " + (f"clean={lr[0][0]}, body_findings={lr[0][1]}, event hash {lr[0][2]}" if lr else "**missing**")
          + f"; exit code {code(d)}; 'Invalid wall_x' in the log: {'**yes**' if 'Invalid wall_x' in t else 'no'}")
    print(f"- the screenshot of the pause: {os.path.relpath(png, HS) if os.path.exists(png) else '**missing**'}")
    dold = os.path.join(O, "limit", "old_bin_7f6adca_window"); to = log(dold)
    iw = re.findall(r"^.*Invalid wall_x = \S+.*$", to, re.M)
    print(f"- information, the same scene with decision 19's binary bin_7f6adca: `{iw[0].strip() if iw else 'no Invalid wall_x line'}`, exit code {code(dold)}")
    ok["item C: the window pauses at the limit (no exit; the run record clean)"] = (bool(g) and bool(lr) and lr[0][0] == "1" and lr[0][1] == "0"
                                                                                   and code(d) == "0" and "Invalid wall_x" not in t and os.path.exists(png))
    # ---- warnings
    new, old = warnings(os.path.join(B, "00ALLINONE.Wall.txt")), warnings(os.path.join(B, "00ALLINONE.Wall_parent_f72d764.txt"))
    added = new - old
    print(f"\n## 9. Warnings (build_clean.sh's -Wall -Wextra pass of the driver's sources; compared without file and line)\n\n"
          f"this build {sum(new.values())}, the parent f72d764 {sum(old.values())}; new: {sum(added.values())}" + (f": {sorted(added)}" if added else ""))
    for n in ("gen3_m1", "gen3_m2"):
        print(f"{n}: {sum(1 for l in open(os.path.join(B, f'{n}.warnings.txt'), errors='ignore') if 'warning:' in l)} warnings")
    ok["0 new warnings"] = not added
    print("\n## Acceptance\n\n| item | result |\n|---|---|")
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else '**FAIL**'} |")
    print(f"\nDECISION 34 ADDITION ACCEPTANCE: {'PASS' if all(ok.values()) else 'FAIL'}")


if __name__ == "__main__":
    main()
