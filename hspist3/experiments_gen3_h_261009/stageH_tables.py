#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.21; stage H): the tables of stage H, printed from the evidence folder of stageH_evidence.py.
  1. the frozen binaries against the build record
  2. H0 rule 4: gen2 byte identity with 7b08827 (ctrl_min, ctrl_leg; the default build and --engine=gen2)
  3. H1 the M1 and M2 harness outputs against the committed ones
  4. H2 the engine test: the verdict of both builds; the long-double output against the default's (on arm64 long double = double)
  5. H3 the driver: every finished run against the uninterrupted run U of its trajectory: the trace, the psi6(t) file,
     speed_of_sound_psi6.csv, run_params.json, the run record ([EDMD3-HEALTH]) without its three timing fields, and the event
     hash; the checkpoint and restart lines (time, events, hash at the checkpoint)
  6. H4 the refusals: exit code and the STOP line
The acceptance (sec. 4.7.12, stage H): a restarted run is byte-identical to the uninterrupted one (traces, psi6(t), event hash).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_h_261009/stageH_tables.py --bin-dir <frozen binaries> --record <build record> --out <evidence dir>
"""
import argparse, glob, hashlib, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REF = (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
       ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
       ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt"))


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, "rb").read() if os.path.exists(p) else None
def same(p, q): a, b = rd(p), rd(q); return a is not None and b is not None and a == b


def record(d):
    """the run record without its timing fields, the hash, and the checkpoint / restart lines"""
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    h = re.findall(r"^(\[EDMD3-HEALTH\] .*?) run_s=\S+ engine_s=\S+ driver_share=\S+$", lg, re.M)
    hs = re.findall(r" hash=([0-9a-f]{16}) ", h[0]) if h else []
    ck = re.findall(r"^\[EDMD3-(CHECKPOINT\] wrote|RESTART\] read) \S+: (\w+ step \d+), t = (\S+) units, events (\d+), hash ([0-9a-f]{16})", lg, re.M)
    stop = re.findall(r"^STOP \(--engine=gen3\): (.*)$|^(STOP: .*)$", lg, re.M)
    return (h[0] if h else None), (hs[0] if hs else None), ck, [x or y for x, y in stop]


def files(d):
    out = {}
    for pat in ("wall_x_positions_*_run0.csv", "psi6_t_*_run0.csv", "speed_of_sound_psi6.csv", "run_params.json"):
        g = glob.glob(os.path.join(d, pat))
        out[pat] = g[0] if g else None
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--record", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# Stage H (checkpoint and restart) evidence tables, printed by experiments_gen3_h_261009/stageH_tables.py\n")
    rec = open(a.record).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The frozen binaries\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2", "gen3_checkpoint_test", "gen3_checkpoint_test_ld"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '**NO**'} |")
    print("\n## 2. H0, rule 4: gen2 byte identity with 7b08827 after the 00ALLINONE.c change (ctrl_min, ctrl_leg)\n")
    ok0 = True
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "h0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "h0", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "h0", f"report_{tag}.txt"), "w").write(r)
        rows = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        ok0 &= len(rows) == 2 and all(l.split("|")[-3].strip() == "IDENTICAL" and l.split("|")[-2].strip() == "IDENTICAL" for l in rows)
        print(f"{'default build' if tag == 'default' else '--engine=gen2'}: " + "; ".join(
            f"{l.split('|')[1].strip()}: audit vs plain {l.split('|')[-3].strip()}, plain vs ref {l.split('|')[-2].strip()}" for l in rows))
    print(f"\nrule 4: {'IDENTICAL (both cases, both builds)' if ok0 else '**NOT IDENTICAL**'}")
    print("\n## 3. H1: the M1 and M2 harness outputs (the committed ones)\n\n| output | identical |\n|---|---|")
    ok1 = True
    for f, ref in REF:
        e = same(os.path.join(O, 'h1', f), os.path.join(HS, ref)); ok1 &= e
        print(f"| {f} | {'IDENTICAL' if e else '**DIFFERENT or missing**'} |")
    print("\n## 4. H2: the engine test (edmd_core/tests/gen3_checkpoint_test.c)\n")
    ok2 = True
    for tag in ("default", "ld"):
        t = open(os.path.join(O, "h2", tag, "checkpoint_test_output.txt")).read()
        v = re.findall(r"^VERDICT: .*$", t, re.M)
        ok2 &= bool(v) and v[0].startswith("VERDICT: PASS")
        print(f"{tag}: {v[0] if v else '**no verdict line**'}")
    print(f"long-double output = default output, byte for byte: "
          f"{'yes' if same(os.path.join(O, 'h2', 'ld', 'checkpoint_test_output.txt'), os.path.join(O, 'h2', 'default', 'checkpoint_test_output.txt')) else 'no (see the build line)'}")
    print("\n## 5. H3: the driver -- every finished run against the uninterrupted run U of its trajectory\n")
    print("| trajectory | run | exit | checkpoint / restart (t [units], events, hash at that moment) | trace = U | psi6(t) = U | "
          "speed_of_sound_psi6.csv = U | run_params.json = U | run record (without timing) = U | event hash at the end |\n|---|---|---|---|---|---|---|---|---|---|")
    acc_R, acc_C, nR, nC = True, True, 0, 0
    for cell in ("A", "B", "C"):
        base = os.path.join(O, "h3", cell)
        if not os.path.isdir(base): continue
        U = os.path.join(base, "U"); fu = files(U); hu, hashu, _, _ = record(U)
        runs = sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))
        runs = ["U"] + [d for d in runs if d != "U"]
        for r in runs:
            d = os.path.join(base, r)
            ex = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
            h, hs, ck, stop = record(d)
            cks = "; ".join(f"{'wrote' if w.startswith('CHECKPOINT') else 'read'} {s} (t {t}, {e} events, {hh})" for w, s, t, e, hh in ck) or "-"
            if r.startswith("R1") or r == "R2b_chain":   # stopped right after writing: no finished outputs
                print(f"| {cell} | {r} | {ex} | {cks} | (stopped) | | | | | |"); continue
            f = files(d)
            eq = [same(f[k], fu[k]) for k in ("wall_x_positions_*_run0.csv", "psi6_t_*_run0.csv", "speed_of_sound_psi6.csv", "run_params.json")]
            eqh = h is not None and h == hu
            cells_ = ["IDENTICAL" if e else "**DIFFERENT**" for e in eq]
            if r != "U":
                full = all(eq) and eqh and hs == hashu and ex == "0"
                if r.startswith("C_"): acc_C &= full; nC += 1
                else: acc_R &= full; nR += 1
            print(f"| {cell} | {r} | {ex} | {cks} | " + (" | ".join(cells_) if r != "U" else "(reference) | | | ") +
                  f" | {'(reference)' if r == 'U' else ('IDENTICAL' if eqh else '**DIFFERENT**')} | {hs} |")
    print(f"\nH3: every restarted run (R2, R3; {nR}) byte-identical to U in the trace, psi6(t), the psi6 summary, run_params.json, the run "
          f"record and the event hash: {'YES' if acc_R and nR else '**NO**'}; every run that wrote a checkpoint and went on (C; {nC}): "
          f"{'YES' if acc_C and nC else '**NO**'}")
    print("\n## 6. H4: the refusals (each must stop with exit 2 and say why)\n\n| case | exit | message |\n|---|---|---|")
    hb = os.path.join(O, "h4"); okr = True
    for k in range(64):
        nf = os.path.join(hb, f"case{k}.name")
        if not os.path.exists(nf): break
        d = os.path.join(hb, f"case{k}")
        ex = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
        lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
        m = re.findall(r"^(STOP.*)$", lg, re.M)
        okr &= ex == "2" and bool(m)
        print(f"| {open(nf).read().strip()} | {ex} | {m[-1] if m else '**no STOP line**'} |")
    print(f"\nH4: {'every case refused with exit 2 and a reason' if okr else '**NOT every case refused**'}")
    print(f"\nACCEPTANCE (stage H): {'PASS' if acc_R and nR and acc_C and ok0 and ok1 and ok2 and okr else '**NOT MET**'} -- restarted runs "
          f"byte-identical to the uninterrupted ones (traces, psi6(t), event hash); writing a checkpoint does not steer; the engine test "
          f"passes in both builds; rule 4 identical; the M1 and M2 harness outputs unchanged; refusals refuse")


if __name__ == "__main__":
    main()
