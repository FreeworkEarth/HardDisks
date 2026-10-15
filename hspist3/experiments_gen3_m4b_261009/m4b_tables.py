#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.30; plan-author decision 13, item 4): the tables of M4b, the compressed-fluid start, printed from
the evidence folder of m4b_evidence.py. The plan author's acceptance (N = 400 and 1600, targets 0.70 and 0.716), as tested here (stated in
this script's commit, before the run):
  1. rule 4 identical: ctrl_min and ctrl_leg, "plain vs ref" IDENTICAL against 7b08827, the default build and --engine=gen2;
  2. the M1/M2 harness outputs byte-identical to the committed ones;
  3. per cell, from the [EDMD3-M4B] line: the achieved L0 equals the requested one within 1e-9 sigma on both sides, and is recorded;
     N_L = N_R = N/2 with no disk in the other compartment; KE per compartment = N_s kT to rounding after the rescale (abs(KE/N_s kT - 1)
     <= 1e-12); the piston work equals the KE gain to the ledger's rounding scale (abs(KE gain - piston work) <= the ledger's scale_E,
     and the ledger's own residual ratio <= 1); the compression state clean=1 and the run's record clean=1; psi6 at the end of the
     compression printed; the smallest wall and divider gaps after the carry printed (> -1e-6 px, or the run would have stopped);
  4. the same seed bit-identical: "repeat" against "plain": the trace, psi6(t), the psi6 summary, run_params.json, the [EDMD3-M4B] lines
     (without the trailing wall time) and the run record (without its timing fields) identical;
  5. observation does not steer: "obs" (psi6 every 0.25, the virial blocks, the snapshots) against "plain": the event hash, the trace,
     the psi6 summary and the [EDMD3-M4B] lines identical;
  6. checkpoint and restart work through the compression: the checkpoint files of C and R1 byte-identical; R2 (restarted) and C
     identical to U in the trace, psi6(t), the psi6 summary, the run record (without timing) and the event hash; R2 repeats the
     preparation ([EDMD3-M4B] line identical to U's) and reads the checkpoint;
  7. the refusals stop with exit 2 and their message.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m4b_261009/m4b_tables.py --bin-dir <frozen binaries> --record <build record> --out <evidence dir>
"""
import argparse, glob, hashlib, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
GATE = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/cluster/resched_gate_261005"
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REF = (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
       ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
       ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt"))
CELLS = ("N1600_eta0.700", "N1600_eta0.716", "N400_eta0.700", "N400_eta0.716")


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, "rb").read() if p and os.path.exists(p) else None
def same(p, q): a, b = rd(p), rd(q); return a is not None and b is not None and a == b
def log(d): p = os.path.join(d, "run.log"); return open(p, errors="ignore").read() if os.path.exists(p) else ""
def one(d, pat): g = glob.glob(os.path.join(d, pat)); return g[0] if len(g) == 1 else None
def exitcode(d): p = os.path.join(d, "exit_code.txt"); return open(p).read().strip() if os.path.exists(p) else "?"


def m4b_lines(d):
    """the two [EDMD3-M4B] lines, the first without its trailing wall time"""
    L = re.findall(r"^\[EDMD3-M4B\] .*$", log(d), re.M)
    return [re.sub(r"; [0-9.]+ s$", "", l) for l in L]


def m4b(d):
    l = m4b_lines(d)
    if not l: return None
    t = l[0]; g = lambda pat: re.search(pat, t)
    v = dict(line=t, psi6=l[1] if len(l) > 1 else None)
    m = g(r"L0 requested (\S+) sigma, achieved left (\S+) right (\S+) sigma"); v["L0"], v["Ll"], v["Lr"] = (float(x) for x in m.groups())
    m = g(r"piston work left (\S+) right (\S+);"); v["wL"], v["wR"] = (float(x) for x in m.groups())
    m = g(r"KE left/right at the start (\S+) (\S+), before the rescale (\S+) (\S+), after (\S+) (\S+) \(after/N_s kT - 1: (\S+) (\S+)\)")
    v["ke0"], v["ke1"], v["ke2"], v["rel"] = (float(m.group(1)), float(m.group(2))), (float(m.group(3)), float(m.group(4))), \
        (float(m.group(5)), float(m.group(6))), (float(m.group(7)), float(m.group(8)))
    m = g(r"KE gain (\S+) minus piston work (\S+) = (\S+); ledger E - E0 - W (\S+), scale (\S+), ratio (\S+);")
    v["kg"], v["wp"], v["dkw"], v["re"], v["scale"], v["ratio"] = (float(x) for x in m.groups())
    m = g(r"end of the compression t=(\S+) sigma-time \(melt (\S+) \+ compression (\S+)\); events (\d+)"); v["t_end"], v["t_melt"], v["t_comp"] = (float(x) for x in m.groups()[:3]); v["events"] = int(m.group(4))
    v["clean"] = g(r" clean=(\d)").group(1); v["hash"] = g(r" hash=([0-9a-f]{16})").group(1)
    m = g(r"N_L=(\d+) N_R=(\d+) \(moved (\d+)\)"); v["NL"], v["NR"], v["moved"] = (int(x) for x in m.groups())
    m = g(r"walls (\S+), divider (\S+)$"); v["gw"], v["gd"] = (float(x) for x in m.groups())
    m = g(r"contact audit (\d+) events, max gap \[px\] pair (\S+) wall (\S+) divider (\S+) piston (\S+);"); v["aud"] = m.groups()
    m = g(r"eta_0=(\S+) .*?T_melt=(\S+) sigma-time u=(\S+) sigma"); v["eta0"], v["Tm"], v["u"] = m.groups()
    return v


def record(d):
    lg = log(d)
    h = re.findall(r"^(\[EDMD3-HEALTH\] .*?) run_s=\S+ engine_s=\S+ driver_share=\S+$", lg, re.M)
    hs = re.findall(r" hash=([0-9a-f]{16}) ", h[0]) if h else []
    cl = re.findall(r": clean=(\d) ", h[0]) if h else []
    return (h[0] if h else None), (hs[0] if hs else None), (cl[0] if cl else None)


def files(d):
    return {k: one(d, p) for k, p in (("trace", "wall_x_positions_*_run0.csv"), ("psi6(t)", "psi6_t_*_run0.csv"),
                                       ("psi6 summary", "speed_of_sound_psi6.csv"), ("run_params", "run_params.json"))}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--record", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# M4b, the compressed-fluid start: evidence tables (261012 sec. 4.7.30), printed by experiments_gen3_m4b_261009/m4b_tables.py\n")
    rec = open(a.record).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The frozen binaries\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    okb = True
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); okb &= want.get(n) == h; print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '**NO**'} |")
    print("\n## 2. Rule 4: gen2 byte identity with 7b08827 after the 00ALLINONE.c change (ctrl_min, ctrl_leg)\n")
    ok1 = True
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "r4", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "r4", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "r4", f"report_{tag}.txt"), "w").write(r)
        rows = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        ok1 &= len(rows) == 2 and all(l.split("|")[-3].strip() == "IDENTICAL" and l.split("|")[-2].strip() == "IDENTICAL" for l in rows)
        print(f"{'default build' if tag == 'default' else '--engine=gen2'}: " + "; ".join(
            f"{l.split('|')[1].strip()}: audit vs plain {l.split('|')[-3].strip()}, plain vs ref {l.split('|')[-2].strip()}" for l in rows))
    print(f"\nrule 4: {'IDENTICAL (both cases, both builds)' if ok1 else '**NOT IDENTICAL**'}")
    print("\n## 3. The M1 and M2 harness outputs (the committed ones)\n\n| output | identical |\n|---|---|")
    ok2 = True
    for f, ref in REF:
        e = same(os.path.join(O, "h", f), os.path.join(HS, ref)); ok2 &= e
        print(f"| {f} | {'IDENTICAL' if e else '**DIFFERENT or missing**'} |")
    print("\n## 4. The four cells: the compression and the carry (the [EDMD3-M4B] line of each 'plain' run)\n")
    print("| cell | eta_0, T_melt, u | L0 requested [sigma] | achieved - requested left / right [sigma] | N_L / N_R (moved) | KE per compartment: start; "
          "before the rescale; after/N_s kT - 1 | KE gain - piston work (ledger scale) | ledger ratio | end of the compression [sigma-time] (melt + "
          "compression) | events; clean; hash | contact audit: events, max gap pair/wall/div/piston [px] | smallest gap after the carry walls / "
          "divider [px] | psi6 at the end of the compression | run record clean | criteria |\n" + "|---" * 15 + "|")
    ok3 = True
    for c in CELLS:
        d = os.path.join(O, "c", c, "plain"); v = m4b(d); rc, hs, cl = record(d)
        if not v: ok3 = False; print(f"| {c} | **no [EDMD3-M4B] line** (exit {exitcode(d)}) |" + " |" * 13); continue
        crit = (abs(v["Ll"] - v["L0"]) <= 1e-9 and abs(v["Lr"] - v["L0"]) <= 1e-9 and v["moved"] == 0 and v["NL"] == v["NR"]
                and max(abs(x) for x in v["rel"]) <= 1e-12 and abs(v["dkw"]) <= v["scale"] and v["ratio"] <= 1.0 and v["clean"] == "1"
                and cl == "1" and v["psi6"] is not None and v["gw"] > -1e-6 and v["gd"] > -1e-6 and exitcode(d) == "0")
        ok3 &= crit
        ps = re.search(r"global (\S+) local (\S+) neighbours (\S+)", v["psi6"] or "")
        print(f"| {c} | {v['eta0']}, {v['Tm']}, {v['u']} | {v['L0']:.15g} | {v['Ll'] - v['L0']:+.2e} / {v['Lr'] - v['L0']:+.2e} | {v['NL']} / {v['NR']} ({v['moved']}) | "
              f"{v['ke0'][0]:.6f}, {v['ke0'][1]:.6f}; {v['ke1'][0]:.3f}, {v['ke1'][1]:.3f}; {v['rel'][0]:+.1e}, {v['rel'][1]:+.1e} | "
              f"{v['dkw']:+.2e} ({v['scale']:.2e}) | {v['ratio']:.2e} | {v['t_end']:.3f} ({v['t_melt']:g} + {v['t_comp']:.3f}) | {v['events']}; {v['clean']}; {v['hash']} | "
              f"{v['aud'][0]}, {v['aud'][1]}/{v['aud'][2]}/{v['aud'][3]}/{v['aud'][4]} | {v['gw']:.3g} / {v['gd']:.3g} | "
              f"{'global ' + ps.group(1) + ', local ' + ps.group(2) + ', neighbours ' + ps.group(3) if ps else '-'} | {cl} | {'yes' if crit else '**NO**'} |")
    print("\n## 5. The same seed bit-identical ('repeat' against 'plain') and observation does not steer ('obs' against 'plain')\n")
    print("| cell | run | trace | psi6(t) | psi6 summary | run_params | [EDMD3-M4B] lines | run record (no timing) | event hash | criteria |\n" + "|---" * 10 + "|")
    ok4 = ok5 = True
    for c in CELLS:
        U = os.path.join(O, "c", c, "plain"); fu = files(U); ru, hu, _ = record(U); mu = m4b_lines(U)
        for tag in ("repeat", "obs"):
            d = os.path.join(O, "c", c, tag); f = files(d); r, h, _ = record(d); ml = m4b_lines(d)
            e = {k: same(f[k], fu[k]) for k in fu}
            if tag == "repeat":
                crit = all(e.values()) and r == ru and ml == mu and h == hu and exitcode(d) == "0"; ok4 &= crit
            else:
                crit = e["trace"] and e["psi6 summary"] and ml == mu and h == hu and exitcode(d) == "0"; ok5 &= crit
            print(f"| {c} | {tag} | {'IDENTICAL' if e['trace'] else 'different'} | {'IDENTICAL' if e['psi6(t)'] else 'different (psi6 cadence)' if tag == 'obs' else 'different'} | "
                  f"{'IDENTICAL' if e['psi6 summary'] else 'different'} | {'IDENTICAL' if e['run_params'] else 'different'} | {'IDENTICAL' if ml == mu else 'different'} | "
                  f"{'IDENTICAL' if r == ru else 'different'} | {h} {'=' if h == hu else '!='} {hu} | {'yes' if crit else '**NO**'} |")
    print("\n## 6. Checkpoint and restart through the compression (N = 400, eta 0.716; U = the cell's 'plain' run)\n")
    kb = os.path.join(O, "k"); U = os.path.join(O, "c", "N400_eta0.716", "plain"); fu = files(U); ru, hu, _ = record(U); mu = m4b_lines(U)
    ok6 = same(os.path.join(kb, "ck_C.bin"), os.path.join(kb, "ck_R1.bin"))
    print(f"checkpoint files ck_C.bin / ck_R1.bin: {os.path.getsize(os.path.join(kb, 'ck_C.bin')) if os.path.exists(os.path.join(kb, 'ck_C.bin')) else '-'} bytes, "
          f"{'IDENTICAL' if ok6 else '**DIFFERENT or missing**'}\n")
    print("| run | exit | trace | psi6(t) | psi6 summary | run record (no timing) | event hash | [EDMD3-M4B] lines = U's | restart line |\n" + "|---" * 9 + "|")
    for r_ in ("C", "R1", "R2"):
        d = os.path.join(kb, r_); f = files(d); r, h, _ = record(d); ml = m4b_lines(d)
        rs = re.findall(r"^\[EDMD3-RESTART\] .*$", log(d), re.M)
        e = {k: same(f[k], fu[k]) for k in fu}
        if r_ != "R1":
            ok6 &= e["trace"] and e["psi6(t)"] and e["psi6 summary"] and r == ru and h == hu and ml == mu and exitcode(d) == "0"
        print(f"| {r_} | {exitcode(d)} | {'IDENTICAL' if e['trace'] else '-' if r_ == 'R1' else 'different'} | {'IDENTICAL' if e['psi6(t)'] else '-' if r_ == 'R1' else 'different'} | "
              f"{'IDENTICAL' if e['psi6 summary'] else '-' if r_ == 'R1' else 'different'} | {'IDENTICAL' if r == ru else '-' if r_ == 'R1' else 'different'} | {h} | "
              f"{'yes' if ml == mu else 'no'} | {rs[0][:120] if rs else '-'} |")
    print("\n## 7. Refusals (each must stop with exit 2)\n\n| case | exit | message |\n|---|---|---|")
    ok7 = True
    for n in sorted(glob.glob(os.path.join(O, "x", "case*.name"))):
        d = n[:-5]; ec = exitcode(d); lg = log(d)
        m = re.findall(r"^(STOP.*)$", lg, re.M); ok7 &= ec == "2" and bool(m)
        print(f"| {open(n).read().strip()} | {ec} | {m[-1][:160] if m else '**no STOP line**'} |")
    allok = okb and ok1 and ok2 and ok3 and ok4 and ok5 and ok6 and ok7
    print(f"\nACCEPTANCE (M4b, decision 13 item 4): {'PASS' if allok else '**NOT MET**'} -- binaries = build record ({'yes' if okb else 'NO'}); rule 4 "
          f"({'yes' if ok1 else 'NO'}); M1/M2 harness outputs ({'yes' if ok2 else 'NO'}); the four cells: achieved L0, N_L = N_R, KE after the rescale, piston "
          f"work = KE gain, clean, psi6 printed ({'yes' if ok3 else 'NO'}); same seed bit-identical ({'yes' if ok4 else 'NO'}); observation does not steer "
          f"({'yes' if ok5 else 'NO'}); checkpoint/restart through the compression ({'yes' if ok6 else 'NO'}); refusals ({'yes' if ok7 else 'NO'})")


if __name__ == "__main__":
    main()
