#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 17, "after P5"; 261012 sec. 4.7.40): the acceptance table of the PROVISIONAL psi6 per
compartment, from psi6c_evidence.py's folder. The items and criteria, written before the runs:
  1. the binaries: the new build's --version names the code commit and its SHA-256 equals its build record; P4's frozen binary
  2. rule 4: resched_gate_261005/audit_runs_261007.py report against 7b08827's outputs (the scratch reference stage H used):
     IDENTICAL for ctrl_min and ctrl_leg, the default build and --engine=gen2
  3. the M1/M2 harness outputs byte-identical to the committed ones (experiments_gen3_m1_261008, experiments_gen3_m2_261009)
  4. S2, the columns on and off (P4a's start-A seed-1 command held 2000 sigma-time, snapshots every 500):
     a. old (06e7248), off and on: exit 0, one run record with clean=1, the same event hash
     b. old vs off: every data file byte-identical (psi6(t), virial blocks, snapshots, trace, speed_of_sound_psi6.csv)
     c. on vs off: the psi6(t) file's first five columns byte-identical in the header and in every row; every other data file
        byte-identical; every row of the on file has 5 + 4 x (compartments) fields
  5. S3, P4a's start-A seed-1 run replayed with the columns on, against P4's original run (06e7248, the columns off):
     a. P4's original files unchanged since P4 (SHA-256 = experiments_gen3_p4_261010/P4_SHA256_FILES.txt)
     b. exit 0, clean=1, and the event hash = P4's
     c. the psi6(t) file's first five columns = P4's file, byte for byte; the virial blocks and the snapshots byte-identical to P4's
     d. (printed) the new columns every 1e4 sigma-time
  6. S4, --gen3-psi6-snapshots:
     a. on P4a's four snapshot files, at 0, 1e4, ..., 1e5 sigma-time: |Psi_c| within 6e-5 and the orientation within 0.06 deg (mod 60)
        of sec. 4.7.37's Python values (d16_p4a_snapshot_psi6_output.txt, printed to 4 and 1 decimals)
     b. on S2 on's and S3's own snapshot files: at every snapshot time the reader's per-compartment values = the run's columns at
        that time (N exact; |Psi| and local within 2e-6; orientation within 2e-3 deg mod 60: the run computes in px, the reader in
        sigma), and the reader's global = the run's psi6_global within 2e-6
     c. SEPARATION (the plan author's acceptance): in P4a's A seed 1 at 5e4, 6e4, ..., 1e5 sigma-time (the times sec. 4.7.37
        printed) compartment 0 at |orientation| > 25 deg and compartment 1 at |orientation| < 5 deg, both with |Psi| > 0.5; at 0, 1e4,
        ..., 4e4 both within 5 deg of 0
  7. S5, the refusals: exit 2 and the STOP line naming --gen3-psi6-compartments
  ACCEPTANCE: PASS if every item holds; else FAIL with the items that do not.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_psi6c_261010/psi6c_tables.py --bin-dir <clean build dir> --out <evidence dir>
"""
import argparse, csv, glob, hashlib, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import psi6c_evidence as EV
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REF = (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
       ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
       ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt"))
DATA = ("psi6_t_", "virial_blocks_", "snapshots_", "wall_x_positions_", "speed_of_sound_psi6.csv")
PY = os.path.join(HS, "experiments_gen3_d16_261010", "d16_p4a_snapshot_psi6_output.txt")
TIMES = [1e4 * k for k in range(11)]


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, "rb").read() if p and os.path.exists(p) else None
def one(d, pre):
    f = sorted(x for x in (os.listdir(d) if os.path.isdir(d) else []) if x.startswith(pre))
    return os.path.join(d, f[0]) if len(f) == 1 else None
def record(d):
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    recs = re.findall(r"^\[EDMD3-HEALTH\] (.*)$", lg, re.M)
    kv = dict(re.findall(r" (\w+)=(\S+)", " " + recs[0].split(": ", 1)[-1])) if len(recs) == 1 else {}
    ec = open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
    return dict(exit=ec, records=len(recs), clean=kv.get("clean"), hash=kv.get("hash"), build=kv.get("build"))
def odiff(a, b): return abs(((a - b + 30.0) % 60.0) - 30.0)
def rows(p):
    return list(csv.DictReader(open(p))) if p and os.path.exists(p) else []
def rows_nc(p):   # the reader's output: '#' lines (its header line and its last line) skipped
    return list(csv.DictReader(l for l in open(p) if not l.startswith("#"))) if p and os.path.exists(p) else []


def first5(p):
    return [",".join(l.split(",")[:5]) for l in open(p, "rb").read().decode().split("\n") if l]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); new = os.path.join(B, "00ALLINONE"); fails = []
    def item(name, ok, txt):
        print(f"| {name} | {txt} | {'yes' if ok else '**NO**'} |")
        if not ok: fails.append(name)
    print("# psi6 per compartment (PROVISIONAL): the acceptance table, printed by experiments_gen3_psi6c_261010/psi6c_tables.py "
          "(261012 sec. 4.7.40; plan-author decision 17)\n")
    print("| item | result | holds |\n|---|---|---|")
    # 1
    rec = open(os.path.join(B, "build_record.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    ver = subprocess.run([new, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    item("1 the new build", want.get("00ALLINONE") == sha(new) and "dirty" not in ver, f"{ver}; SHA-256 {sha(new)[:16]}... = build record")
    item("1 P4's frozen binary", sha(EV.OLD_BIN) == EV.OLD_SHA, f"{EV.OLD_BIN.split('hspist3/')[-1]} SHA-256 {sha(EV.OLD_BIN)[:16]}...")
    # 2 rule 4
    ok0, txt0 = True, []
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "s0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(EV.GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "s0", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "s0", f"report_{tag}.txt"), "w").write(r)
        rr = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        ok = len(rr) == 2 and all(l.split("|")[-3].strip() == "IDENTICAL" and l.split("|")[-2].strip() == "IDENTICAL" for l in rr); ok0 &= ok
        txt0.append(f"{'default' if tag == 'default' else '--engine=gen2'}: " + ", ".join(f"{l.split('|')[1].strip()} {l.split('|')[-2].strip()}" for l in rr))
    item("2 rule 4 (gen2 vs 7b08827)", ok0, "; ".join(txt0))
    # 3 M1/M2
    ok3 = all(rd(os.path.join(O, "s1", f)) is not None and rd(os.path.join(O, "s1", f)) == rd(os.path.join(HS, ref)) for f, ref in REF)
    item("3 M1/M2 harness outputs", ok3, ", ".join(f"{f}: {'IDENTICAL' if rd(os.path.join(O, 's1', f)) == rd(os.path.join(HS, ref)) else 'DIFFERENT'}" for f, ref in REF))
    # 4 S2
    R2 = {t: record(os.path.join(O, "s2", t)) for t in ("old", "off", "on")}
    ok = all(R2[t]["exit"] == "0" and R2[t]["records"] == 1 and R2[t]["clean"] == "1" for t in R2) and len({R2[t]["hash"] for t in R2}) == 1
    item("4a S2 event hash old = off = on", ok, ", ".join(f"{t}: exit {R2[t]['exit']}, clean {R2[t]['clean']}, hash {R2[t]['hash']}, build {R2[t]['build']}" for t in R2))
    d = {t: os.path.join(O, "s2", t) for t in R2}
    ok = all(one(d["old"], p) and rd(one(d["old"], p)) == rd(one(d["off"], p)) for p in DATA)
    item("4b S2 old vs off: data files", ok, ", ".join(f"{p.rstrip('_')}: {'IDENTICAL' if one(d['old'], p) and rd(one(d['old'], p)) == rd(one(d['off'], p)) else 'DIFFERENT/missing'}" for p in DATA))
    pon, poff = one(d["on"], "psi6_t_"), one(d["off"], "psi6_t_")
    f_on, f_off = first5(pon), first5(poff)
    nf = {len(l.split(",")) for l in open(pon).read().split("\n")[1:] if l}
    others = all(rd(one(d["on"], p)) == rd(one(d["off"], p)) for p in DATA if p != "psi6_t_")
    item("4c S2 on vs off", f_on == f_off and others and nf == {13},
         f"psi6(t) first five columns {'IDENTICAL' if f_on == f_off else 'DIFFERENT'} ({len(f_on)} lines incl. header); fields per row {sorted(nf)}; "
         f"other data files {'IDENTICAL' if others else 'DIFFERENT'}; on header: {open(pon).readline().strip()}")
    # 5 S3
    plist = dict((n.split("/", 1)[1], h) for h, n in re.findall(r"^([0-9a-f]{64})  (data_06e7248/\S+)$",
                                                                open(os.path.join(HS, "experiments_gen3_p4_261010", "P4_SHA256_FILES.txt")).read(), re.M))
    origs = [one(EV.ORIG, p) for p in ("psi6_t_", "virial_blocks_", "snapshots_")]
    ok = all(o and plist.get(os.path.relpath(o, os.path.dirname(os.path.dirname(os.path.dirname(EV.ORIG))))) == sha(o) for o in origs)
    item("5a P4's original files unchanged", ok, ", ".join(f"{os.path.basename(o)[:13]}: {'= P4_SHA256_FILES' if o and plist.get(os.path.relpath(o, os.path.dirname(os.path.dirname(os.path.dirname(EV.ORIG))))) == sha(o) else 'NOT in the list or changed'}" for o in origs))
    r3, rp = record(os.path.join(O, "s3")), record(EV.ORIG)
    item("5b S3 event hash = P4's", r3["exit"] == "0" and r3["clean"] == "1" and r3["hash"] == rp["hash"],
         f"S3: exit {r3['exit']}, clean {r3['clean']}, hash {r3['hash']}, build {r3['build']}; P4: hash {rp['hash']}, build {rp['build']}")
    s3 = os.path.join(O, "s3")
    a5 = first5(one(s3, "psi6_t_")) == first5(one(EV.ORIG, "psi6_t_"))
    b5 = all(rd(one(s3, p)) == rd(one(EV.ORIG, p)) for p in ("virial_blocks_", "snapshots_"))
    item("5c S3 files = P4's", a5 and b5, f"psi6(t) first five columns {'IDENTICAL' if a5 else 'DIFFERENT'}; virial blocks and snapshots {'IDENTICAL' if b5 else 'DIFFERENT'}")
    # 6a, 6c: the reader on P4a's snapshots against the Python values
    py = {}
    for l in open(PY):
        m = re.match(r"\| (A|B) \| (\d) \| [^|]+\| (\d+) \| \d+ / \d+ \| ([0-9.]+) \| ([0-9.]+) \| ([-+0-9.]+) \| ([-+0-9.]+) \|", l)
        if m: py[(m.group(1), int(m.group(2)), float(m.group(3)))] = tuple(float(m.group(k)) for k in range(4, 8))
    ok6a, worst_a, worst_o, n6 = True, 0.0, 0.0, 0; S4 = {}
    for st, r in EV.SNAPS:
        S4[(st, r)] = {round(float(x["t_sigma"])): x for x in rows_nc(os.path.join(O, "s4", f"p4a_{st}{r}", "psi6_compartments.csv"))}
        for t in TIMES:
            x = S4[(st, r)].get(round(t)); p = py.get((st, r, t))
            if x is None or p is None: ok6a = False; continue
            da = max(abs(float(x["psi6_0"]) - p[0]), abs(float(x["psi6_1"]) - p[1]))
            do = max(odiff(float(x["psi6_orient_deg_0"]), p[2]), odiff(float(x["psi6_orient_deg_1"]), p[3]))
            worst_a, worst_o, n6 = max(worst_a, da), max(worst_o, do), n6 + 1
            ok6a &= da <= 6e-5 and do <= 0.06
    item("6a reader on P4a's snapshots vs the Python values", ok6a and n6 == 44, f"{n6} of 44 rows; largest difference |Psi| {worst_a:.1e}, orientation {worst_o:.3f} deg")
    A1 = S4[("A", 1)]
    sep = all(abs(float(A1[round(t)]["psi6_orient_deg_0"])) > 25 and abs(float(A1[round(t)]["psi6_orient_deg_1"])) < 5
              and float(A1[round(t)]["psi6_0"]) > 0.5 and float(A1[round(t)]["psi6_1"]) > 0.5 for t in TIMES[5:] if round(t) in A1)
    pre = all(abs(float(A1[round(t)]["psi6_orient_deg_0"])) < 5 and abs(float(A1[round(t)]["psi6_orient_deg_1"])) < 5 for t in TIMES[:5] if round(t) in A1)
    item("6c SEPARATION in P4a's A seed 1", sep and pre and all(round(t) in A1 for t in TIMES),
         "5e4-1e5: compartment 0 at " + ", ".join(f"{float(A1[round(t)]['psi6_orient_deg_0']):+.1f}" for t in TIMES[5:] if round(t) in A1)
         + " deg, compartment 1 at " + ", ".join(f"{float(A1[round(t)]['psi6_orient_deg_1']):+.1f}" for t in TIMES[5:] if round(t) in A1)
         + "; 0-4e4: both within " + f"{max(max(abs(float(A1[round(t)]['psi6_orient_deg_0'])), abs(float(A1[round(t)]['psi6_orient_deg_1']))) for t in TIMES[:5] if round(t) in A1):.1f} deg of 0")
    # 6b: the reader on the runs' own snapshots against the runs' columns
    for tag, run_d, rd_d in (("S2 on", os.path.join(O, "s2", "on"), os.path.join(O, "s4", "s2_on")), ("S3", s3, os.path.join(O, "s4", "s3"))):
        cols = {round(float(x["t_sigma"]), 6): x for x in rows(one(run_d, "psi6_t_"))}
        rdr = rows_nc(os.path.join(rd_d, "psi6_compartments.csv"))
        ok, n, wa, wo, wg = bool(rdr), 0, 0.0, 0.0, 0.0
        for x in rdr:
            y = cols.get(round(float(x["t_sigma"]), 6))
            if y is None: ok = False; continue
            n += 1
            for c in (0, 1):
                ok &= x[f"N_{c}"] == y[f"N_{c}"]
                wa = max(wa, abs(float(x[f"psi6_{c}"]) - float(y[f"psi6_{c}"])), abs(float(x[f"psi6_local_{c}"]) - float(y[f"psi6_local_{c}"])))
                wo = max(wo, odiff(float(x[f"psi6_orient_deg_{c}"]), float(y[f"psi6_orient_deg_{c}"])))
            wg = max(wg, abs(float(x["psi6_global"]) - float(y["psi6_global"])))
        ok &= wa <= 2e-6 and wo <= 2e-3 and wg <= 2e-6
        item(f"6b reader = the run's columns ({tag})", ok and n == len(rdr), f"{n} of {len(rdr)} snapshot times matched; largest difference |Psi|/local {wa:.1e}, "
             f"orientation {wo:.1e} deg, global {wg:.1e}")
    # 7 refusals
    for tag, msg in (("gen2", "--gen3-psi6-compartments needs --engine=gen3"), ("checkpoint", "not carried by a checkpoint"),
                     ("energy_transfer", "wired to the speed-of-sound experiment only")):
        dd = os.path.join(O, "s5", tag); lg = open(os.path.join(dd, "run.log"), errors="ignore").read()
        ec = open(os.path.join(dd, "exit_code.txt")).read().strip()
        st = [l for l in lg.splitlines() if l.startswith("STOP")]
        item(f"7 refusal: {tag}", ec == "2" and any(msg in l and "--gen3-psi6-compartments" in l for l in st), f"exit {ec}; {st[0] if st else 'no STOP line'}")
    # 5d: the new columns of S3 every 1e4 sigma-time
    print("\n## S3: the new columns of P4a's start-A seed-1 run, every 1e4 sigma-time (the columns on; the trajectory = P4's)\n")
    print("| t [sigma-time] | psi6 global | N_0 | psi6_0 | orientation_0 [deg] | local_0 | N_1 | psi6_1 | orientation_1 [deg] | local_1 |\n" + "|---" * 10 + "|")
    c3 = {round(float(x["t_sigma"])): x for x in rows(one(s3, "psi6_t_")) if x["phase"] == "hold"}
    for t in TIMES:
        x = c3.get(round(t))
        if x: print(f"| {t:.0f} | {x['psi6_global']} | {x['N_0']} | {x['psi6_0']} | {x['psi6_orient_deg_0']} | {x['psi6_local_0']} | {x['N_1']} | {x['psi6_1']} | "
                    f"{x['psi6_orient_deg_1']} | {x['psi6_local_1']} |")
    print("\n## S4: --gen3-psi6-snapshots on P4a's snapshots, every 1e4 sigma-time\n")
    print("| run | t | psi6 global | psi6_0 | orientation_0 | psi6_1 | orientation_1 |\n|---|---|---|---|---|---|---|")
    for st, r in EV.SNAPS:
        for t in TIMES:
            x = S4[(st, r)].get(round(t))
            if x: print(f"| {st} seed {r} | {t:.0f} | {x['psi6_global']} | {x['psi6_0']} | {x['psi6_orient_deg_0']} | {x['psi6_1']} | {x['psi6_orient_deg_1']} |")
    print(f"\nACCEPTANCE: {'PASS' if not fails else '**FAIL**: ' + '; '.join(fails)}")


if __name__ == "__main__":
    main()
