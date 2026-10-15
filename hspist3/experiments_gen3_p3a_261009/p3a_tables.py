#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.26; plan-author decision 12, part 3a): the tables of part 3a, printed from the evidence folder
of p3a_evidence.py. Criteria stated before the run (in this script's commit):
  - observation does not steer: the event hash with --gen3-virial-blocks and/or --gen3-snapshots equals the plain run's, and so do
    the trace, the psi6(t) file and the psi6 summary;
  - the symmetric held cell: per run, Z_L and Z_R (means over the hold blocks, SE = SD / sqrt(blocks)) agree within their SEs:
    abs(Z_L - Z_R) / sqrt(SE_L^2 + SE_R^2) < 2; and in every block the KE-weighted mean of the compartments' Z equals the global Z
    (computed by the engine from its own global sum) to 1e-9 relative. Why 1e-9: the block values are differences of cumulative
    double sums of up to ~1e7 pair terms, whose rounding is ~ u sqrt(n) |W| ~ 1e-11 of a block's W late in a 5000-sigma-time run
    (a first smoke run showed 6.6e-13 absolute in its second block), while one pair term lost or counted twice moves a block's Z by
    ~4e-6 relative, so 1e-9 separates the two by three orders of magnitude either way;
  - rule 4 IDENTICAL, the M1/M2 harness outputs byte-identical, the engine checkpoint test PASS (stage H's tables of the rerun);
  - restarts with both observations on byte-identical to the uninterrupted run (trace, psi6(t), psi6 summary, the virial blocks, the
    snapshots, the event hash); two checkpoints of the same state byte-identical (part 3b), here and in stage H's rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p3a_261009/p3a_tables.py --bin-dir <frozen binaries> --record <build record> --out <evidence dir>
"""
import argparse, csv, glob, hashlib, math, os, re, struct, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, "rb").read() if p and os.path.exists(p) else None
def same(p, q): a, b = rd(p), rd(q); return a is not None and b is not None and a == b
def one(d, pat):
    g = glob.glob(os.path.join(d, pat)); return g[0] if len(g) == 1 else None


def record(d):
    lg = open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    h = re.findall(r"^(\[EDMD3-HEALTH\] .*?) run_s=\S+ engine_s=\S+ driver_share=\S+$", lg, re.M)
    rec = h[0] if h else None
    kv = dict(re.findall(r" (\w+)=(\S+)", rec)) if rec else {}
    return rec, kv.get("hash"), kv


def exitcode(d):
    p = os.path.join(d, "exit_code.txt"); return open(p).read().strip() if os.path.exists(p) else "?"


def blocks(d):
    p = one(d, "virial_blocks_*_run0.csv")
    return list(csv.DictReader(open(p))) if p else []


def mse(v):
    v = [x for x in v if math.isfinite(x)]
    if len(v) < 2: return float("nan"), float("nan"), len(v)
    m = sum(v) / len(v); sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)); return m, sd / math.sqrt(len(v)), len(v)


def read_snaps(p):
    b = open(p, "rb").read()
    assert b[:8] == b"G3SNAP1\n", "magic"
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, T = struct.unpack_from("<4d", b, 16); o = 48; out = []
    while o < len(b):
        t, = struct.unpack_from("<d", b, o); ph, k = struct.unpack_from("<ii", b, o + 8); o += 16
        dx = struct.unpack_from(f"<{k}d", b, o); o += 8 * k
        xy = struct.unpack_from(f"<{2 * N}d", b, o); o += 16 * N
        out.append((t, ph, dx, xy))
    return (N, nd, W, H, R, T), out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--record", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# Decision 12, part 3a: the virial per compartment and the position snapshots (261012 sec. 4.7.26), printed by "
          "experiments_gen3_p3a_261009/p3a_tables.py\n")
    rec = open(a.record).read(); want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    print("## 1. The frozen binaries\n\n| binary | --version | SHA-256 | = build record |\n|---|---|---|---|")
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2", "gen3_checkpoint_test", "gen3_checkpoint_test_ld"):
        v = subprocess.run([os.path.join(B, n), "--version"], capture_output=True, text=True).stdout.splitlines()[0] if n.startswith("00") else "-"
        h = sha(os.path.join(B, n)); print(f"| {n} | {v} | {h} | {'yes' if want.get(n) == h else '**NO**'} |")
    # ---- 2. stage H's tables of the rerun, and its checkpoint pairs
    print("\n## 2. Stage H's evidence runner with these binaries: its tables (stageH_tables.py, unchanged; verbatim)\n\n```")
    r = subprocess.run([sys.executable, os.path.join(HS, "experiments_gen3_h_261009", "stageH_tables.py"), "--bin-dir", B, "--record", a.record,
                        "--out", os.path.join(O, "h")], capture_output=True, text=True)
    print(r.stdout.rstrip() + "\n```")
    accH = "ACCEPTANCE (stage H): PASS" in r.stdout
    print("\nTwo checkpoints of the same state (part 3b), stage H's rerun:\n\n| pair | bytes | identical |\n|---|---|---|")
    okp = nP = 0
    for cell, cks in (("A", ("hold1", "hold300000", "record0", "record25000")), ("B", ("hold60000", "record6000")), ("C", ("hold15000", "record3000"))):
        for t in cks:
            p1, p2 = os.path.join(O, "h", "h3", cell, f"ck_C_{t}.bin"), os.path.join(O, "h", "h3", cell, f"ck_R1_{t}.bin")
            e = same(p1, p2); okp += e; nP += 1
            print(f"| {cell} {t}: C / R1 | {len(rd(p1) or b'')} | {'IDENTICAL' if e else '**DIFFERENT**'} |")
    p1, p2 = os.path.join(O, "h", "h3", "A", "ck_C_record25000.bin"), os.path.join(O, "h", "h3", "A", "ck_R2b_chain.bin")
    e = same(p1, p2); okp += e; nP += 1
    print(f"| A record25000: C / R2b (restarted at hold:300000) | {len(rd(p1) or b'')} | {'IDENTICAL' if e else '**DIFFERENT**'} |")
    # ---- 3. e3: observation does not steer
    print("\n## 3. Observation does not steer (P1's N = 400, eta 0.704 cell, M = 50, held 2000 sigma-time, 100-sigma-time record)\n")
    print("| run | exit | event hash | = plain | trace = plain | psi6(t) = plain | psi6 summary = plain | run-record fields that differ from plain |\n|---|---|---|---|---|---|---|---|")
    E3 = os.path.join(O, "e3"); _, hp, kvp = record(os.path.join(E3, "plain")); ok3 = True
    for tag in ("plain", "vb", "snap", "both"):
        d = os.path.join(E3, tag); _, hh, kv = record(d)
        eq = [same(one(d, pat), one(os.path.join(E3, "plain"), pat)) for pat in ("wall_x_positions_*_run0.csv", "psi6_t_*_run0.csv", "speed_of_sound_psi6.csv")]
        diff = sorted(k for k in set(kv) | set(kvp) if kv.get(k) != kvp.get(k))
        if tag != "plain": ok3 &= hh == hp and all(eq) and exitcode(d) == "0"
        print(f"| {tag} | {exitcode(d)} | {hh} | {'(reference)' if tag == 'plain' else ('yes' if hh == hp else '**NO**')} | " +
              " | ".join('(reference)' if tag == 'plain' else ('IDENTICAL' if x else '**DIFFERENT**') for x in eq) + f" | {', '.join(diff) or '-'} |")
    # ---- 4. e4: the symmetric held cell
    print("\n## 4. The symmetric held cell: Z per compartment over the hold blocks (50 sigma-time each)\n")
    print("| cell | seed | blocks | N_L / N_R | KE_L / KE_R | Z_L (SE) | Z_R (SE) | Z_L - Z_R (SE) | z | Z_global (SE) | max abs(KE-weighted mean - global) / Z | criteria |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
    ok4 = True
    for d in sorted(glob.glob(os.path.join(O, "e4", "*"))):
        bl = [x for x in blocks(d) if x["phase"] == "hold"]
        if not bl: print(f"| {os.path.basename(d)} | no blocks |"); ok4 = False; continue
        zl, sl, n = mse([float(x["Z_0"]) for x in bl]); zr, sr, _ = mse([float(x["Z_1"]) for x in bl]); zg, sg, _ = mse([float(x["Z_global"]) for x in bl])
        dz = zl - zr; sdz = math.sqrt(sl ** 2 + sr ** 2); z = dz / sdz
        wmax = max(abs(float(x["weighted_minus_global"])) / abs(float(x["Z_global"])) for x in bl)
        ok = abs(z) < 2 and wmax <= 1e-9 and exitcode(d) == "0"; ok4 &= ok
        cell, seed = os.path.basename(d).split("_seed")
        print(f"| {cell} | {seed} | {n} | {bl[0]['N_0']} / {bl[0]['N_1']} | {float(bl[-1]['KE_0']):.6g} / {float(bl[-1]['KE_1']):.6g} | {zl:.5f} ({sl:.5f}) | "
              f"{zr:.5f} ({sr:.5f}) | {dz:+.5f} ({sdz:.5f}) | {z:+.2f} | {zg:.5f} ({sg:.5f}) | {wmax:.1e} | {'yes' if ok else '**NO**'} |")
    # ---- 5. e5: the snapshots' content
    print("\n## 5. The snapshots against the initial state (--gen3-dump-initial), and their bookkeeping\n")
    d5 = os.path.join(O, "e5"); sp = one(d5, "snapshots_*_run0.bin"); ok5 = False
    if sp and os.path.exists(os.path.join(d5, "initial_state.txt")):
        (N, nd, W, H, R, Ti), snaps = read_snaps(sp)
        lines = open(os.path.join(d5, "initial_state.txt")).read().split("\n")
        xy0 = [tuple(float.fromhex(v) / 24.0 for v in l.split()[:2]) for l in lines[6:6 + N]]
        t0, ph0, dx0, xy = snaps[0]
        dmax = max(max(abs(xy[2 * i] - xy0[i][0]), abs(xy[2 * i + 1] - xy0[i][1])) for i in range(N))
        left = [sum(1 for i in range(N) if s[3][2 * i] < s[2][0]) for s in snaps]
        times = [round(s[0], 6) for s in snaps]
        ok5 = dmax == 0.0 and len(set(left)) == 1 and abs(times[0]) < 1e-12 and exitcode(d5) == "0"
        print(f"header: N {N}, dividers {nd}, box {W:.6f} x {H:.6f} sigma, radius {R}, interval {Ti} sigma-time")
        print(f"snapshots: {len(snaps)} at t = {times[0]} ... {times[-1]} sigma-time (phases {sorted(set(s[1] for s in snaps))})")
        print(f"first snapshot (t = 0) against the dumped initial state / 24: max abs difference {dmax:.3g} sigma over {N} disks (0 = bit-identical)")
        print(f"disks left of the divider in every snapshot: {sorted(set(left))}")
    else:
        print("**the snapshot file or the initial dump is missing**")
    # ---- 6. e6: a restart with both observations on
    print("\n## 6. A restart with both observations on (e3's command with both flags): every finished run against U\n")
    print("| run | exit | trace | psi6(t) | psi6 summary | virial blocks | snapshots | run record (without timing) | event hash |\n|---|---|---|---|---|---|---|---|---|")
    E6 = os.path.join(O, "e6"); U = os.path.join(E6, "U"); ru, hu, _ = record(U); ok6 = True
    pats = ("wall_x_positions_*_run0.csv", "psi6_t_*_run0.csv", "speed_of_sound_psi6.csv", "virial_blocks_*_run0.csv", "snapshots_*_run0.bin")
    for name in sorted(os.listdir(E6)):
        d = os.path.join(E6, name)
        if not os.path.isdir(d) or name == "U" or name.startswith("R1"): continue
        rr, hh, _ = record(d)
        eq = [same(one(d, p), one(U, p)) for p in pats]
        ok6 &= all(eq) and rr == ru and hh == hu and exitcode(d) == "0"
        print(f"| {name} | {exitcode(d)} | " + " | ".join('IDENTICAL' if x else '**DIFFERENT**' for x in eq) +
              f" | {'IDENTICAL' if rr == ru else '**DIFFERENT**'} | {hh}{'' if hh == hu else ' (U: ' + str(hu) + ')'} |")
    print("\n| checkpoint pair (same state) | bytes | identical |\n|---|---|---|")
    for t in ("hold60000", "record3000"):
        p1, p2 = os.path.join(E6, f"ck_C_{t}.bin"), os.path.join(E6, f"ck_R1_{t}.bin")
        e = same(p1, p2); okp += e; nP += 1; ok6 &= e
        print(f"| {t}: C / R1 | {len(rd(p1) or b'')} | {'IDENTICAL' if e else '**DIFFERENT**'} |")
    print(f"\nACCEPTANCE (part 3a): {'PASS' if ok3 and ok4 and ok5 and ok6 and accH else '**NOT MET**'} -- observation does not steer "
          f"({'yes' if ok3 else 'NO'}); held cells: compartments agree and the weighted mean is the global Z ({'yes' if ok4 else 'NO'}); "
          f"snapshots exact ({'yes' if ok5 else 'NO'}); restarts with the observations byte-identical ({'yes' if ok6 else 'NO'}); "
          f"stage H's tables PASS with these binaries, rule 4 included ({'yes' if accH else 'NO'})")
    print(f"ACCEPTANCE (part 3b, with these binaries): {'PASS' if okp == nP else '**NOT MET**'} -- checkpoint pairs of the same state "
          f"byte-identical: {okp} of {nP}")


if __name__ == "__main__":
    main()
