#!/usr/bin/env python3
"""##CHRIS 2026-10-10 (plan-author decision 18, items 1 and 2; 261012 sec. 4.7.42): the acceptance tables of items 1 and 2, from
d18_evidence.py's folder. The items and criteria, written before the runs:
ITEM 1 -- checkpoint/restart carries --gen3-psi6-compartments
  1.1 the binaries: the new build's --version names its commit; 00ALLINONE, gen3_m1, gen3_m2 and both checkpoint tests have the
      SHA-256 of the build record
  1.2 stage H's tables pass again: experiments_gen3_h_261009/stageH_tables.py, unchanged, on stage H's evidence of the new build
      ends with "ACCEPTANCE (stage H): PASS" (that run holds rule 4, H0, and the M1/M2 harness outputs, H1)
  1.3 a restarted run with the columns on is byte-identical to the uninterrupted one: for cells B and C and each checkpoint, C (written
      and going on) and R2 (restarted) against U: the trace, the psi6(t) file with its new columns, speed_of_sound_psi6.csv,
      run_params.json, and for B the virial blocks and the snapshots, byte for byte; the run record without its timing fields and the
      event hash equal; every run exit 0; U's psi6(t) header carries the per-compartment columns
  1.4 the refusals: a columns-on checkpoint restarted without the flag, and a columns-off checkpoint (stage H's) restarted with it:
      exit 2 and the STOP line naming --gen3-psi6-compartments
ITEM 2 -- M4c, the turned orientation (--gen3-lattice-orientation=0|turned)
  2.1 the default and existing commands, the previous build (0d7f1ca) against the new one, byte for byte:
      a. P4a's start-A seed-1 run (new build, default): the event hash bd25077c5f1c9c9f, and the snapshots (the initial state at t = 0
         among them), the psi6(t) file and the virial blocks equal to P4's originals (whose SHA-256 equal P4_SHA256_FILES.txt)
      b. stage C's grid cell N = 400, eta 0.70 ("long"): the initial state (init.txt), the psi6(t) file and the event hash, old = new
         (and against stage C's own evidence f39e485: information)
      c. the same cell with --gen3-lattice=perpendicular, and P4's M4b command: initial state (where dumped), psi6(t), hash, old = new
  2.2 the turned grid (N = 100, 400, 900, 1600 x eta 0.60, 0.70, 0.716, 0.74, 0.78, 0.85), per cell either INFEASIBLE (exit 2, the
      [EDMD3-M4] INFEASIBLE line) or: the achieved eta = the requested to 1e-12 (relative); N_L = N_R = N/2; the smallest initial gaps
      disk-disk, disk-wall and disk-divider > 0; "same" bit-identical to "long" (init.txt), "other" different; no ties at the load
      (no live event sharing its time exactly with another, none due at once); "long" clean=1 over 400 sigma-time, exit 0; the
      [EDMD3-M4C] line names lines along y in both compartments
  2.3 the orientation column at the load (the psi6(t) row at t = 0):
      a. (the code) the turned start in both compartments at |orientation| > 25 deg, and --gen3-lattice-orientation=0 within 5 deg
         of 0, in every feasible cell;
      b. (the plan author's premise, reported, not counted in the acceptance) the default family (the present M4) within 5 deg of 0;
         where it is not, the cell is listed. Changed before the runs: the replica and P4's own N = 1600 snapshots (read with
         --gen3-psi6-snapshots of 0d7f1ca) showed that the present family, which takes the arrangement with the largest gap of
         either orientation, is the turned crystal at N = 1600, eta 0.74 and 0.78 -- P4's start A there reads +-30 deg at t = 0
  2.4 the feasibility ceiling per N (H = 10 sqrt(N/100)) for the default family, 0 and turned, by a Python replica of g3_lattice_fit
      with the orientation filter (bisection on eta, as stage C), with the gap s at eta 0.74 and 0.78 (replica; the binary's own value
      from the run header where the cell ran); the replica's feasibility verdict must equal the binary's in every grid cell of the
      three families, and the replica's default family must equal stage C's replica (stageC_tables.lat_fit) in every grid cell
  2.5 rule 4 identical and the M1/M2 harness outputs identical (stage H's H0 and H1 on the new build, item 1.2)
  ACCEPTANCE: item 1 PASS if 1.1-1.4 hold; item 2 PASS if 2.1-2.5 hold; else FAIL with the items that do not.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_d18_261010/d18_tables.py --bin-dir <bin dir> --out <evidence dir>
  python3 experiments_gen3_d18_261010/d18_tables.py --selftest      the replica against stage C's, and the orientation classifier
"""
import argparse, glob, hashlib, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
for sub in ("experiments_gen3_h_261009", "experiments_gen3_m4_261009"):
    sys.path.insert(0, os.path.join(HS, sub))
import d18_evidence as EV
import stageH_tables as HT
import stageC_tables as CT
PX, R_PX, TH_PX = CT.PX, CT.R_PX, CT.TH_PX
STAGEC_EV = os.path.join(HS, "experiments_gen3_m4_261009", "evidence_f39e485", "grid", "N400_eta0.70", "long")
P4_HASH = "bd25077c5f1c9c9f"


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, "rb").read() if p and os.path.exists(p) else None
def one(d, pat):
    g = sorted(glob.glob(os.path.join(d, pat))); return g[0] if len(g) == 1 else None
def log(d): return open(os.path.join(d, "run.log"), errors="ignore").read() if os.path.exists(os.path.join(d, "run.log")) else ""
def ec(d): return open(os.path.join(d, "exit_code.txt")).read().strip() if os.path.exists(os.path.join(d, "exit_code.txt")) else "?"
def health(d):
    m = re.search(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) .* hash=([0-9a-f]{16}) ", log(d), re.M)
    return (m.group(1), m.group(2)) if m else (None, None)
def lattice_s(d):
    m = re.search(r"^\[EDMD3-RUN\] .*? left=\[orient (\S+) rows (\S+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) g (\S+) s (\S+) n_vac (\d+)", log(d), re.M)
    return dict(orient=m.group(1), rows=m.group(2), ns=int(m.group(3)), nr=int(m.group(4)), s=float(m.group(8)), nvac=int(m.group(9))) if m else None
def load_row(d):
    p = one(d, "psi6_t_*_run0.csv")
    if not p: return None
    lines = open(p).read().split("\n"); hdr = lines[0].split(","); row = lines[1].split(",")
    return dict(zip(hdr, row))


# ---------------------------------------------------------------- the replica of g3_lattice_fit with the M4c orientation filter
def lat_cryst(o, ns, nr, a_s, a_r):
    along = 1 if ns < 2 else (0 if nr < 2 else (1 if a_r < 2.0 * a_s else 0))
    return 0 if along == (1 if o == 1 else 0) else 1


def lat_fit(n, R, Wx, Hy, cryst=None):
    best = -math.inf
    for o in (0, 1):
        Lr, Ls = (Hy, Wx) if o == 0 else (Wx, Hy)
        for alt in (0, 1):
            for ns in range(1, n + 1):
                nr = -(-n // ns)
                if alt:
                    if ns < 2: continue
                    nr = 2
                    while ((ns + 1) // 2) * nr + (ns // 2) * (nr - 1) < n: nr += 1
                s0 = CT.lat_eval(ns, nr, R, Lr, Ls, 0.0, alt)[0]
                if not s0 > 0: continue
                g = 0.0
                if math.isfinite(s0):
                    for _ in range(200):
                        gn = 0.5 * CT.lat_eval(ns, nr, R, Lr, Ls, g, alt)[0]
                        if not gn > 0: g = 0.0; break
                        if abs(gn - g) <= 1e-15 * (1 + gn): g = gn; break
                        g = gn
                else: g = min(0.5 * Ls - R, 0.5 * Lr - R)
                s, a_s, a_r, gw = CT.lat_eval(ns, nr, R, Lr, Ls, g, alt)
                if not (s > 0 and gw > 0): continue
                if cryst is not None and lat_cryst(o, ns, nr, a_s, a_r) != cryst: continue
                best = max(best, s)
    return best


def box(N, H, eta):
    L0 = N * math.pi / (8 * H * eta); boxW = 2 * L0 * PX; c = 0.5 * boxW
    return N // 2, R_PX, c - 0.5 * TH_PX, H * PX


def gap(N, H, eta, cryst=None):
    n, R, W, Hy = box(N, H, eta); return lat_fit(n, R, W, Hy, cryst)


def eta_max(N, H, cryst=None):
    lo, hi = 0.5, 0.9069
    if not gap(N, H, lo, cryst) > 0: return float("nan")
    if gap(N, H, hi, cryst) > 0: return hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if gap(N, H, mid, cryst) > 0: lo = mid
        else: hi = mid
    return lo


def selftest():
    ok = True
    for g in EV.grid():
        n, R, W, Hy = box(g["N"], g["H"], g["eta"])
        a, b = lat_fit(n, R, W, Hy), CT.lat_fit(n, R, W, Hy)
        e = (a == b) or (math.isinf(a) and math.isinf(b)); ok &= e
        if not e: print(f"selftest replica vs stage C at {g['name']}: {a} vs {b} **DIFFERENT**")
    print(f"selftest: the replica's default family = stage C's replica in all {len(EV.grid())} grid cells: {'ok' if ok else '**NO**'}")
    for o, ns, nr, a_s, a_r, want in ((0, 19, 11, 12.7, 45.4, 0), (0, 11, 19, 20.8, 24.0, 1), (1, 11, 19, 20.8, 24.0, 0), (1, 19, 11, 12.7, 45.4, 1),
                                      (0, 1, 50, 0.0, 24.0, 1), (1, 1, 50, 0.0, 24.0, 0)):
        got = lat_cryst(o, ns, nr, a_s, a_r); e = got == want; ok &= e
        print(f"selftest orientation: orient {'perpendicular' if o else 'parallel'} n_s {ns} n_r {nr} a_s {a_s} a_r {a_r}: {got} (expected {want}) {'ok' if e else '**WRONG**'}")
    print(f"SELFTEST {'PASS' if ok else '**FAIL**'}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir"); ap.add_argument("--out"); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: return selftest()
    B, O = os.path.abspath(a.bin_dir), os.path.abspath(a.out); NEW = os.path.join(B, "00ALLINONE")
    f1, f2 = [], []
    def item(lst, name, ok, txt):
        print(f"| {name} | {txt} | {'yes' if ok else '**NO**'} |")
        if not ok: lst.append(name)
    print("# Decision 18, items 1 and 2: the acceptance tables, printed by experiments_gen3_d18_261010/d18_tables.py (261012 sec. 4.7.42)\n")
    # ================= ITEM 1
    print("## Item 1: checkpoint/restart carries --gen3-psi6-compartments\n\n| item | result | holds |\n|---|---|---|")
    rec = open(os.path.join(B, "build_record.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec, re.M))
    ver = subprocess.run([NEW, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    names = ("00ALLINONE", "gen3_m1", "gen3_m2", "gen3_checkpoint_test", "gen3_checkpoint_test_ld")
    item(f1, "1.1 the binaries", all(want.get(n) == sha(os.path.join(B, n)) for n in names) and "dirty" not in ver,
         f"{ver}; {', '.join(n + (' = record' if want.get(n) == sha(os.path.join(B, n)) else ' **differs**') for n in names)}")
    rH = subprocess.run([sys.executable, os.path.join(HS, "experiments_gen3_h_261009", "stageH_tables.py"), "--bin-dir", B, "--record",
                         os.path.join(B, "build_record.txt"), "--out", os.path.join(O, "stageH")], capture_output=True, text=True)
    open(os.path.join(O, "stageH_tables_output.txt"), "w").write(rH.stdout + rH.stderr)
    acc = [l for l in rH.stdout.splitlines() if l.startswith("ACCEPTANCE (stage H)")]
    r4 = [l for l in rH.stdout.splitlines() if l.startswith("rule 4:")]
    h1 = [l for l in rH.stdout.splitlines() if l.startswith("| m") and "_output.txt" in l]
    hOK = bool(acc) and "PASS" in acc[0] and "NOT MET" not in acc[0]
    r4OK = bool(r4) and "IDENTICAL (both cases, both builds)" in r4[0]; h1OK = len(h1) == 3 and all("IDENTICAL" in l for l in h1)
    item(f1, "1.2 stage H's tables (unchanged script)", hOK,
         (acc[0][:120] if acc else "no ACCEPTANCE line") + f"; {r4[0] if r4 else 'no rule 4 line'}; H1: " + ", ".join(l.split("|")[1].strip() + " " + l.split("|")[2].strip() for l in h1))
    for cell, extra in (("B", EV.RB), ("C", EV.RC)):
        base = os.path.join(O, "ckpt", cell); U = os.path.join(base, "U")
        hdr = open(one(U, "psi6_t_*_run0.csv")).readline().strip() if one(U, "psi6_t_*_run0.csv") else ""
        pats = ["wall_x_positions_*_run0.csv", "psi6_t_*_run0.csv", "speed_of_sound_psi6.csv", "run_params.json"] + (
            ["virial_blocks_*_run0.csv", "snapshots_*_run0.bin"] if cell == "B" else [])
        hu, hashu, _, _ = HT.record(U)
        for d in sorted(glob.glob(os.path.join(base, "C_*")) + glob.glob(os.path.join(base, "R2_*"))):
            hx, hashx, ck, stop = HT.record(d)
            diff = [p.split("_*")[0] for p in pats if not (one(d, p) and rd(one(d, p)) == rd(one(U, p)))]
            ok = ec(d) == "0" and ec(U) == "0" and not diff and hx == hu and hashx == hashu and "psi6_orient_deg_0" in hdr
            item(f1, f"1.3 cell {cell} {os.path.basename(d)} vs U", ok,
                 f"files {'IDENTICAL' if not diff else 'DIFFERENT: ' + ', '.join(diff)}; run record {'=' if hx == hu else '**differs**'}; hash {hashx} (U {hashu}); "
                 f"checkpoint/restart lines {len(ck)}; exit {ec(d)}")
    for tag in ("on_ckpt_off_run", "off_ckpt_on_run"):
        d = os.path.join(O, "ckpt", "refusals", tag); st = [l for l in log(d).splitlines() if l.startswith("STOP")]
        item(f1, f"1.4 refusal: {tag}", ec(d) == "2" and any("--gen3-psi6-compartments" in l for l in st), f"exit {ec(d)}; {st[0] if st else 'no STOP line'}")
    # ================= ITEM 2
    print("\n## Item 2: M4c, the turned orientation\n\n| item | result | holds |\n|---|---|---|")
    md = os.path.join(O, "m4c", "ident")
    plist = dict((n.split("/", 1)[1], h) for h, n in re.findall(r"^([0-9a-f]{64})  (data_06e7248/\S+)$",
                                                                open(os.path.join(HS, "experiments_gen3_p4_261010", "P4_SHA256_FILES.txt")).read(), re.M))
    d = os.path.join(md, "p4a_A1"); cl, hs = health(d)
    pf = [(p, one(EV.P4A1, p)) for p in ("snapshots_*_run0.bin", "psi6_t_*_run0.csv", "virial_blocks_*_run0.csv")]
    orig_ok = all(o and plist.get(os.path.relpath(o, os.path.dirname(os.path.dirname(os.path.dirname(EV.P4A1))))) == sha(o) for _, o in pf)
    same_ok = all(rd(one(d, p)) == rd(o) for p, o in pf)
    item(f2, "2.1a P4a A seed 1, new build, default", hs == P4_HASH and cl == "1" and same_ok and orig_ok,
         f"hash {hs} (P4 {P4_HASH}); snapshots, psi6(t), virial blocks {'= P4' if same_ok else '**differ**'}; P4's files {'= P4_SHA256_FILES' if orig_ok else '**changed**'}")
    for tag, label in (("stageC", "2.1b stage C cell N400 eta0.70 (long)"), ("perp", "2.1c --gen3-lattice=perpendicular"), ("m4b", "2.1c M4b (P4 B seed 0)")):
        do, dn = os.path.join(md, f"{tag}_old"), os.path.join(md, f"{tag}_new")
        (co, ho), (cn, hn) = health(do), health(dn)
        init = (rd(os.path.join(do, "init.txt")) == rd(os.path.join(dn, "init.txt"))) if os.path.exists(os.path.join(do, "init.txt")) else None
        psi = rd(one(do, "psi6_t_*_run0.csv")) == rd(one(dn, "psi6_t_*_run0.csv"))
        txt = f"hash old {ho} new {hn}; psi6(t) {'IDENTICAL' if psi else 'DIFFERENT'}" + ("" if init is None else f"; init.txt {'IDENTICAL' if init else 'DIFFERENT'}")
        if tag == "stageC":
            hc = re.search(r" hash=([0-9a-f]{16}) ", log(STAGEC_EV)); hc = hc.group(1) if hc else None
            ic = rd(os.path.join(STAGEC_EV, "init.txt")) == rd(os.path.join(dn, "init.txt"))
            txt += f"; information, against stage C's evidence (f39e485): init.txt {'IDENTICAL' if ic else 'DIFFERENT'}, hash {hc}"
        item(f2, label, ho == hn and ho is not None and co == "1" and cn == "1" and psi and init is not False, txt)
    # 2.2, 2.3: the turned grid, the load orientations
    rows, okg, okl = [], True, True
    for g in EV.grid():
        gd = os.path.join(O, "m4c", "grid", g["name"]); L, S, Ot = (os.path.join(gd, x) for x in ("long", "same", "other"))
        infe = "[EDMD3-M4] INFEASIBLE" in log(L)
        m4c = re.search(r"^\[EDMD3-M4C\] .*; left (.*?), right (.*)$", log(L), re.M)
        r = dict(g=g, infeasible=infe, m4c=(m4c.group(1), m4c.group(2)) if m4c else None, lat=lattice_s(L))
        if infe:
            r["ok"] = ec(L) == "2" and ec(S) == "2" and ec(Ot) == "2"
        else:
            D = CT.read_dump(os.path.join(L, "init.txt")); x = D["P"][:, 0]
            r["eta"] = D["N"] * math.pi * D["R"] ** 2 / (D["boxW"] * D["boxH"]); r["nl"] = int((x < D["c"]).sum()); r["nr"] = D["N"] - r["nl"]
            r["gaps"] = CT.min_gaps(D)
            r["same"] = rd(os.path.join(L, "init.txt")) == rd(os.path.join(S, "init.txt")); r["other"] = rd(os.path.join(L, "init.txt")) != rd(os.path.join(Ot, "init.txt"))
            tm = re.search(r"\[EDMD3-TIES\] at the load: live events (\d+), sharing their time exactly with another live event (\d+), due at once (\d+)", log(L))
            r["ties"] = (int(tm.group(2)), int(tm.group(3))) if tm else None
            r["clean"], _ = health(L)
            lr = load_row(L); r["orient"] = (float(lr["psi6_orient_deg_0"]), float(lr["psi6_orient_deg_1"])) if lr else None
            r["ok"] = (ec(L) == "0" and ec(S) == "0" and ec(Ot) == "0" and abs(r["eta"] / g["eta"] - 1) < 1e-12 and r["nl"] == r["nr"] == g["N"] // 2 and min(r["gaps"]) > 0 and r["same"] and r["other"]
                       and r["ties"] == (0, 0) and r["clean"] == "1" and r["m4c"] == ("lines along y", "lines along y") and r["orient"] is not None
                       and min(abs(v) for v in r["orient"]) > 25)
        okg &= r["ok"]
        for fam in ("default", "zero"):
            ld = os.path.join(O, "m4c", "load", g["name"], fam)
            r[f"{fam}_infeasible"] = "[EDMD3-M4] INFEASIBLE" in log(ld)
            lr = load_row(ld)
            r[f"{fam}_orient"] = (float(lr["psi6_orient_deg_0"]), float(lr["psi6_orient_deg_1"])) if lr else None
            r[f"{fam}_lat"] = lattice_s(ld)
        rows.append(r)
    print("\n### 2.2 The turned grid (--gen3-lattice-orientation=turned), per cell\n")
    print("| cell | result | lattice (left): orient, rows, n_s x n_r, vacancies, gap s [sigma] | achieved / requested eta - 1 | N_L / N_R | smallest gap disk-disk / wall / "
          "divider [sigma] | same seed | other seed | ties at the load (shared, at once) | clean (400 sigma-time) | orientation at the load L / R [deg] | M4C |\n" + "|---" * 12 + "|")
    for r in rows:
        g = r["g"]; lt = r["lat"]
        latx = f"{lt['orient']}, {lt['rows']}, {lt['ns']} x {lt['nr']}, {lt['nvac']}, {lt['s'] / PX:.4f}" if lt else "-"
        if r["infeasible"]:
            print(f"| {g['name']} | INFEASIBLE (exit 2) | - | - | - | - | - | - | - | - | - | {r['m4c']} |"); continue
        print(f"| {g['name']} | {'ok' if r['ok'] else '**NO**'} | {latx} | {r['eta'] / g['eta'] - 1:+.1e} | {r['nl']} / {r['nr']} | "
              + " / ".join(f"{v:.4f}" for v in r["gaps"]) + f" | {'bit-identical' if r['same'] else '**DIFFERENT**'} | {'differs' if r['other'] else '**SAME**'} | "
              f"{r['ties']} | {r['clean']} | " + (" / ".join(f"{v:+.2f}" for v in r["orient"]) if r["orient"] else "-") + f" | {r['m4c']} |")
    item(f2, "2.2 the turned grid", okg, f"{sum(1 for r in rows if not r['infeasible'])} feasible cells pass every item, "
         f"{sum(1 for r in rows if r['infeasible'])} INFEASIBLE (listed); failing: {', '.join(r['g']['name'] for r in rows if not r['ok']) or 'none'}")
    print("\n### 2.3 The orientation at the load (the psi6(t) row at t = 0), L / R [deg]\n")
    print("| cell | default family | its lattice: orient, n_s x n_r | --gen3-lattice-orientation=0 | its lattice | turned |\n|---|---|---|---|---|---|")
    bad_d, bad_z, bad_t = [], [], []
    for r in rows:
        g = r["g"]; o = lambda v: " / ".join(f"{x:+.2f}" for x in v) if v else "INFEASIBLE"
        lat = lambda lt: f"{lt['orient']}, {lt['ns']} x {lt['nr']}" if lt else "-"
        if r["default_orient"] and max(abs(x) for x in r["default_orient"]) >= 5: bad_d.append(g["name"])
        if r["zero_orient"] and max(abs(x) for x in r["zero_orient"]) >= 5: bad_z.append(g["name"])
        if r.get("orient") and min(abs(x) for x in r["orient"]) <= 25: bad_t.append(g["name"])
        print(f"| {g['name']} | {o(r['default_orient'])} | {lat(r['default_lat'])} | {o(r['zero_orient'])} | {lat(r['zero_lat'])} | {o(r.get('orient'))} |")
    item(f2, "2.3a the orientation at the load (the code)", not (bad_z or bad_t),
         f"turned outside |o| > 25: {', '.join(bad_t) or 'none'}; 0 outside |o| < 5: {', '.join(bad_z) or 'none'}")
    print(f"| 2.3b the premise: the default family at 0 deg (reported, not counted) | default outside |o| < 5: {', '.join(bad_d) or 'none'} | "
          f"{'holds' if not bad_d else 'does not hold in the cells listed'} |")
    # 2.4 ceilings
    print("\n### 2.4 The feasibility ceiling per N (largest nominal eta with a positive gap; Python replica, bisection), and the gap s at eta 0.74 and 0.78 [sigma]\n")
    print("| N | H [sigma] | default: ceiling | 0: ceiling | turned: ceiling | default s at 0.74 / 0.78 | 0: s at 0.74 / 0.78 | turned: s at 0.74 / 0.78 (binary) |\n|---|---|---|---|---|---|---|---|")
    for N in EV.NS_LIST:
        H = 10.0 * math.sqrt(N / 100.0); cols = []
        for cr in (None, 0, 1): cols.append(f"{eta_max(N, H, cr):.4f}")
        ss = []
        for cr in (None, 0, 1):
            vals = [gap(N, H, e, cr) for e in (0.74, 0.78)]
            ss.append(" / ".join(f"{v / PX:.4f}" if v > 0 else "infeasible" for v in vals))
        tb = []
        for e in (0.74, 0.78):
            rr = [r for r in rows if r["g"]["N"] == N and abs(r["g"]["eta"] - e) < 1e-9][0]
            tb.append("infeasible" if rr["infeasible"] else (f"{rr['lat']['s'] / PX:.4f}" if rr["lat"] else "-"))
        print(f"| {N} | {H:g} | {cols[0]} | {cols[1]} | {cols[2]} | {ss[0]} | {ss[1]} | {ss[2]} ({' / '.join(tb)}) |")
    agree, n_agree = True, 0
    for r in rows:
        g = r["g"]
        for cr, feas in ((None, not r["default_infeasible"]), (0, not r["zero_infeasible"]), (1, not r["infeasible"])):
            n_agree += 1; agree &= (gap(g["N"], g["H"], g["eta"], cr) > 0) == feas
    sc = all((lambda a, b: a == b or (math.isinf(a) and math.isinf(b)))(lat_fit(*box(g["N"], g["H"], g["eta"])), CT.lat_fit(*box(g["N"], g["H"], g["eta"])))
             for g in EV.grid())
    item(f2, "2.4 the replica: verdict = the binary's; default = stage C's replica", agree and sc, f"{n_agree} cell-family verdicts {'all equal' if agree else '**not all equal**'}; "
         f"default family = stage C's replica in every cell: {'yes' if sc else '**no**'}")
    item(f2, "2.5 rule 4 and the M1/M2 harness outputs (stage H's H0 and H1 on the new build, item 1.2)", r4OK and h1OK,
         f"{r4[0] if r4 else 'no rule 4 line'}; H1: " + ", ".join(l.split("|")[1].strip() + " " + l.split("|")[2].strip() for l in h1))
    print(f"\nACCEPTANCE item 1: {'PASS' if not f1 else '**FAIL**: ' + '; '.join(f1)}")
    print(f"ACCEPTANCE item 2: {'PASS' if not f2 else '**FAIL**: ' + '; '.join(f2)}")


if __name__ == "__main__":
    main()
