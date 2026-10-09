#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.16): every table of stage C (M4), printed from the evidence folder of stageC_evidence.py.
The acceptance, per cell of the grid (N = 100, 400, 900, 1600 x eta = 0.10, pi/8, 0.60, 0.70, 0.716, 0.78, 0.85, 0.90):
  requested against achieved eta (from the dumped initial state: N pi R^2 / (boxW boxH)), with the legacy 1/48-grid box's eta next to it;
  N per compartment (disks left / right of the divider centre in the dump); the smallest initial surface gap disk-disk (all pairs),
  disk-wall and disk-divider [sigma]; the initial ties ([EDMD3-TIES]: live events sharing their time exactly, events due at once);
  same seed bit-identical (the dumps of "long" and "same"), different seeds different ("other"); health over the 400 sigma-time
  run (the run record, clean=1). Infeasible cells (the binary stops: no arrangement of the stated family has a positive gap) are
  listed with the feasibility limit of the family for that N and H, computed by a Python replica of g3_lattice_fit (information;
  the replica's feasibility verdict is checked against the binary's on every cell).
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m4_261009/stageC_tables.py --bin-dir <frozen binaries> --out <evidence dir> > stageC_tables_output.txt
"""
import argparse, hashlib, math, os, re, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import stageC_evidence as EV
MAIN_HS = EV.MAIN_HS; GATE = EV.GATE
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
REC = re.compile(r"\[EDMD3-HEALTH\] (.*?): clean=(\S+) (.*)")
R_PX, TH_PX, PX = 12.0, float(np.float32(0.05) * np.float32(24.0)), 24.0


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_dump(p):
    L = open(p).read().split("\n")
    assert L[0] == "gen3-initial v1"
    N, NL = int(L[1].split()[1]), int(L[1].split()[3])
    hx = lambda s: float.fromhex(s)
    boxW, boxH, R = hx(L[2].split()[1]), hx(L[3].split()[1]), hx(L[4].split()[1])
    c, th = hx(L[5].split()[1]), hx(L[5].split()[2])
    P = np.array([[hx(v) for v in L[6 + i].split()] for i in range(N)])
    return dict(N=N, NL=NL, boxW=boxW, boxH=boxH, R=R, c=c, th=th, P=P)


def min_gaps(D):
    x, y, R = D["P"][:, 0], D["P"][:, 1], D["R"]
    dd = np.inf
    for i0 in range(0, len(x), 512):                    # all pairs, in blocks
        dx = x[i0:i0 + 512, None] - x[None, :]; dy = y[i0:i0 + 512, None] - y[None, :]
        r = np.sqrt(dx * dx + dy * dy); idx = np.arange(i0, min(i0 + 512, len(x)))
        r[np.arange(len(idx)), idx] = np.inf
        dd = min(dd, r.min())
    wall = min((x - R).min(), (D["boxW"] - R - x).min(), (y - R).min(), (D["boxH"] - R - y).min())
    div = (np.abs(x - D["c"]) - 0.5 * D["th"] - R).min()
    return (dd - 2 * R) / PX, wall / PX, div / PX


# ---------------------------------------------------------------- a Python replica of g3_lattice_fit (information)
def lat_eval(ns, nr, R, Lr, Ls, g, alt=0):
    us, ur = Ls - 2 * (R + g), Lr - 2 * (R + g)
    if not (us >= 0 and ur >= 0): return -math.inf, 0, 0, 0
    if alt and (ns < 2 or nr < 2): return -math.inf, 0, 0, 0
    a_s = us / (ns - 1) if ns > 1 else 0.0
    if alt: a_r = ur / (nr - 1)
    else: a_r = ur / ((nr - 1) + (0.5 if ns > 1 else 0.0)) if nr > 1 else (2 * ur if ns > 1 else 0.0)
    d = math.inf
    if nr > 1: d = a_r
    if ns > 1: d = min(d, math.sqrt(a_s * a_s + 0.25 * a_r * a_r))
    if ns > 2: d = min(d, 2 * a_s)
    gw = g
    if ns == 1: gw = min(g, 0.5 * Ls - R)
    if ns == 1 and nr == 1: gw = min(0.5 * Ls - R, 0.5 * Lr - R)
    return d - 2 * R, a_s, a_r, gw


def lat_fit(n, R, Wx, Hy):
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
                s0 = lat_eval(ns, nr, R, Lr, Ls, 0.0, alt)[0]
                if not s0 > 0: continue
                g = 0.0
                if math.isfinite(s0):
                    for _ in range(200):
                        gn = 0.5 * lat_eval(ns, nr, R, Lr, Ls, g, alt)[0]
                        if not gn > 0: g = 0.0; break
                        if abs(gn - g) <= 1e-15 * (1 + gn): g = gn; break
                        g = gn
                else: g = min(0.5 * Ls - R, 0.5 * Lr - R)
                s, _, _, gw = lat_eval(ns, nr, R, Lr, Ls, g, alt)
                if s > 0 and gw > 0: best = max(best, s)
    return best


def feasible(N, H, eta):
    L0 = N * math.pi / (8 * H * eta); boxW = 2 * L0 * PX; c = 0.5 * boxW
    W = c - 0.5 * TH_PX
    return lat_fit(N // 2, R_PX, W, H * PX) > 0


def eta_max(N, H):
    lo, hi = 0.5, 0.9069                                    # feasible at 0.5 for every N here; not above close packing
    if feasible(N, H, hi): return hi
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if feasible(N, H, mid): lo = mid
        else: hi = mid
    return lo


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out)
    print("# Stage C (M4) evidence tables, printed by experiments_gen3_m4_261009/stageC_tables.py\n")
    print("## 1. The frozen binaries (rule 5)\n")
    rec_txt = open(os.path.join(HERE, "build_record_f39e485.txt")).read()
    want = dict((n, h) for h, n in re.findall(r"^\s+([0-9a-f]{64})\s+(\S+)$", rec_txt.split("## SHA-256")[1], re.M))
    print("| binary | SHA-256 (frozen copy) | = build record |\n|---|---|---|")
    for n in ("00ALLINONE", "gen3_m1", "gen3_m2"):
        got = sha(os.path.join(B, n)); print(f"| {n} | {got} | {'yes' if got == want.get(n) else '**NO**'} |")
    v = subprocess.run([os.path.join(B, "00ALLINONE"), "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    print(f"\n`00ALLINONE --version`: {v} ({'no -dirty' if '-dirty' not in v else '**-dirty**'})")
    print("\n## 2. Rule 4: gen2 byte identity with 7b08827 after the M4 change (ctrl_min, ctrl_leg)\n")
    import shutil
    for tag in ("default", "gen2flag"):
        for case in ("ctrl_min", "ctrl_leg"):
            dst = os.path.join(O, "e0", tag, case, "ref")
            if not os.path.exists(dst): shutil.copytree(os.path.join(E0REF, case, "ref"), dst)
        r = subprocess.run([sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "report", "--out", os.path.join(O, "e0", tag)],
                           capture_output=True, text=True).stdout
        open(os.path.join(O, "e0", f"report_{tag}.txt"), "w").write(r)
        rows = [l for l in r.splitlines() if l.startswith("| ctrl_") and "git " in l]
        print(f"{'default build' if tag == 'default' else '--engine=gen2'}: " + "; ".join(f"{l.split('|')[1].strip()}: audit vs plain {l.split('|')[-3].strip()}, plain vs ref {l.split('|')[-2].strip()}" for l in rows))
    print("\n## 3. The engine change (edmd3_tie_stats, read-only): the M1 and M2 harness outputs, byte for byte\n")
    print("| output | committed output | identical |\n|---|---|---|")
    for f, ref in (("m1_audit_output.txt", "experiments_gen3_m1_261008/m1_audit_output.txt"),
                   ("m2_audit_quick_output.txt", "experiments_gen3_m2_261009/m2_audit_quick_output.txt"),
                   ("m2_audit_output.txt", "experiments_gen3_m2_261009/m2_audit_output.txt")):
        p, q = os.path.join(O, "harness", f), os.path.join(HS, ref)
        print(f"| {f} | {ref} | {'IDENTICAL' if os.path.exists(p) and open(p, 'rb').read() == open(q, 'rb').read() else '**DIFFERENT or missing**'} |")
    print("\n## 4. The acceptance grid (gate dense construction H = 10 sqrt(N/100), L0 = N pi / (8 H eta) exactly; M4 lattice, jitter 0.25; "
          "400 sigma-time with the divider held)\n")
    print("| N | eta | H | L0 [sigma] | lattice left (orient, n_s x n_r, n_vac, s [px], jitter [px]) | achieved eta - requested | eta of the legacy "
          "1/48-grid box - requested | N_L / N_R | min gap disk-disk / wall / divider [sigma] | ties at the load (tied, at once, live) | same seed "
          "identical | other seed different | clean (400 sigma-time) | events | run [s] |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    n_ok = n_inf = n_bad = 0; replica_mismatch = []
    for c in EV.cells():
        d0 = os.path.join(O, "grid", f"N{c['N']}_eta{c['lab'].replace('/', '')}")
        lg = os.path.join(d0, "long", "run.log"); t = open(lg, errors="ignore").read() if os.path.exists(lg) else ""
        inf = "[EDMD3-M4] INFEASIBLE" in t
        rep = feasible(c["N"], c["H"], c["eta"])
        if inf == rep: replica_mismatch.append(c)          # the replica says feasible where the binary stopped, or the reverse
        if inf:
            n_inf += 1
            print(f"| {c['N']} | {c['lab']} | {c['H']:g} | {c['L0']:.6f} | **INFEASIBLE** (the binary stops: no arrangement of the stated family "
                  f"has a positive gap) | | | | | | | | | | |"); continue
        m = re.search(r"left=\[orient (\S+) rows (\S+) n_s (\d+) n_r (\d+) a_s \S+ a_r \S+ g \S+ s (\S+) n_vac (\d+) jitter (\S+) px\]", t)
        lat = f"{m.group(1)} ({m.group(2)} rows), {m.group(3)} x {m.group(4)}, {m.group(6)}, {float(m.group(5)):.4g}, {float(m.group(7)):.3g}" if m else "?"
        D = read_dump(os.path.join(d0, "long", "init.txt"))
        eta_a = D["N"] * math.pi * D["R"] ** 2 / (D["boxW"] * D["boxH"])
        sw = math.floor(2 * np.float32(c["L0"]) * np.float32(24.0))                     # the legacy SIM_WIDTH (float L0, int cast)
        eta_g = D["N"] * math.pi * D["R"] ** 2 / (sw * D["boxH"])
        nl = int((D["P"][:, 0] < D["c"]).sum()); nr = D["N"] - nl
        gd, gwl, gdv = min_gaps(D)
        tm = re.search(r"\[EDMD3-TIES\] at the load: live events (\d+), sharing their time exactly with another live event (\d+), due at once (\d+)", t)
        ties = f"{tm.group(2)}, {tm.group(3)}, {tm.group(1)}" if tm else "no line"
        same = open(os.path.join(d0, "long", "init.txt"), "rb").read() == open(os.path.join(d0, "same", "init.txt"), "rb").read()
        other = open(os.path.join(d0, "long", "init.txt"), "rb").read() != open(os.path.join(d0, "other", "init.txt"), "rb").read()
        r = REC.findall(t); clean = r[-1][1] if r else "no record"
        kv = dict(re.findall(r"(\w+)=(\S+)", r[-1][2])) if r else {}
        ev = sum(int(kv.get(f, 0)) for f in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")) if kv else 0
        ok = bool(abs(eta_a - c["eta"]) <= 1e-12 * c["eta"] and nl == D["NL"] == c["N"] // 2 and nr == c["N"] - c["N"] // 2 and min(gd, gwl, gdv) > 0
              and tm and tm.group(2) == "0" and tm.group(3) == "0" and same and other and clean == "1")
        n_ok += ok; n_bad += not ok
        print(f"| {c['N']} | {c['lab']} | {c['H']:g} | {c['L0']:.6f} | {lat} | {eta_a - c['eta']:+.1e} | {eta_g - c['eta']:+.1e} | {nl} / {nr} | "
              f"{gd:.3g} / {gwl:.3g} / {gdv:.3g} | {ties} | {'yes' if same else '**NO**'} | {'yes' if other else '**NO**'} | {clean} | {ev} | "
              f"{kv.get('run_s', '?')} |")
    print(f"\ncells: {n_ok} pass every item, {n_bad} fail an item, {n_inf} infeasible (no run); the Python replica of the fit agrees with the "
          f"binary's feasibility on {len(EV.cells()) - len(replica_mismatch)} of {len(EV.cells())} cells")
    print("\n## 5. The feasibility limit of the stated lattice family (information; the Python replica, bisection on eta)\n")
    print("| N | H [sigma] | largest nominal eta with a positive gap |\n|---|---|---|")
    for N in EV.NS_LIST:
        H = 10.0 * math.sqrt(N / 100.0); print(f"| {N} | {H:g} | {eta_max(N, H):.4f} |")


if __name__ == "__main__":
    main()
