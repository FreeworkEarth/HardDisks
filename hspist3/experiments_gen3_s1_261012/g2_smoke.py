#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.10, GATE G2): THE CELL SMOKE ON THE MAC. "Every cell of
B, D and S (both starts): one seed, lightest mass, hold 50 sigma-time, 5 periods. Printed per cell: recorded eta and lengths equal the
requested ones to 1e-12, N_L = N_R, ties at the load = 0, the orientation, clean=1."
  the lines   per cell of the registered task list (tasks_stage1_koa_261012.txt) its line of the lightest mass (alpha 0.5 in B and D, 1 in
              S) and s = 0: 110 lines; the cell's registered L0, H, N_s, stride, Delta_psi, T_snap and start, unchanged;
  the worker  the REGISTERED worker cluster/gen3_stage1_261012/s1_worker_koa.sh, copied with exactly two changes (each found exactly once,
              printed): --wall-hold-steps=$hold -> --wall-hold-steps=3000 (50 sigma-time) and --target-oscillations=200 ->
              --target-oscillations=5; everything else -- the command, the health rule, the file moves -- is the registered worker's;
  the binary  the clean Mac build of f72d764 (experiments_gen3_d18_261010/bin_f72d764, its build record's SHA-256 checked first);
  the checks  from each run's log: clean=1; the run record's exact_box = 1, its eta against the requested eta (relative <= 1e-12), its box
              width / 24 against 2 L0 and its height / 24 against H (relative <= 1e-12), N_L = N_R = N_s; [EDMD3-TIES] at the load: 0;
              the orientation (the run record's left/right fields, and from them the crystal orientation by g3_lat_cryst of f72d764, M4c:
              close-packed lines along x or along y -- the registration's "orientation the default chose", per cell; in block S the
              [EDMD3-M4C] line must name the requested start and the derived orientation must be it);
              and, beyond the gate's list, the binary's f_pred (printed %.8g) against the generator's replication rounded to 8 digits, and
              its planned steps against max(10000, ceil(5 / nu_pred / dt)) (the command's --oscillation-min-steps); the files the worker moved (trace, psi6, virial; the snapshot file where T_snap > 0).
SDL's dummy driver; at most 8 processes; --pilot K runs the first K lines only (the data estimate before the full run).
usage (from hspist3/): python3 experiments_gen3_s1_261012/g2_smoke.py --out <data dir> [--pilot K] > experiments_gen3_s1_261012/g2_smoke_output.txt
"""
import argparse, glob, hashlib, math, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage1_261012")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import gen_stage1_koa_261012 as GEN
P = os.path.join(HS, "cluster", "gen3_stage1_261012")
BIN = os.path.join(HS, "experiments_gen3_d18_261010", "bin_f72d764", "00ALLINONE")
BIN_SHA = "d0a8c782c36269c151b64bc264d0c6c4b92b06eef71b41d3fec30c35264b7a07"
CHANGES = [("--wall-hold-steps=$hold", "--wall-hold-steps=3000"), ("--target-oscillations=200", "--target-oscillations=5")]
SIDE = re.compile(r"(left|right)=\[orient (parallel|perpendicular) rows (alternate|equal) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) ")


def cryst(orient, ns, nr, a_s, a_r):
    """f72d764 00ALLINONE.c's g3_lat_cryst (M4c): 0 = close-packed lines along x, 1 = along y (the turned crystal)"""
    along_rows = 1 if ns < 2 else (0 if nr < 2 else int(a_r < 2.0 * a_s))
    return 0 if along_rows == int(orient == "perpendicular") else 1


RUN = re.compile(r"^\[EDMD3-RUN\] .*?box=(\S+) x (\S+) px exact_box=(\d) eta=(\S+) N=(\d+) N_L=(\d+) N_R=(\d+) seeding=(\S+) left=\[orient (\S+ \S+ \S+)[^]]*\] right=\[orient (\S+ \S+ \S+)", re.M)


def lines():
    reg = GEN.registered() if hasattr(GEN, "registered") else None
    L = open(os.path.join(P, "tasks_stage1_koa_261012.txt")).read().splitlines(); out, seen = [], set()
    rows, _ = GEN.plan(); meta = {(r["cid"], r["M"]): r for r in rows}
    light = {}
    for r in rows: light[r["cid"]] = min(light.get(r["cid"], 10 ** 9), r["M"])
    for l in L:
        f = l.split(); cid = f[1].split("/")[2]; M = int(f[2]); s = int(f[3])
        if cid in seen or M != light[cid] or s != 0: continue
        seen.add(cid); out.append((l, meta[(cid, M)]))
    return out


def run(args):
    i, line, worker, out = args
    env = dict(os.environ, HD_BIN=BIN, HD_DATA=out, SDL_VIDEODRIVER="dummy", SDL_RENDER_DRIVER="software")
    t = time.time()
    p = subprocess.run(["perl", "-e", "alarm 1800; exec @ARGV", "bash", worker] + line.split(), env=env, capture_output=True, text=True)
    return i, p.returncode, time.time() - t, (p.stdout + p.stderr).strip()[-300:]


def check(out, line, m):
    f = line.split(); rel, M = f[1], int(f[2]); d = os.path.join(out, rel)
    log = open(os.path.join(d, "run.log"), errors="replace").read() if os.path.exists(os.path.join(d, "run.log")) else ""
    r = dict(cid=m["cid"], M=M, start=m["start"])
    h = re.search(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", log, re.M); r["clean"] = h.group(1) if h else "-"
    g = RUN.search(log)
    if not g: r["ok"] = False; r["why"] = "no run record"; return r
    W, Hp, xb, eta, N, NL, NR, seed, ol, orr = g.groups()
    W, Hp, eta = float(W), float(Hp), float(eta)
    r["eta_rel"] = abs(eta - m["eta"]) / m["eta"]; r["w_rel"] = abs(W / 24 - 2 * m["L0"]) / (2 * m["L0"]); r["h_rel"] = abs(Hp / 24 - m["H"]) / m["H"]
    r["NLR"] = (int(NL), int(NR)); r["exact"] = int(xb); r["orient"] = f"{ol} / {orr}"
    sides = {m_.group(1): cryst(m_.group(2), int(m_.group(4)), int(m_.group(5)), float(m_.group(6)), float(m_.group(7)))
             for m_ in SIDE.finditer(re.search(r"^\[EDMD3-RUN\] .*$", log, re.M).group(0))}
    r["cryst"] = (sides.get("left"), sides.get("right"))
    t = re.search(r"^\[EDMD3-TIES\] at the load: live events (\d+), sharing their time exactly with another live event (\d+)", log, re.M)
    r["ties"] = int(t.group(2)) if t else -1
    c = re.search(r"^\[EDMD3-M4C\] lattice orientation requested: (\S+) .*?\); (.*)$", log, re.M)
    r["m4c"] = f"{c.group(1)}: {c.group(2)}" if c else "-"
    fp = re.search(r"f_pred=(\S+)\s+samples=(\d+)", log)
    r["fpred_rel"] = abs(float(fp.group(1)) - m["nup"]) / m["nup"] if fp else math.nan
    r["fpred_str"] = fp is not None and fp.group(1) == f"{m['nup']:.8g}"          # the binary prints %.8g: the replication to 8 digits
    r["steps_ok"] = fp is not None and int(fp.group(2)) == max(10000, math.ceil((5.0 / m["nup"]) / GEN.DTF))   # --oscillation-min-steps=10000
    have = lambda pat: bool(glob.glob(os.path.join(d, pat)))
    r["files"] = (have("wall_x_positions_*_run0.csv"), have("psi6_t_*_run0.csv.gz"), have("virial_blocks_*_run0.csv.gz"), have("snapshots_*_run0.bin.gz"))
    snap_need = float(f[12]) > 0
    r["ok"] = (r["clean"] == "1" and r["exact"] == 1 and r["eta_rel"] <= 1e-12 and r["w_rel"] <= 1e-12 and r["h_rel"] <= 1e-12 and r["NLR"] == (m["Ns"], m["Ns"])
               and r["ties"] == 0 and (m["start"] == "default" or (r["m4c"].startswith(m["start"] + ":") and r["cryst"] == ((0, 0) if m["start"] == "0" else (1, 1)))) and r["fpred_str"] and r["steps_ok"]
               and all(r["files"][:3]) and (r["files"][3] or not snap_need))
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--pilot", type=int, default=0); ap.add_argument("--check-only", action="store_true"); a = ap.parse_args()
    if hashlib.sha256(open(BIN, "rb").read()).hexdigest() != BIN_SHA: sys.exit("STOP: the f72d764 Mac binary is not the recorded one")
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    w = open(os.path.join(P, "s1_worker_koa.sh")).read()
    for old, new in CHANGES:
        if w.count(old) != 1: sys.exit(f"STOP: '{old}' occurs {w.count(old)} times in the registered worker")
        w = w.replace(old, new)
    worker = os.path.join(out, "s1_worker_smoke.sh"); open(worker, "w").write(w)
    L = lines()
    if a.pilot: L = L[:a.pilot]
    print("# Gate G2, the cell smoke on the Mac (sec. 4.8.10), printed by experiments_gen3_s1_261012/g2_smoke.py\n")
    print(f"binary: the clean Mac build of f72d764 (SHA-256 {BIN_SHA[:16]}...); the registered worker with exactly these changes: "
          + "; ".join(f"'{o}' -> '{n}'" for o, n in CHANGES) + f"; {len(L)} cells (the lightest mass, s = 0)\n")
    t0 = time.time(); bad = []
    if a.check_only:
        print("--check-only: the runs of an earlier call are checked again (nothing runs)")
    else:
        with ThreadPoolExecutor(max_workers=8) as ex:
            res = list(ex.map(run, [(i, l, worker, out) for i, (l, m) in enumerate(L)]))
        bad = [x for x in res if x[1] != 0]
        print(f"{len(res)} runs in {time.time() - t0:.0f} s; worker exit codes other than 0: {len(bad)} " + ("; ".join(f'{L[x[0]][1]["cid"]}: {x[3]}' for x in bad[:5]) if bad else ""))
    print("\n| cell | M | start | clean | exact box | eta rel. diff | width rel. diff | height rel. diff | N_L, N_R | ties at the load | orientation (run record) | "
          "[EDMD3-M4C] | f_pred rel. diff | planned steps | files: trace, psi6, virial, snapshot | ok |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    R = [check(out, l, m) for l, m in L]
    for r in R:
        if "why" in r: print(f"| {r['cid']} | {r['M']} | {r['start']} | {r['clean']} | {r['why']} | | | | | | | | | | | **FAIL** |"); continue
        print(f"| {r['cid']} | {r['M']} | {r['start']} | {r['clean']} | {r['exact']} | {r['eta_rel']:.1e} | {r['w_rel']:.1e} | {r['h_rel']:.1e} | {r['NLR'][0]}, {r['NLR'][1]} | "
              f"{r['ties']} | {r['orient']}; crystal: {' / '.join('lines along x (0 deg)' if c == 0 else ('lines along y (turned)' if c == 1 else '-') for c in r['cryst'])} | {r['m4c']} | {r['fpred_rel']:.1e} ({'= the replication to 8 digits' if r['fpred_str'] else '**differs**'}) | "
              f"{'= max(10000, ceil(5/nu_pred/dt))' if r['steps_ok'] else '**differs**'} | "
              f"{', '.join('yes' if x else 'no' for x in r['files'])} | {'PASS' if r['ok'] else '**FAIL**'} |")
    print(f"\nG2: {'PASS' if all(r['ok'] for r in R) and not bad else 'FAIL'} ({sum(r['ok'] for r in R)} of {len(R)} cells)")
    du = subprocess.run(["du", "-sk", out], capture_output=True, text=True).stdout.split()[0]
    print(f"data: {int(du) / 1024:.1f} MB in {out} ({int(du) / 1024 / max(1, len(R)):.2f} MB per cell)")


if __name__ == "__main__":
    main()
