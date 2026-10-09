#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.20; stage G): PREPARED, NOT RUN. Parts (c) and (d) of g3_gate_koa.sbatch.
(c) the KOA rate table of gen3: N = 100, 400, 900, 1600 x eta = pi/8, 0.70, 0.78; the gate's dense construction H = 10 sqrt(N/100),
    L0 = N pi / (8 H eta) exact (--gen3-exact-box), M4 seeding, 500 sigma-time held, 100 sigma-time released (M = 300), one run per
    cell, ONE AT A TIME (the rate is a timing): events, run_s, engine_s, driver share, events/s from each run record.
(d) the long-double spot check: the koa-ld binary against the koa binary, same commands and seeds: N = 100 at pi/8 (M = 300) and
    N = 400 at eta 0.70 (M = 500), 4 seeds each, the T-prime protocol (200 predicted periods, M4 seeding, exact box); per run: clean,
    events, the event hash, the largest contact gap and gap / u_t (the long-double quantum is 2^-50, so gap / u_t is quoted against
    each build's own u_t, the engine's value in its [EDMD3-GAP] line), and the argmax nu (reduce_B.py's estimator, information).
usage (from ~/harddisks_gen3/hspist3, inside the job): python3 cluster/gen3_koa_261009/g3_rate_ld_koa.py --bin B --bin-ld B_LD --out OUT [--jobs 8]
"""
import argparse, math, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import numpy as np
BASE = 20261219      # never used before (stage G)


def sos(binp, d, N, eta, M, seed, hold, record=None):
    H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta); ns = N // 2
    nu = T.kr_cs(eta) * T.x_of(M, L0)
    c = [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
         f"--particles={N}", f"--particles-boxes={ns},{ns}", f"--height={H:.17g}", "--particle-radius=0.5", "--wall-thickness=0.05",
         "--wall-thickness-vis=0.05", f"--lengths={L0:.17g}", f"--wall-masses={M}", "--repeats=1", f"--seed={BASE}",
         f"--wall-hold-steps={hold}", "--fixed-dt=0.4", "--oscillation-safety=1.0", "--oscillation-min-steps=10",
         "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={T.d_stride(nu)}", f"--speed-sound-run-dir={d}",
         f"--speed-sound-exact-seed={seed}", "--engine=gen3", "--gen3-seeding=lattice", "--gen3-exact-box"]
    return c + ([f"--record-sigma-time={record}"] if record else ["--target-oscillations=200"])


def run(cmd, d):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cmd) + "\n")
    return subprocess.Popen(cmd, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT, env=dict(os.environ, HD_CONTACT_AUDIT="1"))


def record(d):
    t = open(os.path.join(d, "run.log"), errors="ignore").read()
    m = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) (.*)$", t, re.M)
    if not m: return None
    kv = dict(re.findall(r"(\w+)=(\S+)", m[-1][1])); kv["clean"] = m[-1][0]
    kv["events"] = sum(int(kv[f]) for f in ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band"))
    g = re.search(r"\[EDMD3-GAP\] .*?max gap \[px\] pair (\S+), wall (\S+), divider (\S+)", t)
    kv["gap"] = max(float(x) for x in g.groups()) if g else math.nan
    u = re.search(r"\[EDMD3-GAP\] .*?time quantum u_t (\S+) units", t); kv["u_t"] = float(u.group(1)) if u else math.nan   # the engine's own
    return kv


def argmax_nu(d):
    import glob
    fs = glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv"))
    if not fs: return math.nan
    from paper1_populate_cs_err_20261002 import TD, X_EDGE
    t, x, nup = T._load(fs[0]); dt = (t[-1] - t[0]) / (len(t) - 1); n = T._prefix(t, nup, TD)
    if n is None: return math.nan
    P, df = T._spectrum(x[:n], dt); k = int(round(TD / X_EDGE)); return (k + int(np.argmax(P[k:]))) * df


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--bin-ld", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=8); a = ap.parse_args()
    O = os.path.abspath(a.out)
    print("== (c) the KOA rate table of gen3 (one run at a time)\n")
    print("| N | eta | events | run [s] | engine [s] | driver share | events/s | clean |\n|---|---|---|---|---|---|---|---|")
    k = 0
    for N in (100, 400, 900, 1600):
        for lab, eta in (("pi/8", math.pi / 8), ("0.70", 0.70), ("0.78", 0.78)):
            d = os.path.join(O, "rate", f"N{N}_{lab.replace('/', '')}")
            p = run(sos(a.bin, d, N, eta, 300, T.run_seed(BASE, 0, k, 0), 30000, 100), d); p.wait(); k += 1
            r = record(d)
            if not r: print(f"| {N} | {lab} | no run record (exit {p.returncode}) | | | | | |"); continue
            print(f"| {N} | {lab} | {r['events']} | {float(r['run_s']):.1f} | {float(r['engine_s']):.1f} | {float(r['driver_share']):.3f} | "
                  f"{r['events'] / float(r['run_s']):.3g} | {r['clean']} |", flush=True)
    print("\n== (d) the long-double spot check (koa-ld against koa, same commands and seeds)\n")
    jobs = []
    for N, eta, M in ((100, math.pi / 8, 300), (400, 0.70, 500)):
        for r in range(4):
            s = T.run_seed(BASE, 1, N, r)
            for tag, b in (("koa", a.bin), ("koa-ld", a.bin_ld)):
                jobs.append((N, eta, M, r, tag, os.path.join(O, "ld", f"N{N}_M{M}", f"seed{r}_{tag}"), sos(b, None, N, eta, M, s, 2000)))
    running = []
    for N, eta, M, r, tag, d, c in jobs:
        while len(running) >= a.jobs:
            running = [x for x in running if x.poll() is None]; time.sleep(0.5)
        c = [x if not x.startswith("--speed-sound-run-dir=") else f"--speed-sound-run-dir={d}" for x in c]
        running.append(run(c, d))
    for x in running: x.wait()
    print("| N | M | seed | build | clean | events | event hash | max contact gap [px] | u_t | gap / u_t [px/unit] | argmax nu |\n|---|---|---|---|---|---|---|---|---|---|---|")
    for N, eta, M, r, tag, d, c in jobs:
        q = record(d)
        if not q: print(f"| {N} | {M} | {r} | {tag} | no run record | | | | | | |"); continue
        print(f"| {N} | {M} | {r} | {tag} | {q['clean']} | {q['events']} | {q['hash']} | {q['gap']:.3g} | {q['u_t']:.3g} | {q['gap'] / q['u_t']:.3g} | "
              f"{argmax_nu(d):.6f} |")


if __name__ == "__main__":
    main()
