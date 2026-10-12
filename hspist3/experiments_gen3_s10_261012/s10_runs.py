#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (decision 38 item 5; Report S1-0, gate G0 of the recorded sec. 4.8.10, items b, c and f): the code facts
checked by running the code, on the Mac, with the clean Mac build of f72d764 (experiments_gen3_d18_261010/bin_f72d764, its build
record's SHA-256 is checked first). Nothing here is a stage 1 trajectory: test seeds, three small cells, short probes.

  B1  block B, every flag of sec. 4.8.3.2 at once: N = 100, eta = 0.30 (L0 off the 1/48 grid), alpha = 0.5, the registered hold
      (T_eq = 1000 sigma-time = 60,000 steps) and record (200 predicted periods), psi6 per compartment every Delta_psi, virial
      blocks of 50 sigma-time, snapshots every (T_eq + T_rec)/100; HD_KE_TRACE=1 and HD_CONTACT_AUDIT=1 as in Test G's B command.
  S0, ST  block S, the same with --gen3-lattice-orientation=0 and =turned: N = 100, eta = 0.60, alpha = 1.
  AF1-AF4  Part A's held path (the worker's AF command, conf_worker.sh) under gen3, 600 held steps and a 60-step tail:
      AF1 the registered start asked for (--gen3-seeding=lattice), AF2 the exact box asked for (--gen3-exact-box),
      AF3 the legacy start at the eta = 0.30 cell (L0 off the grid), divider at L0,
      AF4 the legacy start at an L0 on the 1/48 grid (13.083333..., 628 px wide), divider off the 1/24 grid (L0 + 0.0123 sigma).
  AF5, AF6  the dense end of Part A's grid, L0 on the 1/48 grid next to the eta = 0.69 cell (eta 0.6905): AF5 N = 100 (L0 = 273/48),
      the held path's start there; AF6 N = 1600 (L0 = 1092/48), 3000 held steps (50 sigma-time), the event log's rate per kind.
  AF5v, AF6v  AF5 and AF6 without --quiet and with 60 held steps, only to read the seeder's own lines (--quiet hides them).
Every run: SDL dummy driver, at most 3 at once, a time limit of 600 s each.
usage (from hspist3/): python3 experiments_gen3_s10_261012/s10_runs.py      -> data_f72d764/<case>/ (the tables: s10_tables.py)"""
import hashlib, math, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

BIN = os.path.join(HS, "experiments_gen3_d18_261010", "bin_f72d764", "00ALLINONE")
BIN_SHA = "d0a8c782c36269c151b64bc264d0c6c4b92b06eef71b41d3fec30c35264b7a07"   # its build_record.txt
OUT = os.path.join(HERE, "data_f72d764")
DT = 0.4 / 24.0
ENV = dict(os.environ, SDL_VIDEODRIVER="dummy", SDL_RENDER_DRIVER="software")


def nu_pred_binary(N, H, L0, M):
    """the binary's speed_sound_henderson_frequency_estimate (f72d764 00ALLINONE.c:16019-16039) as the B loop calls it (:17059-17065)
    under --gen3-exact-box: eta from L0 exactly, L_eff = (float) L0 - 1"""
    import numpy as np
    eta = N * math.pi * 0.25 / (2.0 * L0 * H)
    a = 0.125; num = 1.0 + eta + 3.0 * a * eta * eta - a * eta ** 3
    cs = math.sqrt(2.0 * num / (1.0 - eta) ** 3)
    return cs * T.k_root(M / N) / (2.0 * math.pi * (float(np.float32(L0)) - 1.0))


def b_case(name, N, eta, alpha, seed, extra):
    H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta); Ns = N // 2; M = int(round(alpha * N))
    Teq = 1000.0 if eta <= 0.60 + 1e-12 else 5000.0
    nup = nu_pred_binary(N, H, L0, M); Trec = 200.0 / nup
    dpsi = max(1.0, Trec / 2000.0); Tsnap = (Teq + Trec) / 100.0
    nukr = T.kr_cs(eta) * T.k_root(M / N) / (2.0 * math.pi * (L0 - 1.025)); stride = T.d_stride(nukr)
    d = os.path.join(OUT, name)
    argv = [BIN, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={N}", f"--particles-boxes={Ns},{Ns}", f"--height={H:.17g}", "--particle-radius=0.5", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", f"--lengths={L0:.17g}", f"--wall-masses={M}", "--repeats=1", f"--seed={seed}",
            f"--wall-hold-steps={int(round(60 * Teq))}", "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={stride}",
            f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}", "--engine=gen3", "--gen3-seeding=lattice", "--gen3-exact-box",
            "--gen3-psi6-compartments", f"--psi6-every={dpsi:.17g}", "--gen3-virial-blocks=50", f"--gen3-snapshots={Tsnap:.17g}"] + extra
    env = dict(ENV, HD_KE_TRACE="1", HD_CONTACT_AUDIT="1")
    meta = dict(N=N, eta=eta, alpha=alpha, H=H, L0=L0, M=M, Teq=Teq, nu_pred=nup, Trec=Trec, dpsi=dpsi, Tsnap=Tsnap, nu_kr=nukr, stride=stride)
    return name, d, argv, env, meta


def af_case(name, L0, xw, extra, hold=600, post=60, N=100):
    H = 10.0 * math.sqrt(N / 100.0); Ns = N // 2
    d = os.path.join(OUT, name); seed = 7700100 + int(name[2])
    argv = [BIN, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0", "--seed-drift-order=drift-first",
            f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv", "--trace-every=60", f"--particles={N}",
            f"--particles-boxes={Ns},{Ns}", "--particle-radius=0.5", f"--l0={L0:.17g}", f"--height={H:g}", "--num-walls=1",
            f"--wall-positions={xw:.17g}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
            "--eff-output=wall-ke", f"--wall-hold-steps={hold}", f"--steps={post}", "--fixed-dt=0.4", "--kbt1", f"--seed={seed}",
            "--engine=gen3"] + extra
    env = dict(ENV, HD_PISTON_EVENTS=os.path.join(d, "ev.csv"))
    return name, d, argv, env, dict(L0=L0, xw=xw, hold=hold, post=post)


def run(case):
    name, d, argv, env, meta = case
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "command.txt"), "w") as fh:
        fh.write(" ".join(f"{k}={v}" for k, v in env.items() if k.startswith("HD_")) + " " + " ".join(argv) + "\n")
    with open(os.path.join(d, "meta.txt"), "w") as fh:
        fh.write(repr(meta) + "\n")
    t = time.time()
    with open(os.path.join(d, "run.log"), "w") as fh:
        p = subprocess.run(["perl", "-e", "alarm 600; exec @ARGV"] + argv, cwd=d, env=env, stdout=fh, stderr=subprocess.STDOUT)
    with open(os.path.join(d, "exit_code.txt"), "w") as fh:
        fh.write(f"{p.returncode}\n")
    return name, p.returncode, time.time() - t


def main():
    if hashlib.sha256(open(BIN, "rb").read()).hexdigest() != BIN_SHA: sys.exit("STOP: the f72d764 Mac binary is not the recorded one")
    L0_30 = 100 * math.pi / (8 * 10 * 0.30)
    cases = [b_case("B1", 100, 0.30, 0.5, 7700001, []),
             b_case("S0", 100, 0.60, 1.0, 7700002, ["--gen3-lattice-orientation=0"]),
             b_case("ST", 100, 0.60, 1.0, 7700003, ["--gen3-lattice-orientation=turned"]),
             af_case("AF1", L0_30, L0_30, ["--gen3-seeding=lattice"]),
             af_case("AF2", L0_30, L0_30, ["--gen3-exact-box"]),
             af_case("AF3", L0_30, L0_30, []),
             af_case("AF4", 628.0 / 48.0, 628.0 / 48.0 + 0.0123, []),
             af_case("AF5", 273.0 / 48.0, 273.0 / 48.0, []),
             af_case("AF6", 1092.0 / 48.0, 1092.0 / 48.0, [], hold=3000, post=60, N=1600),
             af_case("AF5v", 273.0 / 48.0, 273.0 / 48.0, [], hold=60, post=60),
             af_case("AF6v", 1092.0 / 48.0, 1092.0 / 48.0, [], hold=60, post=60, N=1600)]
    for c in cases:
        if c[0].endswith("v"): c[2].remove("--quiet")
    done = [c for c in cases if os.path.exists(c[1])]                  # a case run before is skipped, never overwritten
    for c in done: print(f"{c[0]}: run before, skipped")
    cases = [c for c in cases if c not in done]
    with ThreadPoolExecutor(max_workers=3) as ex:
        for name, rc, s in ex.map(run, cases):
            print(f"{name}: exit {rc}, {s:.1f} s", flush=True)


if __name__ == "__main__":
    main()
