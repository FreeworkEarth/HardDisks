#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.16; stage C of the programme of sec. 4.7.12): the M4 evidence runs, every one with the frozen
binaries of the committed tree f39e485 (build_clean.sh; M4 = 3c072fb + the second lattice form f39e485), each run in its own folder, at most --jobs processes at once (default 12;
programme rule 1: at most 12 simulation processes in all, so it starts only after Test G has finished).
  e0        rule 4: gen2 byte identity after the M4 change of 00ALLINONE.c (ctrl_min, ctrl_leg; default and --engine=gen2)
  harness   rule 8 / the engine change (edmd3_tie_stats, read-only): gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit
  grid      the acceptance grid: N = 100, 400, 900, 1600 x eta = 0.10, pi/8, 0.60, 0.70, 0.716, 0.78, 0.85, 0.90; the gate's dense
            construction H = 10 sqrt(N/100), L0 = N pi / (8 H eta) exactly (--gen3-exact-box, %.17g), M4 seeding (--gen3-seeding=
            lattice, the default jitter 0.25), divider held. Per cell three runs: "long" (seed s1: 400 sigma-time held, 1 sigma-time
            released, the initial state dumped), "same" (s1 again, 1 sigma-time held: the initial state must be bit-identical), "other"
            (s2: it must differ). Seeds run_seed(20261216, 0, cell index, 0 and 1); base 20261216 never used before.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m4_261009/stageC_evidence.py --bin-dir <frozen binaries> --out <evidence dir> [--only e0,harness,grid] [--jobs 4]
"""
import argparse, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
MAIN_HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
GATE = os.path.join(MAIN_HS, "cluster", "resched_gate_261005")
sys.path.insert(0, os.path.join(MAIN_HS, "validation"))
import tests_20260913 as T
NS_LIST = (100, 400, 900, 1600)
ETAS = (("0.10", 0.10), ("pi/8", math.pi / 8), ("0.60", 0.60), ("0.70", 0.70), ("0.716", 0.716), ("0.78", 0.78), ("0.85", 0.85), ("0.90", 0.90))
BASE = 20261216


def cells():
    out = []
    for N in NS_LIST:
        for lab, eta in ETAS:
            k = len(out); H = 10.0 * math.sqrt(N / 100.0); L0 = N * math.pi / (8.0 * H * eta)
            out.append(dict(k=k, N=N, lab=lab, eta=eta, H=H, L0=L0, s1=T.run_seed(BASE, 0, k, 0), s2=T.run_seed(BASE, 0, k, 1)))
    return out


def cmd(binp, d, c, seed, hold):
    NS = c["N"] // 2
    return [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={c['N']}", f"--particles-boxes={NS},{NS}", f"--height={c['H']:.17g}", "--particle-radius=0.5", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", f"--lengths={c['L0']:.17g}", "--wall-masses=300", "--repeats=1", f"--seed={BASE}",
            f"--wall-hold-steps={hold}", "--fixed-dt=0.4", "--record-sigma-time=1", "--oscillation-min-steps=10", "--speed-sound-log-stride=60",
            f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}", "--engine=gen3", "--gen3-seeding=lattice", "--gen3-exact-box",
            f"--gen3-dump-initial={d}/init.txt"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--only"); ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    only = set(a.only.split(",")) if a.only else None
    want = lambda k: only is None or k in only
    BIN = os.path.join(B, "00ALLINONE")
    print(f"binaries: {B}; out: {O}; {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", flush=True)
    print(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout, flush=True)
    jobs = []
    if want("e0"):
        w = os.path.join(B, "bin_engine_gen2.sh")
        with open(w, "w") as fh: fh.write(f"#!/bin/bash\nexec {BIN} \"$@\" --engine=gen2\n")
        os.chmod(w, 0o755)
        for tag, b in (("default", BIN), ("gen2flag", w)):
            jobs.append((f"e0 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "e0", tag),
                                       "--cases", "ctrl_min,ctrl_leg", "--jobs", "2"], os.path.join(O, "e0"), os.path.join(O, "e0", f"run_{tag}.txt"), {}, 2, False))
    if want("harness"):
        h = os.path.join(O, "harness")
        jobs.append(("harness", f"'{B}/gen3_m1' audit > m1_audit_output.txt && '{B}/gen3_m2' audit --quick > m2_audit_quick_output.txt && "
                                f"'{B}/gen3_m2' audit > m2_audit_output.txt", h, os.path.join(h, "run.txt"), {}, 1, True))
    if want("grid"):
        for c in cells():
            for tag, seed, hold in (("long", c["s1"], 24000), ("same", c["s1"], 60), ("other", c["s2"], 60)):
                d = os.path.join(O, "grid", f"N{c['N']}_eta{c['lab'].replace('/', '')}", tag)
                jobs.append((f"grid N{c['N']} {c['lab']} {tag}", cmd(BIN, d, c, seed, hold), d, os.path.join(d, "run.log"), {"HD_CONTACT_AUDIT": "1"}, 1, False))
    run, done, t0 = [], [], time.time()
    def reap():
        for j in list(run):
            p, name, ts, wgt = j
            if p.poll() is not None:
                run.remove(j); done.append((name, p.returncode)); print(f"  done {name}: rc {p.returncode}, {time.time() - ts:.0f} s", flush=True)
    for name, c, d, log, env, wgt, shell in jobs:
        while sum(x[3] for x in run) + wgt > a.jobs: reap(); time.sleep(0.3)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write((c if shell else " ".join(c)) + "\n")
        p = subprocess.Popen(c, cwd=d, stdout=open(log, "w"), stderr=subprocess.STDOUT, env=dict(os.environ, **env), shell=shell)
        run.append((p, name, time.time(), wgt))
    while run: reap(); time.sleep(0.3)
    print(f"finished {time.strftime('%H:%M:%S')} after {time.time() - t0:.0f} s; non-zero exits: " +
          (", ".join(f"{n} {rc}" for n, rc in done if rc != 0) or "none"), flush=True)


if __name__ == "__main__":
    main()
