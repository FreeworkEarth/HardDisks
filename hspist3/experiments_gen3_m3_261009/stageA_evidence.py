#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.14; stage A of the programme of sec. 4.7.12): the M3 evidence runs, every one with the
frozen binaries of a committed tree (build_clean.sh; programme rule 5), each run in its own folder, at most 12 processes at once.
  e0        acceptance 1 / rule 4: gen2 byte identity, ctrl_min and ctrl_leg through cluster/resched_gate_261005/audit_runs_261007.py
            (default build, and --engine=gen2 through a wrapper), the 7b08827 references copied next to them for its report
  harness   rule 8 after the engine fix: gen3_m1 audit, gen3_m2 audit --quick, gen3_m2 audit (compared byte for byte with the
            committed outputs); gen3_body_rule_test 20000 (amendment b); gen3_band_edge_test and gen3_band_edge_ties (amendment a)
  replay    acceptances 2 and 3: run_replays.sh (16 harness cells x 4 reader / stop cadences through the driver)
  prod      amendment d: N = 400, eta 0.70, free divider M = 500, 2e4 sigma-time record, contact audit, schedule audit every 1e6 events
  a3        stage A3: the energy-transfer loop (the A-fixed protocol) under gen3 at the default cadence, --validator-every=1, and
            --validator-every=1 --trace-every=1
  acc3      acceptance 3 (speed of sound): one trajectory at the default cadence, at every cadence maximal, and minimal; acc6:
            the default run repeated (same-seed determinism)
  acc4      acceptance 4: the gate's cases under gen3 (free_M50, free_M500, free_M1500, free_M2000, dense_M50, dense_M2000, ctrl_min,
            afix), plain runs, contact audit
  acc5      acceptance 5, information: gen2 and gen3 on one fluid cell (N = 400 at pi/8, H = L0 = 20, M = 300), 8 seeds each, the
            same seeds (run_seed(20261091, 5, 0, r)), reduce_B.py's argmax estimator
  a4        stage A4, alone after everything else (timings): both loops at N = 400 and 1600, pi/8 and eta 0.70
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_m3_261009/stageA_evidence.py --bin-dir <frozen binaries> --out <evidence dir> [--only e0,harness,...]
"""
import argparse, math, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
MAIN_HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
GATE = os.path.join(MAIN_HS, "cluster", "resched_gate_261005")
sys.path.insert(0, os.path.join(MAIN_HS, "validation")); sys.path.insert(0, MAIN_HS); sys.path.insert(0, GATE)
import audit_runs_261007 as AR          # the gate's commands (sos_cmd, afix_cmd, case_cmd)
import tests_20260913 as T
E0REF = "/private/tmp/claude-501/-Users-chrisharing-Desktop-CCS-complex-coupled-systems-Repo-HardDisks/91cb08ec-0599-4faa-a6af-d5a7ca834255/scratchpad/e0head/out"
MAXP = 12


def geo(N, eta):
    """The gate's dense construction: H = 10 sqrt(N/100), NS = N/2 per side, L0 from eta on the 1/24-sigma grid (floor), as the
    gate's dense cell (269/48 at N = 100) and the production run (269/24 at N = 400); pi/8: L0 = 2 NS / H exactly."""
    H = 10.0 * math.sqrt(N / 100.0); NS = N // 2
    L0 = 2.0 * NS / H if abs(eta - math.pi / 8) < 1e-12 else math.floor(NS * math.pi * 0.25 / (H * eta) * 24.0) / 24.0
    return H, NS, L0


def sos(binp, d, N, eta, M, seed, extra=(), exact_seed=None, target=200):
    H, NS, L0 = geo(N, eta)
    c = [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
         f"--particles={2 * NS}", f"--particles-boxes={NS},{NS}", f"--height={H:g}", "--particle-radius=0.5", "--wall-thickness=0.05",
         "--wall-thickness-vis=0.05", f"--lengths={L0:.6f}", f"--wall-masses={M}", "--repeats=1", f"--seed={seed}", "--wall-hold-steps=2000",
         "--fixed-dt=0.4", f"--target-oscillations={target}", "--oscillation-safety=1.0", "--oscillation-min-steps=10000",
         "--oscillation-max-steps=400000000", "--speed-sound-log-stride=60", f"--speed-sound-run-dir={d}"]
    if exact_seed is not None: c.append(f"--speed-sound-exact-seed={exact_seed}")
    return c + list(extra)


def et(binp, d, N, eta, seed, hold, steps, extra=()):
    H, NS, L0 = geo(N, eta)
    return [binp, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0", "--seed-drift-order=drift-first",
            f"--energy-transfer-summary={d}/summary_{seed}.csv", f"--energy-transfer-trace={d}/tr_{seed}.csv", "--trace-every=600",
            f"--particles={2 * NS}", f"--particles-boxes={NS},{NS}", "--particle-radius=0.5", f"--l0={L0:.6f}", f"--height={H:.6f}",
            "--num-walls=1", f"--wall-positions={L0:.6f}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", "--eff-output=wall-ke", f"--wall-hold-steps={hold}", f"--steps={steps}", "--fixed-dt=0.4",
            "--kbt1", f"--seed={seed}"] + list(extra)


class Pool:
    def __init__(self): self.run = []; self.done = []
    def wait_slots(self, need=1):
        while True:
            self.reap()
            if sum(w for _, _, _, w in self.run) + need <= MAXP: return
            time.sleep(0.5)
    def reap(self):
        for j in list(self.run):
            p, name, t0, w = j
            if p.poll() is not None:
                self.run.remove(j); self.done.append((name, p.returncode, time.time() - t0)); print(f"  done {name}: rc {p.returncode}, {time.time() - t0:.0f} s", flush=True)
    def start(self, name, cmd, cwd, log, env=None, weight=1, shell=False):
        self.wait_slots(weight); os.makedirs(cwd, exist_ok=True)
        with open(os.path.join(cwd, "command.txt"), "w") as fh: fh.write((cmd if shell else " ".join(cmd)) + "\n")
        p = subprocess.Popen(cmd, cwd=cwd, stdout=open(log, "w"), stderr=subprocess.STDOUT, env=dict(os.environ, **(env or {})), shell=shell)
        self.run.append((p, name, time.time(), weight))
    def join(self):
        while self.run: self.reap(); time.sleep(0.5)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin-dir", required=True); ap.add_argument("--out", required=True); ap.add_argument("--only")
    a = ap.parse_args(); B = os.path.abspath(a.bin_dir); O = os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    only = set(a.only.split(",")) if a.only else None
    want = lambda k: only is None or k in only
    BIN = os.path.join(B, "00ALLINONE")
    print(f"binaries: {B}; out: {O}; {time.strftime('%Y-%m-%d %H:%M:%S %Z')}", flush=True)
    print(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout, flush=True)
    P = Pool()
    if want("e0"):
        w = os.path.join(B, "bin_engine_gen2.sh")
        with open(w, "w") as fh: fh.write(f"#!/bin/bash\nexec {BIN} \"$@\" --engine=gen2\n")
        os.chmod(w, 0o755)
        for tag, b in (("default", BIN), ("gen2flag", w)):
            P.start(f"e0 {tag}", [sys.executable, os.path.join(GATE, "audit_runs_261007.py"), "run", "--bin", b, "--out", os.path.join(O, "e0", tag),
                                  "--cases", "ctrl_min,ctrl_leg", "--jobs", "4"], os.path.join(O, "e0"), os.path.join(O, "e0", f"run_{tag}.txt"), weight=4)
    if want("harness"):
        h = os.path.join(O, "harness")
        cmd = (f"'{B}/gen3_m1' audit > m1_audit_output.txt && '{B}/gen3_m2' audit --quick > m2_audit_quick_output.txt && "
               f"'{B}/gen3_m2' audit > m2_audit_output.txt && '{B}/gen3_body_rule_test' 20000 > body_rule_test_output.txt; "
               f"'{B}/gen3_band_edge_test' > band_edge_output.txt; echo band_edge_test rc $?; '{B}/gen3_band_edge_ties' > band_edge_ties_output.txt; "
               f"echo band_edge_ties rc $?")
        P.start("harness", cmd, h, os.path.join(h, "run.txt"), shell=True)
    if want("replay"):
        r = os.path.join(O, "replay")
        P.start("replay", f"bash '{HERE}/run_replays.sh' '{BIN}' '{B}' '{r}/states' > replay_output.txt", r, os.path.join(r, "run.txt"), shell=True)
    if want("prod"):
        d = os.path.join(O, "prod")
        P.start("prod", sos(BIN, d, 400, 0.70, 500, 20261009, ["--engine=gen3", "--record-sigma-time=20000", "--gen3-audit-every=1000000"]),
                d, os.path.join(d, "run.log"), env={"HD_CONTACT_AUDIT": "1"})
    if want("a3"):
        for tag, ex in (("default", []), ("validator1", ["--validator-every=1"]), ("validator1_trace1", ["--validator-every=1", "--trace-every=1"])):
            d = os.path.join(O, "a3", tag)
            P.start(f"a3 {tag}", AR.afix_cmd(BIN, d) + ["--engine=gen3"] + ex, d, os.path.join(d, "run.log"),
                    env={"HD_PISTON_EVENTS": os.path.join(d, "ev_9700.csv")})
    if want("acc3"):
        seed = T.run_seed(20261013, 0, T.A1_MASSES.index(300), 0)
        for tag, ex in (("default", []), ("max_cadence", ["--validator-every=1", "--psi6-every=0.05", "--speed-sound-log-stride=1"]),
                        ("min_cadence", ["--validator-every=1000000", "--psi6-every=1000", "--speed-sound-log-stride=600"]), ("default_repeat", [])):
            d = os.path.join(O, "acc3", tag)
            P.start(f"acc3 {tag}", sos(BIN, d, 100, math.pi / 8, 300, 20261013, ["--engine=gen3"] + ex, exact_seed=seed), d, os.path.join(d, "run.log"))
    if want("acc4"):
        for case in ("free_M50", "free_M500", "free_M1500", "free_M2000", "dense_M50", "dense_M2000", "ctrl_min", "afix"):
            d = os.path.join(O, "acc4", case)
            env = {"HD_CONTACT_AUDIT": "1"}
            if case == "afix": env["HD_PISTON_EVENTS"] = os.path.join(d, "ev_9700.csv")
            os.makedirs(d, exist_ok=True)
            P.start(f"acc4 {case}", AR.case_cmd(case, BIN, d, "plain") + ["--engine=gen3"], d, os.path.join(d, "run.log"), env=env)
    if want("acc5"):
        for r in range(8):
            s = T.run_seed(20261091, 5, 0, r)
            for eng in ("gen2", "gen3"):
                d = os.path.join(O, "acc5", eng, "cell", "m_300", f"r{r}")
                P.start(f"acc5 {eng} r{r}", sos(BIN, d, 400, math.pi / 8, 300, 20261091, [f"--engine={eng}"], exact_seed=s), d, os.path.join(d, "run.log"))
    P.join()
    if want("acc5"):      # reduce_B.py reads <cell>/m_<M>/wall_x_positions_*: the traces are moved up from the run folders
        for eng in ("gen2", "gen3"):
            m = os.path.join(O, "acc5", eng, "cell", "m_300")
            for r in range(8):
                src = os.path.join(m, f"r{r}")
                for f in os.listdir(src):
                    if f.startswith("wall_x_positions_"): shutil.move(os.path.join(src, f), os.path.join(m, f.replace("_run0.csv", f"_run{r}.csv")))
            subprocess.run([sys.executable, os.path.join(MAIN_HS, "cluster", "confinement_20261013", "reduce_B.py"), os.path.join(O, "acc5", eng, "cell")],
                           stdout=open(os.path.join(O, "acc5", f"reduce_{eng}.txt"), "w"), stderr=subprocess.STDOUT)
    if want("a4"):        # alone: one process at a time
        for loop in ("sos", "et"):
            for N in (400, 1600):
                for eta, lab in ((math.pi / 8, "pi8"), (0.70, "070")):
                    d = os.path.join(O, "a4", f"{loop}_N{N}_{lab}")
                    c = (sos(BIN, d, N, eta, 300, 20261092, ["--engine=gen3", "--record-sigma-time=2000"]) if loop == "sos"
                         else et(BIN, d, N, eta, 9711, 60000, 1200, ["--engine=gen3"]))
                    P.start(f"a4 {loop} N{N} {lab}", c, d, os.path.join(d, "run.log")); P.join()
    print(f"finished {time.strftime('%H:%M:%S')}; exit codes: " + ", ".join(f"{n} {rc}" for n, rc, _ in P.done), flush=True)


if __name__ == "__main__":
    main()
