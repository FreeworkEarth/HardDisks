#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.19; stage F of the plan-author programme of sec. 4.7.12): PILOT P2, drift with an equilibrated
start. A pilot, no paper use. Question: does the +7.5 / +10.4 / +7.8 % first-half to second-half frequency drift of A1 v2 (sec. 4.7.11)
survive a settled start?
gen3, M4 seeding (--gen3-seeding=lattice, jitter 0.25), N = 100, the A1 v2 geometry: H = 10 and the A1 v2 box exactly -- the binary's
1/48-sigma grid box of each A1 v2 cell (SIM_WIDTH = (int)(2 L0 x 24) px with A1 v2's L0 = 5.6913, 5.5702, 5.531, 5.4923), passed as
an exact length (--gen3-exact-box), so eta_true = 0.6905, 0.7060, 0.7113, 0.7167 as in A1 v2. M = 50, 300, 2000; 12 seeds per cell,
run_seed(20261218, 0, cell index, r) (base never used before). Hold T_eq = min(max(1e4, 10 tau_psi6), 5e4) sigma-time with tau_psi6
the mean integrated autocorrelation time of P1 at N = 100 and the nearest eta of P1's grid (0.6905 -> 0.70, 0.7060 -> 0.704,
0.7113 -> 0.712, 0.7167 -> 0.716; from p1_tables_output.txt, read by --teq); then release and record 200 predicted periods (the
speed-of-sound protocol: --target-oscillations=200, the campaign stride T.d_stride of the KR prediction). --psi6-every=1 (as P1;
the default 0.25 would cost 4 x the disk over holds of up to 5e4 sigma-time). Each run in its own folder, at most 12 processes,
the psi6(t) file compressed after its run (SHA-256 of the uncompressed file first). Rule 10: a dirty or missing run record is
reported and kept, never rerun.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p2_261009/p2_run.py --bin <frozen 00ALLINONE> --expect-bin-sha <sha256> --teq <p1_tables_output.txt> --out <root>
"""
import argparse, hashlib, math, os, re, subprocess, sys, time

MAIN_HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
sys.path.insert(0, os.path.join(MAIN_HS, "validation"))
import tests_20260913 as T
import numpy as np
A1V2 = ((0.6905, 5.6913, 0.70), (0.7060, 5.5702, 0.704), (0.7113, 5.531, 0.712), (0.7167, 5.4923, 0.716))   # eta_true, A1 v2 L0, P1 eta
MASSES, NSEED, BASE, H = (50, 300, 2000), 12, 20261218, 10.0


def grid_L0(L0):
    """The A1 v2 box: the binary's SIM_WIDTH = (int)(2 * L0_UNITS * 24) with L0_UNITS a float; the exact length is SIM_WIDTH / 48."""
    return math.floor(float(np.float32(2.0) * np.float32(L0) * np.float32(24.0))) / 48.0


def teq_table(path):
    """T_eq per P1 eta at N = 100 from p1_tables_output.txt: min(max(1e4, 10 x mean tau), 5e4)."""
    out = {}
    for l in open(path):
        if l.startswith("## psi6 LOCAL"): break             # the global table (the primary measure) only
        m = re.match(r"\| 100 \| (\d\.\d+) \| .*?\| ([\d.e+-]+|nan); [^|]*\| (yes|\*\*TAU NOT RESOLVED\*\*) \|", l)
        if m:
            tau = float(m.group(2))
            out[round(float(m.group(1)), 3)] = (min(max(1e4, 10 * tau), 5e4) if math.isfinite(tau) else 5e4, tau, m.group(3))
    return out


def cells(teq):
    out = []
    for k, (eta_t, L0a, eta_p1) in enumerate(A1V2):
        T_eq, tau, res = teq[round(eta_p1, 3)]
        L0 = grid_L0(L0a)
        for M in MASSES:
            nu = T.kr_cs(eta_t) * T.x_of(M, L0)
            out.append(dict(k=len(out), eta_t=eta_t, L0=L0, L0a=L0a, M=M, T_eq=T_eq, tau=tau, tau_res=res, eta_p1=eta_p1,
                            stride=T.d_stride(nu), seeds=[T.run_seed(BASE, 0, len(out), r) for r in range(NSEED)],
                            name=f"eta{eta_t:.4f}_M{M}"))
    return out


def cmd(binp, d, c, seed):
    hold = int(round(c["T_eq"] * 60.0))
    return [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
            "--particles=100", "--particles-boxes=50,50", f"--height={H:g}", "--particle-radius=0.5", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", f"--lengths={c['L0']:.17g}", f"--wall-masses={c['M']}", "--repeats=1", f"--seed={BASE}",
            f"--wall-hold-steps={hold}", "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={c['stride']}",
            f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}", "--engine=gen3", "--gen3-seeding=lattice", "--gen3-exact-box",
            "--psi6-every=1"]


def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--expect-bin-sha", required=True)
    ap.add_argument("--teq", required=True); ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=12)
    a = ap.parse_args()
    if a.jobs > 12: sys.exit("STOP: at most 12 processes (programme rule 1)")
    binp, O = os.path.abspath(a.bin), os.path.abspath(a.out); os.makedirs(O, exist_ok=True)
    if sha(binp) != a.expect_bin_sha: sys.exit("STOP: binary SHA-256")
    teq = teq_table(a.teq)
    C = cells(teq)
    log = open(os.path.join(O, "runner.log"), "a")
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    build = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S %Z')} P2: binary {binp} ({a.expect_bin_sha}; {build}); T_eq from {a.teq}")
    for c in C: say(f"  {c['name']}: L0 {c['L0']:.17g} (A1 v2 {c['L0a']}), T_eq {c['T_eq']:g} sigma-time (P1 eta {c['eta_p1']}: tau {c['tau']:.4g}, "
                    f"{c['tau_res']}), stride {c['stride']}")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())
    jobs = [(c, r, s, os.path.join(O, c["name"], f"seed{r}")) for c in sorted(C, key=lambda c: -c["T_eq"] * (1 + c["M"] / 1000)) for r, s in enumerate(c["seeds"])]
    run, t0 = [], time.time()
    def reap():
        for j in list(run):
            p, c, r, d, ts = j
            if p.poll() is None: continue
            run.remove(j)
            lg = open(os.path.join(d, "run.log"), errors="ignore").read()
            recs = re.findall(r"^\[EDMD3-HEALTH\] .*?: clean=(\S+) ", lg, re.M); nb = len(re.findall(r"^\[EDMD3\] built #", lg, re.M))
            bad = "" if recs == ["1"] and nb == 1 and p.returncode == 0 else f"exit {p.returncode}, run record {recs}, builds {nb}"
            for f in [f for f in os.listdir(d) if f.startswith("psi6_t_") and f.endswith(".csv")]:
                with open(os.path.join(os.path.dirname(d), ".sha256_uncompressed"), "a") as fh: fh.write(f"{sha(os.path.join(d, f))}  seed{r}/{f}\n")
                subprocess.run(["gzip", "-9", os.path.join(d, f)], check=True)
            say(f"  {c['name']} seed{r}: exit {p.returncode}, {time.time() - ts:.0f} s" + (f"  **FINDING: {bad}**" if bad else ""))
    for c, r, s, d in jobs:
        while len(run) >= a.jobs: reap(); time.sleep(0.5)
        os.makedirs(d, exist_ok=True); cm = cmd(binp, d, c, s)
        with open(os.path.join(d, "command.txt"), "w") as fh: fh.write(" ".join(cm) + "\n")
        run.append((subprocess.Popen(cm, cwd=d, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT), c, r, d, time.time()))
    while run: reap(); time.sleep(0.5)
    say(f"# done {time.strftime('%Y-%m-%d %H:%M:%S %Z')} after {time.time() - t0:.0f} s")
    say(subprocess.run(["df", "-h", O], capture_output=True, text=True).stdout.strip())


if __name__ == "__main__":
    main()
