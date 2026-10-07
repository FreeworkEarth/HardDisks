#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10, plan-author decision parts D1, E0 and E2): deterministic tests of the minimal
rescheduling -- the schedule-equivalence audit (--resched-audit), the contact audit (HD_CONTACT_AUDIT=1), the proof that the
audit does not steer (outputs with and without it byte-identical), and E0 (outputs of the new binary with the audit OFF
byte-identical to those of the reference binary, build 73fc07f). Diagnostics only: nothing from these runs enters a figure.

Cases (all N = 100, fixed seeds):
  free_M50, free_M500, free_M1500, free_M2000  epi8_H_H10_L10 (eta = pi/8, H = L_0 = 10, N_s = 50), the speed-of-sound protocol
      of the campaign (divider held 2000 steps, then free), 200 periods, the smoke pilot's seed of that mass
      (run_seed(20261013, 0, m, 0)); free divider: every divider hit changes the divider velocity (epoch path)
  afix          the A-fixed cell epi8_H_H10_L10 at x_0, seed 9700, the campaign protocol (divider held 312000 steps, then
      1200 released); held divider: no epoch change, only the hitting disk is rescheduled
  dense_M50, dense_M2000   eta ~ 0.70 (H = 10, L_0 = 269/48 = 5.604167, N_s = 50), the speed-of-sound protocol, 25 periods
  ctrl_min, ctrl_leg   the smoke trajectory (M = 50, 25 periods) with the audit after EVERY event (mode 2), minimal and
      legacy: the legacy row shows the |dt| that rounding alone gives for events kept in the heap
  free_M50_long  (KOA, E2) as free_M50 over 600 periods, so the lightest mass has >= 1e5 divider events
  afix_leg       (KOA, E0) the A-fixed trajectory on the legacy path
Variants per case: "audit" (--resched-audit or =all, HD_CONTACT_AUDIT=1), "plain" (no audit flag, HD_CONTACT_AUDIT=1), and
with --ref-bin "ref" (the reference binary, HD_CONTACT_AUDIT=1). Compared files: speed-of-sound wall trace and psi6 file;
energy-transfer ev_, tr_ and red_ files.
RULE (plan author, written before running): any missing or extra event, or |dt| > 1e-9, in a mode-1 audit of the minimal
path = a defect in the schedule logic; zero over all runs = the minimal schedule equals the legacy schedule on these
trajectories. The mode-2 controls are information.
usage (from hspist3/):
  python3 cluster/resched_gate_261005/audit_runs_261007.py run --bin <new> [--ref-bin <ref>] --out <dir> [--jobs 8] [--cases a,b]
  python3 cluster/resched_gate_261005/audit_runs_261007.py report --out <dir>
"""
import argparse, glob, math, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
sys.path.insert(0, CL); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T

AUD = re.compile(r"\[EDMD-AUDIT\] mode (\d): audited events (\d+); matched comparisons (\d+); missing (\d+); extra (\d+); "
                 r"\|dt\| > 1e-9 (\d+); duplicate live events disagreeing (\d+); max \|dt\| over matched (\S+); "
                 r"\|dt\| > 1e-9 and > 1e-10 of the horizon (\d+); max \|dt\|/horizon (\S+)")
CON = re.compile(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), outer walls (\S+), "
                 r"divider (\S+), pistons (\S+)\n")
CASES = ["free_M50", "free_M500", "free_M1500", "free_M2000", "afix", "dense_M50", "dense_M2000", "ctrl_min", "ctrl_leg",
         "free_M50_long", "afix_leg"]   # the last two for KOA (E0/E2): M = 50 over 600 periods (>= 1e5 divider events), A-fixed legacy


def sos_cmd(binp, M, out, L0, H, NS, target, eta):
    """The speed-of-sound command of cluster/confinement_pilot.py cmd(), with the cell as parameters."""
    mi = T.A1_MASSES.index(M); x = T.x_of(M, L0); nu = T.kr_cs(eta) * x; seed = T.run_seed(20261013, 0, mi, 0)
    return [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first",
            "--edmd-acc=0", f"--particles={2 * NS}", f"--particles-boxes={NS},{NS}", f"--height={H}", "--particle-radius=0.5",
            "--wall-thickness=0.05", "--wall-thickness-vis=0.05", f"--lengths={L0:.4f}", f"--wall-masses={M}", "--repeats=1",
            "--seed=20261013", "--wall-hold-steps=2000", "--fixed-dt=0.4", f"--target-oscillations={target}",
            "--oscillation-safety=1.0", "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000",
            f"--speed-sound-log-stride={T.d_stride(nu)}", f"--speed-sound-run-dir={out}", f"--speed-sound-exact-seed={seed}"]


def afix_cmd(binp, d):
    return [binp, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0",
            "--seed-drift-order=drift-first", f"--energy-transfer-summary={d}/summary_9700.csv",
            f"--energy-transfer-trace={d}/tr_9700.csv", "--trace-every=600", "--particles=100", "--particles-boxes=50,50",
            "--particle-radius=0.5", "--l0=10.000000", "--height=10.000000", "--num-walls=1", "--wall-positions=10.000000",
            "--wall-mass-factors=1000000000", "--wall-thickness=0.05", "--wall-thickness-vis=0.05", "--eff-output=wall-ke",
            "--wall-hold-steps=312000", "--steps=1200", "--fixed-dt=0.4", "--kbt1", "--seed=9700"]


def case_cmd(case, binp, d, variant):
    eta8 = math.pi / 8; L0d = 269 / 48
    if case == "free_M50_long": c = sos_cmd(binp, 50, d, 10.0, 10.0, 50, 600, eta8)
    elif case.startswith("free_M"): c = sos_cmd(binp, int(case[6:]), d, 10.0, 10.0, 50, 200, eta8)
    elif case.startswith("dense_M"): c = sos_cmd(binp, int(case[7:]), d, L0d, 10.0, 50, 25, 50 * math.pi * 0.25 / (10.0 * L0d))
    elif case in ("afix", "afix_leg"): c = afix_cmd(binp, d)
    else: c = sos_cmd(binp, 50, d, 10.0, 10.0, 50, 25, eta8)
    if case in ("ctrl_leg", "afix_leg"): c = c + ["--legacy-resched"]
    if variant == "audit": c = c + (["--resched-audit=all"] if case.startswith("ctrl_") else ["--resched-audit"])
    return c


def run(a):
    a.out = os.path.abspath(a.out); os.makedirs(a.out, exist_ok=True)
    cases = a.cases.split(",") if a.cases else CASES
    jobs = []
    for case in cases:
        for variant, binp in [("audit", a.bin), ("plain", a.bin)] + ([("ref", a.ref_bin)] if a.ref_bin else []):
            d = os.path.join(a.out, case, variant)
            if os.path.exists(d): sys.exit(f"{d} exists -- not overwriting")
            os.makedirs(d); binp = os.path.abspath(binp)
            c = case_cmd(case, binp, d, variant)
            open(os.path.join(d, "version.txt"), "w").write(subprocess.run([binp, "--version"], capture_output=True, text=True).stdout)
            env = dict(os.environ, HD_KE_TRACE="1", HD_CONTACT_AUDIT="1")
            if case.startswith("afix"): env["HD_PISTON_EVENTS"] = os.path.join(d, "ev_9700.csv")
            open(os.path.join(d, "command.txt"), "w").write(" ".join(c) + "\n")
            jobs.append((d, c, env))
    procs = []; t0 = time.time()
    for d, c, env in jobs:
        while len([p for p in procs if p[1].poll() is None]) >= a.jobs: time.sleep(0.5)
        procs.append((d, subprocess.Popen(c, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT, env=env, cwd=d), time.time()))
    for d, p, ts in procs:
        rc = p.wait(); print(f"{os.path.relpath(d, a.out)}: exit {rc}")
    t1 = 312000 * 0.4 / 24.0
    for case in cases:
        if not case.startswith("afix"): continue
        for d in glob.glob(os.path.join(a.out, case, "*")):
            subprocess.run([sys.executable, os.path.join(CL, "confinement_20261013", "reduce_AF.py"), f"{d}/ev_9700.csv", f"{d}/tr_9700.csv",
                            f"{d}/red_9700.csv", "200", f"{t1:.9f}"])
    print(f"wall time {time.time() - t0:.0f} s")


def outputs(d):
    fs = sorted(glob.glob(os.path.join(d, "wall_x_positions_*_run0.csv"))) + sorted(glob.glob(os.path.join(d, "speed_of_sound_psi6.csv")))
    return fs + [os.path.join(d, f) for f in ("ev_9700.csv", "tr_9700.csv", "red_9700.csv") if os.path.exists(os.path.join(d, f))]


def same(d1, d2):
    f1, f2 = outputs(d1), outputs(d2)
    if [os.path.basename(f) for f in f1] != [os.path.basename(f) for f in f2] or not f1: return "files differ"
    return "IDENTICAL" if all(subprocess.run(["cmp", "-s", x, y]).returncode == 0 for x, y in zip(f1, f2)) else "**DIFFERENT**"


def report(a):
    from resched_audit_bruteforce_261007 import report as bf
    a.out = os.path.abspath(a.out)
    print(f"## Schedule-equivalence audit and contact audit ({a.out})\n")
    print("| case | version (audit run) | mode | audited events | matched | missing | extra | abs(dt) > 1e-9 | duplicate live disagreeing | "
          "max abs(dt) matched | abs(dt) > 1e-9 and > 1e-10 of horizon | max abs(dt)/horizon | max contact gap [px] (dd, wall, div, piston) | "
          "audit vs plain | plain vs ref |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    viol = 0; viol_rel = 0; evlines = []
    for case in [c for c in CASES if os.path.isdir(os.path.join(a.out, c))]:
        da, dp, dr = (os.path.join(a.out, case, v) for v in ("audit", "plain", "ref"))
        log = open(os.path.join(da, "run.log"), errors="ignore").read(); m = AUD.search(log); c = CON.search(log)
        ver = open(os.path.join(da, "version.txt")).read().split("\n")[0].replace("00ALLINONE  ", "")
        if not m:
            print(f"| {case} | {ver} | **no [EDMD-AUDIT] line** | | | | | | | | | | |"); viol += 1; continue
        mode, n, cmp_, miss, ext, dtn, dup, mx, dtr, mr = m.groups()
        if mode == "1" and case not in ("ctrl_leg", "afix_leg"): viol += int(miss) + int(ext) + int(dtn); viol_rel += int(miss) + int(ext) + int(dtr)
        evlines += [f"{case}: {l}" for l in log.split("\n") if l.startswith("[EDMD-AUDIT-EV]")]
        cg = ", ".join(f"{float(v):.1e}" for v in c.groups()[1:]) if c else "**none**"
        print(f"| {case} | {ver} | {mode} | {n} | {cmp_} | {miss} | {ext} | {dtn} | {dup} | {float(mx):.2e} | {dtr} | {float(mr):.2e} | {cg} | {same(da, dp)} | "
              f"{same(dp, dr) if os.path.isdir(dr) else 'n/a'} |")
    print(f"\nRULE (mode-1 audits of the minimal path): missing + extra + (abs(dt) > 1e-9) summed over all runs = {viol} -> "
          f"{'ZERO: the minimal schedule equals the legacy schedule on these trajectories' if viol == 0 else 'NON-ZERO: defect in the schedule logic by the rule'}")
    print(f"same with the relative criterion (abs(dt) > 1e-9 AND > 1e-10 of the prediction horizon; information, not registered): {viol_rel}")
    print("\n### Every reported event, recomputed at 60 digits (validation/resched_audit_bruteforce_261007.py)\n")
    if evlines:
        print("(the first 50 reported events of each run; the engine counts all of them, see the table above)\n")
        for case in sorted({l.split(": ", 1)[0] for l in evlines}, key=CASES.index):
            print(f"**{case}**"); bf([l.split(": ", 1)[1] for l in evlines if l.startswith(case + ": ")], table=False)
        print("\nfirst 12 events in full:\n"); bf([l.split(": ", 1)[1] for l in evlines[:12]])
    else:
        print("(no [EDMD-AUDIT-EV] line in any run)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "report"])
    ap.add_argument("--bin"); ap.add_argument("--ref-bin"); ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--cases"); a = ap.parse_args()
    sys.exit({"run": run, "report": report}[a.what](a))
