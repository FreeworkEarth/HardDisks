#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.13, third plan-author decision of 2026-10-07, item 3): the ASan/UBSan runs of gate v3.
A SCRATCH build of 7b08827 (asan.sbatch: gcc -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer, written outside the
clone, version line "00ALLINONE  git 7b08827  target asan-scratch") runs four trajectories, each with --resched-audit (mode 1):
  smoke_min, smoke_leg   the smoke trajectory (M = 50, 25 periods; audit_runs_261007.sos_cmd, as its ctrl_ cases), minimal and
                         --legacy-resched
  tT300_min, tT300_leg   the first M = 300 trajectory of Test T (tasks_T_epi8_H_H10_L10.txt, r = 0, both policy lines), with the
                         B-mode command of cluster/confinement_20261013/conf_worker.sh, built from that line
HD_KE_TRACE=1 and HD_CONTACT_AUDIT=1 as in Test T. ASAN_OPTIONS and UBSAN_OPTIONS below (LeakSanitizer ON, the default on Linux).
REQUIREMENT (plan author, item 3): zero sanitizer reports, exit 0, audit missing = extra = 0, in every run. Any report = stop and
paste. A "sanitizer report" here is ANY output line containing "Sanitizer" or "runtime error" (ASan, UBSan, LSan, their SUMMARY
lines and their warnings) -- counted, never filtered.
usage (from hspist3/):
  python3 cluster/resched_gate_261005/asan_runs_261007.py run --bin <scratch asan binary> --out <new dir> [--jobs 4]
  python3 cluster/resched_gate_261005/asan_runs_261007.py report --out <dir>
"""
import argparse, math, os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import audit_runs_261007 as AR

CASES = ["smoke_min", "smoke_leg", "tT300_min", "tT300_leg"]
TASKS_T = os.path.join(HERE, "tasks_T_epi8_H_H10_L10.txt")
ASAN_OPTIONS = "detect_leaks=1:halt_on_error=1:abort_on_error=0:print_summary=1"
UBSAN_OPTIONS = "print_stacktrace=1:print_summary=1"
SAN = re.compile(r"Sanitizer|runtime error")


def tT_line(pol):
    """The first M = 300 line of Test T for this policy (r = 0)."""
    for l in open(TASKS_T):
        f = l.split()
        if f[0] == "B" and f[2] == "300" and f[3] == "0" and f[10] == pol: return f
    sys.exit(f"STOP: no M = 300, r = 0, {pol} line in {TASKS_T}")


def worker_B_cmd(binp, f, d):
    """conf_worker.sh, mode B, verbatim (B <rel> <M> <r> <seed> <L0> <H> <Ns> <stride> <base> <policy>), run directory d."""
    _, rel, M, r, seed, L0, H, NS, stride, base, pol = f
    c = [binp, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first",
         "--edmd-acc=0", f"--particles={2 * int(NS)}", f"--particles-boxes={NS},{NS}", f"--height={H}", "--particle-radius=0.5",
         "--wall-thickness=0.05", "--wall-thickness-vis=0.05", f"--lengths={L0}", f"--wall-masses={M}", "--repeats=1", f"--seed={base}",
         "--wall-hold-steps=2000", "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0",
         "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={stride}",
         f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}"]
    return c + (["--legacy-resched"] if pol == "legacy" else [])


def case_cmd(case, binp, d):
    if case.startswith("smoke_"):
        c = AR.sos_cmd(binp, 50, d, 10.0, 10.0, 50, 25, math.pi / 8) + (["--legacy-resched"] if case == "smoke_leg" else [])
    else:
        c = worker_B_cmd(binp, tT_line("minimal" if case == "tT300_min" else "legacy"), d)
    return c + ["--resched-audit"]


def run(a):
    a.out = os.path.abspath(a.out); binp = os.path.abspath(a.bin)
    if os.path.exists(a.out): sys.exit(f"STOP: {a.out} exists -- not overwriting")
    ver = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout
    if "target asan" not in ver.split("\n")[0] and not os.environ.get("HD_ASAN_RUNNER_TEST"):
        sys.exit(f"STOP: {binp} is not a sanitizer scratch build: {ver.splitlines()[:1]}")
    env = dict(os.environ, HD_KE_TRACE="1", HD_CONTACT_AUDIT="1", ASAN_OPTIONS=ASAN_OPTIONS, UBSAN_OPTIONS=UBSAN_OPTIONS)
    procs = []
    for case in CASES:
        d = os.path.join(a.out, case); os.makedirs(d)
        c = case_cmd(case, binp, d)
        open(os.path.join(d, "version.txt"), "w").write(ver)
        open(os.path.join(d, "command.txt"), "w").write(" ".join(c) + "\n")
        open(os.path.join(d, "env.txt"), "w").write(f"ASAN_OPTIONS={ASAN_OPTIONS}\nUBSAN_OPTIONS={UBSAN_OPTIONS}\nHD_KE_TRACE=1\nHD_CONTACT_AUDIT=1\n")
        while len([p for p in procs if p[1].poll() is None]) >= a.jobs: time.sleep(0.5)
        procs.append((d, subprocess.Popen(c, stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT, env=env, cwd=d), time.time()))
    for d, p, t0 in procs:
        rc = p.wait(); dt = time.time() - t0
        open(os.path.join(d, "rc.txt"), "w").write(f"{rc} {dt:.0f}\n"); print(f"{os.path.basename(d)}: exit {rc} ({dt:.0f} s)")


def report(a):
    a.out = os.path.abspath(a.out)
    print(f"## ASan/UBSan runs of gate v3 ({a.out})\n")
    vers = {open(os.path.join(a.out, c, "version.txt")).read().split("\n")[0] for c in CASES if os.path.exists(os.path.join(a.out, c, "version.txt"))}
    print(f"binary: {sorted(vers)}\n{open(os.path.join(a.out, CASES[0], 'env.txt')).read() if os.path.exists(os.path.join(a.out, CASES[0], 'env.txt')) else ''}")
    print("| run | exit | wall [s] | sanitizer report lines | audit mode | audited events | missing | extra | max contact [px] | ok |\n"
          "|---|---|---|---|---|---|---|---|---|---|")
    ok_all = len(vers) == 1 and all("target asan" in v for v in vers)
    for case in CASES:
        d = os.path.join(a.out, case)
        if not os.path.exists(os.path.join(d, "rc.txt")):
            print(f"| {case} | **not run / not finished** | | | | | | | | **NO** |"); ok_all = False; continue
        rc, wall = open(os.path.join(d, "rc.txt")).read().split()
        log = open(os.path.join(d, "run.log"), errors="ignore").read()
        san = [l for l in log.split("\n") if SAN.search(l)]
        m = AR.AUD.search(log); c = AR.CON.search(log)
        cmax = max(float(v) for v in c.groups()[1:]) if c else float("nan")
        ok = rc == "0" and not san and m is not None and int(m.group(4)) == 0 and int(m.group(5)) == 0
        ok_all &= ok
        print(f"| {case} | {rc} | {wall} | {len(san)} | {m.group(1) if m else '**none**'} | {m.group(2) if m else ''} | "
              f"{m.group(4) if m else ''} | {m.group(5) if m else ''} | {cmax:.1e} | {'yes' if ok else '**NO**'} |")
        for l in san[:20]: print(f"    {case}: {l}")
    print(f"\nASAN (decision 3, item 3): {'CLEAN -- zero sanitizer reports, exit 0, audit missing = extra = 0 in all four runs' if ok_all else 'NOT CLEAN -- stop and paste this report'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "report"])
    ap.add_argument("--bin"); ap.add_argument("--out", required=True); ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args(); sys.exit({"run": run, "report": report}[a.what](a))
