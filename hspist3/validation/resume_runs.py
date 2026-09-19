#!/usr/bin/env python3
"""##CHRIS 2026-09-13: per-RUN resumable executor for single-(eta, N, M) speed-of-sound cells.

Why this exists. A campaign cell runs its repeats one after another inside a single
00ALLINONE process, and the cell is only marked .done when every repeat has finished.
On 2026-09-13 all processes stopped around 04:19 (battery), leaving the ten N = 2500 cells
of the A2 top-up with 4-7 of their 10 repeats complete. Re-running those cells would redo
every completed repeat. This executor runs exactly the missing repeats instead, and it makes
every later cell (the alpha = 2 cells) resumable at the level of single runs.

Exactness. The runner derives the seed of repeat r as speed_sound_run_seed(base, l, m, r),
a SplitMix hash. With one length and one mass per invocation, l = m = 0. run_seed() below
replicates that function and is checked against every seed recorded in run.log, and against
the Seed column of every existing trace, before anything is launched. Each missing repeat is
executed alone with --repeats=1 and --speed-sound-exact-seed=<its seed>. The runner names that
single output run0; the trace is verified (last Time reaches Planned_Duration, Seed column
equals the expected seed) and renamed to run<r>. Its stdout is appended to the cell's run.log
with the run index rewritten, so the strict health contract and the T_i audit are parsed
exactly as for an uninterrupted cell. Existing complete traces are never modified; an
existing trace that fails its check stops the executor.

Equivalence of a lone exact-seed run with the same repeat inside a multi-repeat process is
not assumed: --gate reproduces an already-completed top-up run (r >= 1, so it follows other
runs in its original process) and requires a byte-identical trace and an identical health
count before any resume is allowed.

Power. New runs start only while the machine is on AC power, unless the file given with
--allow-battery-flag exists. Runs already in flight are never interrupted.
"""
import argparse
import csv
import filecmp
import glob
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
HSPIST = os.path.dirname(HERE)
BIN = "./00ALLINONE"
MASK64 = (1 << 64) - 1
MASSES = (50, 200, 500, 1000, 2000)
_print_lock = threading.Lock()
_file_lock = threading.Lock()


def say(msg):
    with _print_lock:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S %Z')}] {msg}", flush=True)


def run_seed(base, l, m, r):
    """Exact replica of speed_sound_run_seed() in 00ALLINONE.c."""
    x = base & MASK64
    x = (x + 0x9E3779B97F4A7C15 * (l + 1)) & MASK64
    x = (x + 0xBF58476D1CE4E5B9 * (m + 1)) & MASK64
    x = (x + 0x94D049BB133111EB * (r + 1)) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    x ^= x >> 31
    return (x ^ (x >> 32)) & 0xFFFFFFFF


def k_root(alpha):
    lo, hi = 1e-12, math.pi / 2 - 1e-12
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if math.cos(mid) / math.sin(mid) - alpha * mid > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def geometry(eta, N):
    """Same arithmetic and formatting as the campaign shell scripts."""
    fac = math.sqrt(N / 100)
    return f"{3.926990816987241 / eta * fac:.6f}", f"{10 * fac:.6f}"


def eta_tag(eta):
    return "eta_" + ("%.2f" % eta).replace(".", "p")


def rel(c):
    return "/".join(c["dir"].rstrip("/").split("/")[-3:])


# ---------------------------------------------------------------- plans
def plan_topup(root, part):
    """validation/run_A2_topup.sh: seed = 21260912 + 100000*ei + 1000*ni + mi."""
    cells = []
    for ei, eta in enumerate((0.10, 0.30)):
        for ni, N in enumerate((900, 1600, 2500)):
            if part == "b" and N != 2500:
                continue
            if part == "a" and N == 2500:
                continue
            for mi, M in enumerate(MASSES):
                L0, H = geometry(eta, N)
                cells.append(dict(dir=f"{root}/{eta_tag(eta)}/N{N}/m_{M}", eta=eta, N=N, M=M,
                                  L0=L0, H=H, base=21260912 + 100000 * ei + 1000 * ni + mi,
                                  repeats=10 if N == 2500 else 25, max_steps=40000000,
                                  origin="validation/run_A2_topup.sh"))
    return cells


def plan_alpha(root, repeats):
    """validation/run_A2_alpha2.sh: seed = 21860912 + i, i in script loop order."""
    cells, i = [], 0
    for eta in (0.10, 0.30):
        for N, M in ((1600, 3200), (2500, 5000)):
            L0, H = geometry(eta, N)
            cells.append(dict(dir=f"{root}/{eta_tag(eta)}/N{N}/m_{M}", eta=eta, N=N, M=M,
                              L0=L0, H=H, base=21860912 + i, repeats=repeats,
                              max_steps=80000000, origin="validation/run_A2_alpha2.sh"))
            i += 1
    return cells


def argv_for(c, run_dir, repeats, exact_seed=None):
    a = [BIN, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1",
         "--seed-drift-order=drift-first", "--edmd-acc=0",
         f"--particles={c['N']}", f"--particles-boxes={c['N'] // 2},{c['N'] // 2}",
         f"--height={c['H']}", "--particle-radius=0.5",
         "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
         f"--lengths={c['L0']}", f"--wall-masses={c['M']}", f"--repeats={repeats}",
         f"--seed={c['base']}",
         "--wall-hold-steps=2000", "--fixed-dt=0.4",
         "--target-oscillations=25", "--oscillation-safety=1.5",
         "--oscillation-min-steps=10000", f"--oscillation-max-steps={c['max_steps']}",
         "--speed-sound-log-stride=auto", f"--speed-sound-run-dir={run_dir}"]
    if exact_seed is not None:
        a.append(f"--speed-sound-exact-seed={exact_seed}")
    return a


# ---------------------------------------------------------------- checks
def trace_path(cdir, M, r):
    g = glob.glob(os.path.join(cdir, f"wall_x_positions_L0_*_wallmassfactor_{M}_run{r}.csv"))
    if not g:
        return None
    return g[0] if len(g) == 1 else "AMBIGUOUS"


def trace_ok(path, expect_seed=None):
    try:
        with open(path, newline="") as fh:
            rd = csv.reader(fh)
            hdr = next(rd)
            first = next(rd)
            last = first
            for row in rd:
                if row:
                    last = row
        it, ipd, isd = hdr.index("Time"), hdr.index("Planned_Duration"), hdr.index("Seed")
        planned, lt = float(first[ipd]), float(last[it])
        if not (lt >= 0.999 * planned):
            return False, f"incomplete: last Time {lt} < Planned_Duration {planned}"
        if expect_seed is not None and int(float(first[isd])) != int(expect_seed):
            return False, f"Seed column {first[isd]} != expected {expect_seed}"
        return True, "ok"
    except Exception as e:  # noqa: BLE001
        return False, f"unreadable ({e})"


def norm_tok(t):
    if t.startswith("--speed-sound-run-dir="):
        return "--speed-sound-run-dir=" + os.path.normpath(os.path.abspath(t.split("=", 1)[1]))
    return t


def verify_cell(c):
    problems = []
    log = os.path.join(c["dir"], "run.log")
    if os.path.exists(log):
        txt = open(log, errors="replace").read()
        for M, r, s in re.findall(r"Running: L0 = [\d.]+, M = (\d+)\*m, run = (\d+), seed = (\d+)", txt):
            if int(M) != c["M"]:
                problems.append(f"logged mass {M} != {c['M']}")
            if run_seed(c["base"], 0, 0, int(r)) != int(s):
                problems.append(f"run {r}: logged seed {s} != replicated {run_seed(c['base'], 0, 0, int(r))}")
    cmdf = os.path.join(c["dir"], "00_COMMAND.md")
    if os.path.exists(cmdf):
        parts = open(cmdf, errors="replace").read().split("```sh", 1)
        if len(parts) == 2:
            line = next((ln for ln in parts[1].splitlines() if ln.strip() and not ln.startswith("```")), "")
            rec = [norm_tok(t) for t in shlex.split(line)]
            exp = [norm_tok(t) for t in argv_for(c, c["dir"], c["repeats"])]
            if rec != exp:
                diff = [(a, b) for a, b in zip(rec, exp) if a != b]
                problems.append(f"recorded command differs from template: {diff[:3]} "
                                f"(tokens {len(rec)} vs {len(exp)})")
    return problems


def missing_runs(c):
    miss = []
    for r in range(c["repeats"]):
        p = trace_path(c["dir"], c["M"], r)
        if p == "AMBIGUOUS":
            raise SystemExit(f"ambiguous trace for run {r} in {c['dir']}")
        if p is None:
            miss.append(r)
            continue
        ok, why = trace_ok(p, run_seed(c["base"], 0, 0, r))
        if not ok:
            raise SystemExit(f"existing trace {p} fails its check ({why}); refusing to touch it")
    return miss


def cost(c):
    """Relative cost of one run: particles x steps, steps ~ L_eff / K(alpha)."""
    return c["N"] * (float(c["L0"]) - 1.0) / k_root(c["M"] / c["N"])


# ---------------------------------------------------------------- power
def on_ac():
    try:
        out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True, timeout=20).stdout
        return "AC Power" in out
    except Exception:  # noqa: BLE001
        return True


def wait_power(flag):
    waited = False
    while not (on_ac() or (flag and os.path.exists(flag))):
        if not waited:
            say(f"on battery -- holding new runs until AC power (or: touch {flag})")
            waited = True
        time.sleep(60)
    if waited:
        say("power OK, continuing")


# ---------------------------------------------------------------- execution
def execute(c, r, flag, results):
    seed = run_seed(c["base"], 0, 0, r)
    wait_power(flag)
    tmp = os.path.join(c["dir"], f".resume_run{r}")
    if os.path.isdir(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    say(f"start  {rel(c)} run {r} seed {seed}")
    t0 = time.time()
    with open(os.path.join(tmp, "stdout.log"), "w") as out:
        rc = subprocess.call(argv_for(c, tmp, 1, seed), cwd=HSPIST, stdout=out,
                             stderr=subprocess.STDOUT, env=dict(os.environ, HD_KE_TRACE="1"))
    dt = time.time() - t0
    tr = glob.glob(os.path.join(tmp, f"wall_x_positions_L0_*_wallmassfactor_{c['M']}_run0.csv"))
    ok, why = False, f"rc={rc}, run0 traces={len(tr)}"
    if rc == 0 and len(tr) == 1:
        ok, why = trace_ok(tr[0], seed)
    if not ok:
        dead = os.path.join(c["dir"], f".failed_run{r}_{time.strftime('%Y%m%d_%H%M%S')}")
        os.rename(tmp, dead)
        say(f"FAIL   {rel(c)} run {r}: {why} (kept in {os.path.basename(dead)})")
        results.append((c["dir"], r, "fail", why, dt, 0))
        return
    final = os.path.join(c["dir"], os.path.basename(tr[0])[: -len("run0.csv")] + f"run{r}.csv")
    with _file_lock:
        if os.path.exists(final):
            say(f"FAIL   {rel(c)} run {r}: {os.path.basename(final)} appeared meanwhile; not overwriting")
            results.append((c["dir"], r, "fail", "target exists", dt, 0))
            return
        os.rename(tr[0], final)
        txt = open(os.path.join(tmp, "stdout.log"), errors="replace").read()
        txt = re.sub(r"run = 0,", f"run = {r},", txt)
        txt = re.sub(r"run=0 seed=", f"run={r} seed=", txt)
        health = len(re.findall(r"EDMD-HEALTH", txt))
        stamp = time.strftime("%Y-%m-%d %H:%M:%S %Z")
        with open(os.path.join(c["dir"], "run.log"), "a") as fh:
            fh.write(f"\n##RESUME {stamp} run {r} seed {seed} (exact-seed single run, {dt:.0f} s)\n")
            fh.write(txt)
        p6 = os.path.join(tmp, "speed_of_sound_psi6.csv")
        cell_p6 = os.path.join(c["dir"], "speed_of_sound_psi6.csv")
        if os.path.exists(p6):
            lines = open(p6).read().splitlines()
            if lines:
                if not os.path.exists(cell_p6):
                    with open(cell_p6, "w") as fh:
                        fh.write(lines[0] + "\n")
                with open(cell_p6, "a") as fh:
                    for ln in lines[1:]:
                        fh.write(ln + "\n")
        cmdf = os.path.join(c["dir"], "00_COMMAND.md")
        if not os.path.exists(cmdf):
            with open(cmdf, "w") as fh:
                fh.write("# Command\n\n"
                         "Equivalent single invocation of this cell. It was executed run by run with\n"
                         "`--repeats=1 --speed-sound-exact-seed=<seed>` by validation/resume_runs.py;\n"
                         "see 00_RESUME.md for the per-run seeds.\n\n```sh\n"
                         + " ".join(argv_for(c, c["dir"], c["repeats"])) + "\n```\n")
        rpath = os.path.join(c["dir"], "00_RESUME.md")
        new = not os.path.exists(rpath)
        with open(rpath, "a") as fh:
            if new:
                fh.write(f"# Per-run execution record\n\n"
                         f"Origin: `{c['origin']}`. Base seed {c['base']}; run r uses "
                         f"speed_sound_run_seed(base, l=0, m=0, r), replicated in validation/resume_runs.py and "
                         f"verified against run.log. Each run below was executed alone with --repeats=1 and "
                         f"--speed-sound-exact-seed, renamed run0 -> run<r>, and its log appended to run.log "
                         f"with the run index rewritten.\n\n"
                         f"| run | seed | seconds | health lines | finished |\n|---|---|---|---|---|\n")
            fh.write(f"| {r} | {seed} | {dt:.0f} | {health} | {stamp} |\n")
    shutil.rmtree(tmp, ignore_errors=True)
    say(f"done   {rel(c)} run {r}  {dt / 60:.0f} min  health lines={health}")
    results.append((c["dir"], r, "ok", "", dt, health))


def safe_execute(c, r, flag, results):
    try:
        execute(c, r, flag, results)
    except Exception as e:  # noqa: BLE001
        say(f"FAIL   {rel(c)} run {r}: exception {e!r}")
        results.append((c["dir"], r, "fail", repr(e), 0.0, 0))


# ---------------------------------------------------------------- gate
def gate(cell_dir, r, out_dir):
    cell_dir = os.path.normpath(os.path.abspath(cell_dir))
    root = os.path.dirname(os.path.dirname(os.path.dirname(cell_dir)))
    c = next((x for x in plan_topup(root, "all") if os.path.normpath(x["dir"]) == cell_dir), None)
    if c is None:
        say(f"GATE: {cell_dir} is not a top-up cell")
        return 2
    probs = verify_cell(c)
    if probs:
        say("GATE: provenance problems: " + "; ".join(probs))
        return 2
    seed = run_seed(c["base"], 0, 0, r)
    orig = trace_path(c["dir"], c["M"], r)
    if orig in (None, "AMBIGUOUS"):
        say(f"GATE: no unique original trace for run {r}")
        return 2
    ok, why = trace_ok(orig, seed)
    if not ok:
        say(f"GATE: original trace fails check: {why}")
        return 2
    out_dir = os.path.abspath(out_dir)
    shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir)
    say(f"GATE: re-running {rel(c)} run {r} alone with exact seed {seed}")
    t0 = time.time()
    with open(os.path.join(out_dir, "stdout.log"), "w") as out:
        rc = subprocess.call(argv_for(c, out_dir, 1, seed), cwd=HSPIST, stdout=out,
                             stderr=subprocess.STDOUT, env=dict(os.environ, HD_KE_TRACE="1"))
    new = glob.glob(os.path.join(out_dir, f"wall_x_positions_L0_*_wallmassfactor_{c['M']}_run0.csv"))
    same = rc == 0 and len(new) == 1 and filecmp.cmp(orig, new[0], shallow=False)
    olog = open(os.path.join(c["dir"], "run.log"), errors="replace").read()
    oh = len(re.findall(rf"EDMD-HEALTH\] L0=[\d.]+ M={c['M']} run={r} seed={seed}:", olog))
    nh = len(re.findall(r"EDMD-HEALTH", open(os.path.join(out_dir, "stdout.log"), errors="replace").read()))
    say(f"GATE: rc={rc}  {time.time() - t0:.0f} s  trace byte-identical={same}  "
        f"health lines original={oh} reproduction={nh}")
    passed = same and oh == nh
    say("GATE PASS" if passed else "GATE FAIL")
    return 0 if passed else 1


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", choices=["topup-a", "topup-b", "alpha"])
    ap.add_argument("--root")
    ap.add_argument("--repeats", type=int, default=10, help="alpha plan only")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--allow-battery-flag", default="")
    ap.add_argument("--plan-only", action="store_true")
    ap.add_argument("--gate")
    ap.add_argument("--gate-run", type=int, default=1)
    ap.add_argument("--gate-out")
    a = ap.parse_args()
    os.chdir(HSPIST)
    if a.gate:
        sys.exit(gate(a.gate, a.gate_run, a.gate_out))
    if not (a.plan and a.root):
        ap.error("--plan and --root are required unless --gate is given")
    root = os.path.abspath(a.root)
    cells = plan_topup(root, a.plan[-1]) if a.plan.startswith("topup") else plan_alpha(root, a.repeats)

    bad, jobs = False, []
    for c in cells:
        probs = verify_cell(c)
        if probs:
            bad = True
            say(f"PROVENANCE {rel(c)}: " + "; ".join(probs))
        miss = missing_runs(c)
        say(f"{rel(c)}: base {c['base']}  complete {c['repeats'] - len(miss)}/{c['repeats']}  missing {miss}")
        jobs += [(cost(c), c, r) for r in miss]
    if bad:
        raise SystemExit("provenance check failed; nothing launched")
    jobs.sort(key=lambda j: -j[0])
    say(f"{len(jobs)} runs to execute, {a.jobs} at a time, longest first")
    if a.plan_only or not jobs:
        return

    results = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(safe_execute, c, r, a.allow_battery_flag, results) for _, c, r in jobs]
        for fu in futs:
            fu.result()

    for c in cells:
        if os.path.isdir(c["dir"]) and not missing_runs(c):
            open(os.path.join(c["dir"], ".done"), "a").close()
    fails = [x for x in results if x[2] != "ok"]
    health = sum(x[5] for x in results)
    say(f"finished: {len(results) - len(fails)} ok, {len(fails)} failed, {health} health lines in new runs; "
        f"cells complete {sum(1 for c in cells if os.path.exists(os.path.join(c['dir'], '.done')))}/{len(cells)}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
