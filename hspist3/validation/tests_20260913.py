#!/usr/bin/env python3
"""##CHRIS 2026-09-13: overnight tests A (record length), B (settling time), C (stride) and the
conditional launch D (A1 v2), as specified in the go of 2026-09-13. Deadline for the report:
2026-09-14 08:00 HST.

Estimator (0000_PLAN_OVERALL/ALL_MARKDOWNS/260913_method_roman_fft_COWORK.md): per trajectory,
mean-subtracted FFT of the divider position from release, frequency = centre of the largest bin at
f > 0, no transient drop, no cuts. Per mass: mean over seeds (primary), median as an extra column.
c_s = slope of nu against x_M = K(M/2N_side)/(2 pi (L0 - 2r)) through the origin; error = scatter
(1 sigma) of the per-mass values nu_M/x_M. Only discard: the strict health contract.
Extra estimators, reported beside it: 3-bin parabolic interpolation of log P around the maximum
(no interpolation when the maximum is bin 1, whose lower neighbour is the removed mean), and a
Whittle-likelihood fit of the driven damped-oscillator line shape S(u) = A/((u^2-u0^2)^2+(g u)^2)+B
to the seed-averaged periodogram, in the band 0.5-2x the median per-seed peak.

Execution: every trajectory is one 00ALLINONE process with --repeats=1 and
--speed-sound-exact-seed. Seeds follow the single-invocation convention
speed_sound_run_seed(base, l=0, m=mass index, r). Test B uses the SAME seeds at every hold, so the
hold comparison is paired. Runs are resumable one by one; new runs start only on AC power or with
the battery at or above BATTERY_FLOOR percent. A and B record exactly the target number of
predicted oscillations (--oscillation-safety=1.0), so each analysis cut is a prefix of the record.
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HSP = os.path.dirname(HERE)
REPO = os.path.dirname(HSP)
sys.path.insert(0, HSP)

ROOT = os.path.join(HSP, "experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN")
TROOT = os.path.join(ROOT, "tests_20260913")
DROOT = os.path.join(ROOT, "A1v2_20260914")
STATE = os.path.join(ROOT, "_orchestration_tests_20260913")
MD = os.path.join(REPO, "0000_PLAN_OVERALL", "ALL_MARKDOWNS")
# ##CHRIS 2026-09-19: output moved into the per-paper folders. PLOTS is now where new figures and
# CSVs are WRITTEN; 260909_plots is frozen and still holds every historical file, so reads fall back
# to it (and to the sorted copies) through plot_path() below. Nothing was moved or deleted.
PLOTS_FROZEN = os.path.join(MD, "260909_plots")
PAPER1 = os.path.join(os.path.dirname(MD), "paper1_speedofsound", "experiments")
PAPER2 = os.path.join(os.path.dirname(MD), "paper2_energytransfer", "experiments")
PLOTS = os.path.join(PAPER1, "final")
os.makedirs(PLOTS, exist_ok=True)


def plot_path(name, write=False):
    """Resolve a figure/CSV by name. Writes go to PLOTS; reads search the sorted copies and then
    the frozen archive, so scripts written before 2026-09-19 keep working unchanged."""
    if write:
        return os.path.join(PLOTS, name)
    for d in (PLOTS, os.path.join(PAPER1, "estimator_tests"), os.path.join(PAPER1, "archive"),
              os.path.join(PAPER2, "final"), PLOTS_FROZEN):
        q = os.path.join(d, name)
        if os.path.exists(q):
            return q
    return os.path.join(PLOTS, name)
STATUS = os.path.join(MD, "260913_tests_STATUS.md")
REPORT = os.path.join(MD, "260913_tests_REPORT.md")
A1_DIR = os.path.join(ROOT, "campaign_r25_psi6_20260823")
LOWETA_DIR = os.path.join(ROOT, "routeA_lowdensity_20260912")

N_TOTAL, N_SIDE, RDISK, H = 100, 50, 0.5, "10.0"
# ##CHRIS 2026-09-18: ONE source of truth for the divider thickness. It is passed to the binary in
# argv_for() and consumed by x_of() below, so the launcher and the estimator cannot disagree -- which
# is exactly how L_eff came to be wrong (the estimator assumed a zero-thickness divider while every
# run had one 0.05 sigma thick). A disk centre in one compartment spans r .. L0 - t/2 - r, so
# L_eff = L0 - 2r - t/2. Changing this constant changes both the runs and the analysis together.
WALL_T = 0.05
DT_SIGMA = 0.4 / 24.0
HOLD_NOW = 2000
ROMAN_MASSES = [20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
A1_MASSES = [50, 100, 200, 300, 500, 750, 1000, 1500, 2000]
B_MASSES = [100, 200, 500, 1000, 2000]
HOLDS = [2000, 20000, 200000]
CUTS = [25, 50, 100, 200, 500, 1000]
L0S = [("20.0", 20.0), ("7.5", 7.5)]
ROMAN_TABLE = {20.0: (2.20, 0.02), 7.5: (5.99, 0.09)}
SEEDS, JOBS = 25, 10
A_BASE, B_BASE, D_BASE = 22060913, 22160913, 22260914
MAX_STEPS = 400000000
BATTERY_FLOOR = 35
D_SAMPLES_PER_PERIOD = 32
D_TARGET_CAP = 200
D_WALL_LIMIT_S = 10 * 3600
D_FIRST_ETAS = 15
MASK64 = (1 << 64) - 1
HEALTH_RE = re.compile(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): forced_advance=(\d+) "
                       r"wall_clamp_repairs=(\d+) overlap_repairs=(\d+) wall_overdue=(\d+)")
COUNTERS = ("forced_advance", "clamp_repair", "overlap_repair", "wall_overdue")

_lock = threading.Lock()
QUIET = False


# ------------------------------------------------------------------ small utilities
def now():
    return time.strftime("%Y-%m-%d %H:%M:%S %Z")


def status(msg):
    if QUIET:
        print(f"[selftest] {msg}", flush=True)
        return
    with _lock:
        new = not os.path.exists(STATUS)
        with open(STATUS, "a") as fh:
            if new:
                fh.write("# Overnight tests 2026-09-13: status log\n\n"
                         "One line per step, newest at the bottom. Written automatically by "
                         "`hspist3/validation/tests_20260913.py`. The full report is "
                         "`260913_tests_REPORT.md` in this folder.\n\n")
            fh.write(f"- **{now()}** {msg}\n")
        print(f"[{now()}] {msg}", flush=True)


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


def l_eff(L0, wall_t=None):
    """Centre-accessible length of one compartment: L0 - 2r - t/2 (##CHRIS 2026-09-18)."""
    t = WALL_T if wall_t is None else wall_t
    return L0 - 2 * RDISK - 0.5 * t


def assert_wall_thickness(run_dir):
    """##CHRIS 2026-09-18: refuse to analyse data whose recorded thickness contradicts WALL_T.

    The failure that motivated this was silent: the analysis assumed one geometry and the data had
    another, and nothing in either file said so. Where a campaign records its command (00_COMMAND.md,
    or a run.log that carries the argv) this reads the flag back and raises on a mismatch. Campaigns
    that record nothing return None -- absence of evidence, reported as such, never assumed to agree.
    """
    import glob as _glob, re as _re
    for name in ("00_COMMAND.md", "run.log"):
        for f in _glob.glob(os.path.join(run_dir, "**", name), recursive=True)[:1]:
            m = _re.search(r"--wall-thickness=([0-9.]+)", open(f, errors="replace").read())
            if m:
                got = float(m.group(1))
                if abs(got - WALL_T) > 1e-12:
                    raise SystemExit(f"wall thickness mismatch: {run_dir} ran with t = {got}, "
                                     f"the estimator is built for t = {WALL_T}. Fix one or the other.")
                return got
    return None


def x_of(M, L0, wall_t=None):
    return k_root(M / (2.0 * N_SIDE)) / (2 * math.pi * l_eff(L0, wall_t))


def kr_cs(eta):
    import plot_speed_of_sound_edmd as sos
    a = np.array([eta]); h = 1e-5
    Z = sos.Z_kolafa_rottner_2006(a)
    dZ = (sos.Z_kolafa_rottner_2006(a + h) - sos.Z_kolafa_rottner_2006(a - h)) / (2 * h)
    return float(sos.cs_adiabatic_2d_monatomic(Z, dZ, a, kbt=1, m=1)[0])


def eta_of(L0):
    return N_TOTAL * math.pi * RDISK * RDISK / (2 * L0 * float(H))


def lpt(durations, slots):
    load = [0.0] * slots
    for d in sorted(durations, reverse=True):
        i = load.index(min(load)); load[i] += d
    return max(load) if load else 0.0


def tag(s):
    return s.replace(".", "p")


# ------------------------------------------------------------------ geometry and timing from A1
def a1_leaf_table():
    """Per A1 density: L0 string, predicted frequency at M = 50, seconds per integration step
    measured on the A1 leaf (leaf wall time / total planned steps incl. hold)."""
    rows = []
    leaves = sorted(glob.glob(f"{A1_DIR}/eta_*")) + sorted(glob.glob(f"{LOWETA_DIR}/eta_*"))
    for d in leaves:
        if not os.path.isdir(d) or not os.path.exists(f"{d}/speed_of_sound_psi6.csv"):
            continue
        cmd = open(f"{d}/00_COMMAND.md", errors="replace").read()
        L0s = re.search(r"--lengths=([\d.]+)", cmd).group(1)
        reps = int(re.search(r"--repeats=(\d+)", cmd).group(1))
        f50 = sorted(glob.glob(f"{d}/wall_x_positions_*wallmassfactor_50_run0.csv"))[0]
        with open(f50) as fh:
            r0 = next(csv.DictReader(fh))
        nu50 = float(r0["Predicted_Frequency"]); eta = float(r0["eta"])
        periods = float(r0["Planned_Duration"]) * nu50
        k50 = k_root(0.5)
        total_steps = 0.0
        for M in A1_MASSES:
            nup = nu50 * k_root(M / 100.0) / k50
            total_steps += reps * (periods / (nup * DT_SIGMA) + HOLD_NOW)
        wall = os.path.getmtime(f"{d}/speed_of_sound_psi6.csv") - os.stat(d).st_birthtime
        rows.append(dict(leaf=os.path.basename(d), eta=eta, L0=L0s, nu50=nu50,
                         s_per_step=wall / total_steps, runs=reps * len(A1_MASSES), wall=wall))
    rows.sort(key=lambda z: z["eta"])
    return rows


def nu_pred(nu50, M):
    return nu50 * k_root(M / 100.0) / k_root(0.5)


def leaf_for_L0(table, L0):
    return min(table, key=lambda z: abs(float(z["L0"]) - L0))


# ------------------------------------------------------------------ running trajectories
def power_ok():
    try:
        out = subprocess.run(["pmset", "-g", "batt"], capture_output=True, text=True, timeout=20).stdout
    except Exception:  # noqa: BLE001
        return True, "pmset unavailable"
    if "AC Power" in out:
        return True, "AC power"
    m = re.search(r"(\d+)%", out)
    pct = int(m.group(1)) if m else 100
    return pct >= BATTERY_FLOOR, f"battery {pct}%"


_power_waiting = [False]


def gate_power():
    while True:
        ok, why = power_ok()
        with _lock:
            waiting = _power_waiting[0]
        if ok:
            if waiting:
                with _lock:
                    _power_waiting[0] = False
                status(f"power OK again ({why}); new runs continue")
            return
        if not waiting:
            with _lock:
                _power_waiting[0] = True
            status(f"PAUSED starting new runs: {why}, below the {BATTERY_FLOOR}% floor and not on AC. "
                   f"Runs resume by themselves once the charger is connected.")
        time.sleep(60)


def argv_for(job, run_dir):
    return ["./00ALLINONE", "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1",
            "--seed-drift-order=drift-first", "--edmd-acc=0",
            f"--particles={N_TOTAL}", f"--particles-boxes={N_SIDE},{N_SIDE}",
            f"--height={H}", f"--particle-radius={RDISK}",
            f"--wall-thickness={WALL_T}", f"--wall-thickness-vis={WALL_T}",
            f"--lengths={job['L0']}", f"--wall-masses={job['M']}", "--repeats=1",
            f"--seed={job['base']}", f"--wall-hold-steps={job['hold']}", "--fixed-dt=0.4",
            f"--target-oscillations={job['target']}", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", f"--oscillation-max-steps={MAX_STEPS}",
            f"--speed-sound-log-stride={job['stride']}", f"--speed-sound-run-dir={run_dir}",
            f"--speed-sound-exact-seed={job['seed']}"]


def trace_check(path, seed=None):
    try:
        with open(path, newline="") as fh:
            rd = csv.reader(fh)
            hdr = next(rd); first = next(rd)
        with open(path, "rb") as fh:
            fh.seek(0, 2); size = fh.tell(); fh.seek(max(0, size - 8192))
            last = fh.read().decode(errors="replace").strip().splitlines()[-1].split(",")
        it, ipd, isd = hdr.index("Time"), hdr.index("Planned_Duration"), hdr.index("Seed")
        planned, lt = float(first[ipd]), float(last[it])
        if not lt >= 0.999 * planned:
            return False, f"incomplete: last Time {lt} < {planned}"
        if seed is not None and int(float(first[isd])) != int(seed):
            return False, f"Seed column {first[isd]} != {seed}"
        return True, "ok"
    except Exception as e:  # noqa: BLE001
        return False, f"unreadable ({e})"


def final_trace(job):
    g = glob.glob(os.path.join(job["cell"], f"wall_x_positions_L0_*_wallmassfactor_{job['M']}_run{job['r']}.csv"))
    return g[0] if len(g) == 1 else None


def execute(job, ledger):
    existing = final_trace(job)
    if existing and trace_check(existing, job["seed"])[0]:
        return "skip"
    gate_power()
    cell, r, M = job["cell"], job["r"], job["M"]
    tmp = os.path.join(cell, f".run{r}")
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    t0 = time.time()
    argv = argv_for(job, tmp)
    # ##CHRIS 2026-09-18: record the command beside the data. The speed-of-sound runs wrote only a
    # run.log, which does not carry the argv, so the geometry a campaign actually used could not be
    # read back from its own directory -- that is what made the L_eff error invisible for a week.
    # The energy-transfer harness already writes this file; now the Paper 1 harness does too.
    with open(os.path.join(cell, "00_COMMAND.md"), "w") as fh:
        fh.write(f"# Command\n\n- Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n"
                 f"- CWD: `{HSP}`\n- Estimator geometry: L_eff = L0 - 2r - t/2 with "
                 f"r = {RDISK}, t = {WALL_T}\n\n```sh\n" + shlex.join(argv) + "\n```\n")
    with open(os.path.join(tmp, "stdout.log"), "w") as out:
        rc = subprocess.call(argv, cwd=HSP, stdout=out, stderr=subprocess.STDOUT,
                             env=dict(os.environ, HD_KE_TRACE="1"))
    dt = time.time() - t0
    txt = open(os.path.join(tmp, "stdout.log"), errors="replace").read()
    counters = dict.fromkeys(COUNTERS, 0)
    m = HEALTH_RE.search(txt)
    if m:
        counters = dict(zip(COUNTERS, (int(m.group(i)) for i in (4, 5, 6, 7))))
    tr = glob.glob(os.path.join(tmp, f"wall_x_positions_L0_*_wallmassfactor_{M}_run0.csv"))
    ok, why = False, f"rc={rc}, traces={len(tr)}"
    if rc == 0 and len(tr) == 1:
        ok, why = trace_check(tr[0], job["seed"])
    rec = dict(test=job["test"], cell=os.path.relpath(cell, ROOT), L0=job["L0"], M=M, r=r, seed=job["seed"],
               hold=job["hold"], target=job["target"], stride=job["stride"], rc=rc, ok=ok, why=why,
               seconds=round(dt, 1), finished=now(), **counters)
    if ok:
        final = os.path.join(cell, os.path.basename(tr[0])[: -len("run0.csv")] + f"run{r}.csv")
        with _lock:
            os.replace(tr[0], final)
            txt2 = re.sub(r"run = 0,", f"run = {r},", txt)
            txt2 = re.sub(r"run=0 seed=", f"run={r} seed=", txt2)
            with open(os.path.join(cell, "run.log"), "a") as fh:
                fh.write(f"\n##RUN {now()} run {r} seed {job['seed']} ({dt:.0f} s)\n")
                fh.write(txt2)
            cmdf = os.path.join(STATE, f"command_{job['test']}.txt")
            if not os.path.exists(cmdf):
                with open(cmdf, "w") as fh:
                    fh.write("HD_KE_TRACE=1 " + shlex.join(argv_for(job, cell)) + "\n")
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        os.rename(tmp, os.path.join(cell, f".failed_run{r}_{time.strftime('%Y%m%d_%H%M%S')}"))
    with _lock:
        with open(ledger, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
    return "ok" if ok else "fail"


def run_jobs(test, jobs, heartbeat=None):
    ledger = os.path.join(STATE, f"ledger_{test}.jsonl")
    done = {"ok": 0, "fail": 0, "skip": 0}
    stop = threading.Event()
    t_start = time.time()

    def beat():
        while not stop.wait(1800):
            status(f"{test} progress: {done['ok'] + done['skip']}/{len(jobs)} trajectories complete, "
                   f"{done['fail']} failed, {(time.time() - t_start) / 3600:.1f} h elapsed")
            if heartbeat:
                try:
                    heartbeat()
                except Exception as e:  # noqa: BLE001
                    status(f"report refresh failed: {e!r}")

    threading.Thread(target=beat, daemon=True).start()

    def one(job):
        try:
            res = execute(job, ledger)
        except Exception as e:  # noqa: BLE001
            res = "fail"
            status(f"{test}: exception on {os.path.relpath(job['cell'], ROOT)} run {job['r']}: {e!r}")
        with _lock:
            done[res] += 1
        return res

    with ThreadPoolExecutor(max_workers=JOBS) as ex:
        list(ex.map(one, jobs))
    stop.set()
    return done


# ------------------------------------------------------------------ job lists
def cost_seconds(table, L0, M, periods, hold):
    leaf = leaf_for_L0(table, L0)
    nup = nu_pred(leaf["nu50"], M)
    return (periods / (nup * DT_SIGMA) + hold) * leaf["s_per_step"]


def jobs_A(table):
    J = []
    for li, (L0s, L0) in enumerate(L0S):
        for mi, M in enumerate(ROMAN_MASSES):
            base = A_BASE + li
            for r in range(SEEDS):
                J.append(dict(test="A", cell=os.path.join(TROOT, "A_length", f"L0_{tag(L0s)}", f"m_{M}"),
                              L0=L0s, L0f=L0, M=M, mi=mi, r=r, base=base, seed=run_seed(base, 0, mi, r),
                              hold=HOLD_NOW, target=1000, stride="auto",
                              cost=cost_seconds(table, L0, M, 1000, HOLD_NOW)))
    first = [j for j in J if j["L0"] == "20.0" and j["M"] == 500 and j["r"] < 10]
    rest = sorted([j for j in J if not (j["L0"] == "20.0" and j["M"] == 500 and j["r"] < 10)], key=lambda j: -j["cost"])
    return first + rest


def jobs_B(table, LB):
    J = []
    for li, (L0s, L0) in enumerate(L0S):
        for hold in HOLDS:
            for mi, M in enumerate(B_MASSES):
                base = B_BASE + li
                for r in range(SEEDS):
                    J.append(dict(test="B", cell=os.path.join(TROOT, "B_hold", f"L0_{tag(L0s)}", f"hold_{hold}", f"m_{M}"),
                                  L0=L0s, L0f=L0, M=M, mi=mi, r=r, base=base, seed=run_seed(base, 0, mi, r),
                                  hold=hold, target=LB, stride="auto",
                                  cost=cost_seconds(table, L0, M, LB, hold)))
    return sorted(J, key=lambda j: -j["cost"])


def d_stride(nup):
    return max(1, int(math.floor((1.0 / (nup * DT_SIGMA)) / D_SAMPLES_PER_PERIOD)))


def jobs_D(table, TD):
    J = []
    for ei, leaf in enumerate(table):
        L0 = float(leaf["L0"])
        for mi, M in enumerate(A1_MASSES):
            base = D_BASE + ei
            nup = nu_pred(leaf["nu50"], M)
            steps = TD / (nup * DT_SIGMA)
            for r in range(SEEDS):
                J.append(dict(test="D", cell=os.path.join(DROOT, leaf["leaf"], f"m_{M}"), L0=leaf["L0"], L0f=L0,
                              M=M, mi=mi, r=r, base=base, seed=run_seed(base, 0, mi, r), hold=HOLD_NOW,
                              target=TD, stride=d_stride(nup), eta=leaf["eta"], ei=ei,
                              cost=(steps + HOLD_NOW) * leaf["s_per_step"], rows=steps / d_stride(nup)))
    return J


# ------------------------------------------------------------------ analysis (worker side)
def _load(path):
    import pandas as pd
    with open(path, newline="") as fh:
        r0 = next(csv.DictReader(fh))
    d = pd.read_csv(path, usecols=["Time", "Displacement(σ)"])
    return d["Time"].to_numpy(float), d["Displacement(σ)"].to_numpy(float), float(r0["Predicted_Frequency"])


def _prefix(t, nup, N):
    Tc = N / nup
    if t[-1] < 0.999 * Tc:
        return None
    return int(np.searchsorted(t, Tc * (1 + 1e-9), side="right"))


def _spectrum(x, dt):
    xx = x - x.mean()
    return np.abs(np.fft.rfft(xx)) ** 2, 1.0 / (len(xx) * dt)


def _peak(P, df):
    k = 1 + int(np.argmax(P[1:]))
    delta = 0.0
    if 2 <= k < len(P) - 1 and P[k - 1] > 0 and P[k + 1] > 0:
        a, b, c = math.log(P[k - 1]), math.log(P[k]), math.log(P[k + 1])
        den = a - 2 * b + c
        if den < 0:
            delta = max(-0.5, min(0.5, 0.5 * (a - c) / den))
    return k, k * df, (k + delta) * df


def whittle_center(Ibar, df, fc, R):
    from scipy.optimize import minimize
    idx = np.arange(len(Ibar))
    f = idx * df
    sel = (idx >= 1) & (f >= 0.5 * fc) & (f <= 2.0 * fc)
    if sel.sum() < 8 or not (fc > 0):
        return float("nan")
    u = f[sel] / fc
    I = Ibar[sel] / Ibar[sel].max()
    B0 = max(float(np.percentile(I, 10)), 1e-8)
    g0 = 0.08
    A0 = max(float(I.max()) - B0, 1e-8) * g0 * g0

    def nll(p):
        lA, u0, lg, lB = p
        if not (0.5 < u0 < 2.0):
            return 1e30
        S = math.exp(lA) / ((u * u - u0 * u0) ** 2 + (math.exp(lg) * u) ** 2) + math.exp(lB)
        if not np.all(np.isfinite(S)) or np.any(S <= 0):
            return 1e30
        return R * float(np.sum(np.log(S) + I / S))

    best = None
    for u0s in (1.0, float(u[int(np.argmax(I))])):
        res = minimize(nll, [math.log(A0), u0s, math.log(g0), math.log(B0)], method="Nelder-Mead",
                       options=dict(maxiter=8000, xatol=1e-8, fatol=1e-10))
        if best is None or res.fun < best.fun:
            best = res
    u0 = float(best.x[1])
    return u0 * fc if 0.5 < u0 < 2.0 else float("nan")


def analyse_cell(task):
    """task = (key, L0, M, [(r, path, discarded)], cuts). Returns per-cut per-run peaks and the
    seed-averaged spectrum fit. Discarded (health) runs are carried along, marked, never used."""
    key, L0, M, runs, cuts = task
    per_cut = {N: {"r": [], "nu": [], "par": [], "k": [], "disc": [], "P": []} for N in cuts}
    dfs = {}
    for r, path, disc in runs:
        try:
            t, x, nup = _load(path)
        except Exception:  # noqa: BLE001
            continue
        dt = (t[-1] - t[0]) / (len(t) - 1)
        for N in cuts:
            n = _prefix(t, nup, N)
            if n is None or n < 16:
                continue
            P, df = _spectrum(x[:n], dt)
            k, nu, par = _peak(P, df)
            c = per_cut[N]
            c["r"].append(r); c["nu"].append(nu); c["par"].append(par); c["k"].append(k); c["disc"].append(bool(disc))
            if not disc:
                c["P"].append(P); dfs.setdefault(N, df)
    out = {"key": key, "L0": L0, "M": M, "cuts": {}}
    for N in cuts:
        c = per_cut[N]
        used = [i for i, d in enumerate(c["disc"]) if not d]
        nus = [c["nu"][i] for i in used]
        fit = float("nan")
        avg_peak = float("nan")
        if len(c["P"]) >= 1:
            L = min(len(p) for p in c["P"])
            Ibar = np.mean([p[:L] for p in c["P"]], axis=0)
            # peak of the seed-averaged spectrum: average every seed's |FFT|^2 first, then its largest bin at f > 0
            avg_peak = (1 + int(np.argmax(Ibar[1:]))) * dfs[N]
            if len(c["P"]) >= 3:
                fit = whittle_center(Ibar, dfs[N], float(np.median(nus)), len(c["P"]))
        out["cuts"][N] = dict(
            r=c["r"], nu=c["nu"], par=c["par"], k=c["k"], disc=c["disc"],
            n_used=len(used), n_disc=len(c["disc"]) - len(used),
            nu_mean=float(np.mean(nus)) if nus else float("nan"),
            nu_sd=float(np.std(nus, ddof=1)) if len(nus) > 1 else float("nan"),
            nu_median=float(np.median(nus)) if nus else float("nan"),
            par_mean=float(np.mean([c["par"][i] for i in used])) if used else float("nan"),
            nu_mean_incl_disc=float(np.mean(c["nu"])) if c["nu"] else float("nan"),
            drift_frac=float(np.mean([c["k"][i] <= 3 for i in used])) if used else float("nan"),
            fit_f0=fit, avg_peak=avg_peak)
    return out


def analyse(tasks, workers=JOBS):
    with Pool(workers) as p:
        return p.map(analyse_cell, tasks, chunksize=1)


def slope(xs, ys):
    xs = np.array(xs, float); ys = np.array(ys, float)
    ok = np.isfinite(ys) & np.isfinite(xs)
    if ok.sum() < 3:
        return float("nan"), float("nan"), int(ok.sum())
    xs, ys = xs[ok], ys[ok]
    return float((xs * ys).sum() / (xs * xs).sum()), float(np.std(ys / xs, ddof=1)), int(ok.sum())


def health_of(path_log):
    bad = set()
    if os.path.exists(path_log):
        for m in HEALTH_RE.finditer(open(path_log, errors="replace").read()):
            if any(int(m.group(i)) for i in (4, 5, 6, 7)):
                bad.add(int(m.group(2)))
    return bad


def cell_runs(cell, M):
    bad = health_of(os.path.join(cell, "run.log"))
    runs = []
    for p in glob.glob(os.path.join(cell, f"wall_x_positions_L0_*_wallmassfactor_{M}_run*.csv")):
        r = int(re.search(r"_run(\d+)\.csv$", p).group(1))
        if trace_check(p)[0]:
            runs.append((r, p, r in bad))
    return sorted(runs)


# ------------------------------------------------------------------ TEST A analysis + figure
ESTIMATORS = [("roman_mean", "Román: mean of peak bins", "#2a78d6", "o"),
              ("roman_median", "median of peak bins", "#eb6834", "s"),
              ("parabolic", "3-bin parabolic peak", "#1baf7a", "^"),
              ("spectrum_fit", "seed-averaged spectrum fit", "#eda100", "D")]


def summarize(results, L0, cuts, masses):
    rows, per_mass = [], []
    for N in cuts:
        xs, ym, yd, yp, yf, yi, ya = [], [], [], [], [], [], []
        used = disc = 0
        for M in masses:
            res = next((q for q in results if q["L0"] == L0 and q["M"] == M), None)
            if res is None or N not in res["cuts"]:
                continue
            c = res["cuts"][N]
            xs.append(x_of(M, L0)); ym.append(c["nu_mean"]); yd.append(c["nu_median"])
            yp.append(c["par_mean"]); yf.append(c["fit_f0"]); yi.append(c["nu_mean_incl_disc"]); ya.append(c.get("avg_peak", float("nan")))
            used += c["n_used"]; disc += c["n_disc"]
            per_mass.append(dict(L0=L0, cut=N, M=M, n_used=c["n_used"], n_disc=c["n_disc"], nu_mean=c["nu_mean"],
                                 nu_sd=c["nu_sd"], nu_sem=c["nu_sd"] / math.sqrt(c["n_used"]) if c["n_used"] > 1 else float("nan"),
                                 nu_median=c["nu_median"], parabolic_mean=c["par_mean"], spectrum_fit_f0=c["fit_f0"],
                                 drift_fraction_bins_1_3=c["drift_frac"], avg_spectrum_peak=c.get("avg_peak", float("nan")), x_M=x_of(M, L0)))
        cm, sm, nm = slope(xs, ym); cd, sd_, _ = slope(xs, yd); cp, sp, _ = slope(xs, yp)
        cf, sf, nf = slope(xs, yf); ci, _, _ = slope(xs, yi); ca, sa, _ = slope(xs, ya)
        rows.append(dict(L0=L0, cut=N, roman_mean=cm, roman_mean_scatter=sm, roman_median=cd, parabolic=cp,
                         spectrum_fit=cf, spectrum_fit_masses=nf, avg_spectrum_peak=ca, avg_spectrum_peak_scatter=sa,
                         n_masses=nm, trajectories_used=used,
                         trajectories_discarded=disc, roman_mean_incl_discarded=ci))
    return rows, per_mass


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def figure_A(rows, cuts, out_base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    INK, INK2, MUTED, GRID, AXIS, SURF, KRC = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#e34948"
    L0vals = [L0 for _, L0 in L0S]
    fig, axes = plt.subplots(1, len(L0vals), figsize=(7.4 * len(L0vals), 5.2), facecolor=SURF)
    fig.subplots_adjust(left=0.06, right=0.86, bottom=0.22, top=0.85, wspace=0.62)
    for ax, L0 in zip(np.atleast_1d(axes), L0vals):
        ax.set_facecolor(SURF)
        eta = eta_of(L0); kr = kr_cs(eta); rom, rom_e = ROMAN_TABLE[L0]
        sub = sorted([r for r in rows if r["L0"] == L0], key=lambda r: r["cut"])
        vals = [rom + rom_e, rom - rom_e, kr]
        labels = []
        for key, label, color, marker in ESTIMATORS:
            pts = [(r["cut"], r[key], r["roman_mean_scatter"]) for r in sub if np.isfinite(r[key])]
            if not pts:
                continue
            xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
            vals += ys
            if key == "roman_mean":
                ax.errorbar(xs, ys, yerr=[q[2] for q in pts], color=color, lw=2, marker=marker, ms=7, capsize=3,
                            elinewidth=1, mec=SURF, mew=1.5, zorder=4)
            else:
                ax.plot(xs, ys, color=color, lw=2, marker=marker, ms=7, mec=SURF, mew=1.5, zorder=3)
            labels.append((ys[-1], label, color))
        ax.axhspan(rom - rom_e, rom + rom_e, color=GRID, alpha=0.8, lw=0, zorder=1)
        ax.axhline(rom, color=INK, lw=1.2, ls="--", zorder=2)
        ax.axhline(kr, color=KRC, lw=1.2, zorder=2)
        labels.append((rom, f"Román 2002: {rom:.2f} ± {rom_e:.2f}", INK))
        labels.append((kr, f"Kolafa–Rottner: {kr:.4f}", KRC))
        lo, hi = min(vals), max(vals)
        span = hi - lo if hi > lo else abs(hi) * 0.02 + 1e-3
        ax.set_ylim(lo - 0.06 * span, hi + 0.06 * span)
        ax.set_xscale("log"); ax.set_xticks(cuts); ax.set_xticks([], minor=True)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlim(min(cuts) / 1.3, max(cuts) * 1.3)
        ax.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_color(AXIS)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.set_xlabel("record length used, oscillations from release", color=INK2, fontsize=10)
        ax.set_ylabel("speed of sound c_s", color=INK2, fontsize=10)
        ax.set_title(f"L0 = {L0:g} σ   (η = {eta:.3f})", color=INK, fontsize=11, loc="left")
        # labels in the right margin, outside the plot, de-collided; a short connector in the
        # series colour carries identity, the text stays in ink
        y0, y1 = ax.get_ylim()
        items = sorted(labels, key=lambda z: z[0])
        fys = []
        for yv, _, _ in items:
            fy = min(max((yv - y0) / (y1 - y0), 0.0), 1.0)
            if fys and fy - fys[-1] < 0.075:
                fy = fys[-1] + 0.075
            fys.append(fy)
        if fys and fys[-1] > 1.0:
            shift = fys[-1] - 1.0
            fys = [f - shift for f in fys]
        for (yv, text, color), fy in zip(items, fys):
            fa = min(max((yv - y0) / (y1 - y0), 0.0), 1.0)
            ax.annotate(text, xy=(1.0, fa), xycoords="axes fraction", xytext=(1.05, fy), textcoords="axes fraction",
                        color=INK2, fontsize=8.5, va="center", ha="left", annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color=color, lw=1.4, shrinkA=0, shrinkB=1))
    handles = [Line2D([0], [0], color=c, lw=2, marker=m, ms=7, mec=SURF, mew=1.5, label=l) for _, l, c, m in ESTIMATORS]
    handles += [Line2D([0], [0], color=INK, lw=1.2, ls="--", label="Román 2002 Table I, grey band = its error"),
                Line2D([0], [0], color=KRC, lw=1.2, label="Kolafa–Rottner 2006")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5, labelcolor=INK2,
               bbox_to_anchor=(0.46, 0.0))
    fig.suptitle("TEST A: c_s against record length, the same trajectories cut to each length\n"
                 "N = 100, 12 Román masses, 25 seeds each, drift-first, strict health contract. "
                 "Error bar on the Román mean = 1σ scatter of per-mass c_s.",
                 color=INK, fontsize=10.5)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", dpi=200, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)


def analysis_A(out_prefix, cuts=CUTS, masses=ROMAN_MASSES, cell_dir_fn=None):
    tasks = []
    for L0s, L0 in L0S:
        for M in masses:
            cell = cell_dir_fn(L0s, L0, M) if cell_dir_fn else os.path.join(TROOT, "A_length", f"L0_{tag(L0s)}", f"m_{M}")
            runs = cell_runs(cell, M)
            if runs:
                tasks.append(((L0, M), L0, M, runs, cuts))
    results = analyse(tasks)
    rows, per_mass = [], []
    for _, L0 in L0S:
        rr, pm = summarize(results, L0, cuts, masses)
        rows += rr; per_mass += pm
    write_csv(f"{out_prefix}_cs_vs_record_length.csv", rows)
    write_csv(f"{out_prefix}_per_mass.csv", per_mass)
    figure_A(rows, cuts, f"{out_prefix}_cs_vs_record_length")
    # record length for B: shortest cut whose Roman mean changes by < 0.3 % to the next cut
    choice = {}
    for _, L0 in L0S:
        seq = [r for r in rows if r["L0"] == L0 and np.isfinite(r["roman_mean"])]
        seq.sort(key=lambda r: r["cut"])
        pick = None
        for a, b in zip(seq, seq[1:]):
            if abs(b["roman_mean"] - a["roman_mean"]) / a["roman_mean"] < 0.003:
                pick = a["cut"]; break
        choice[L0] = pick
    picks = [p for p in choice.values() if p is not None]
    # ##CHRIS 2026-09-13 11:00: if either L0 never settles, fall back to 200 oscillations, the default named in
    # the go ("use 200 if A is not yet analysed"). The earlier fallback of 1000 was an unspecified choice, and
    # interim TEST A data showed the slow-drift peak overtaking the resonance bin at 500-1000 oscillations.
    LB = max(picks) if len(picks) == len(choice) else 200
    r20 = sorted([r for r in rows if r["L0"] == 20.0 and np.isfinite(r["roman_mean"])], key=lambda r: r["cut"])
    ch = abs(r20[-1]["roman_mean"] - r20[-2]["roman_mean"]) / r20[-2]["roman_mean"] if len(r20) >= 2 else float("nan")
    summary = dict(rows=rows, per_mass=per_mass, choice_per_L0={str(k): v for k, v in choice.items()}, LB=LB,
                   longest_pair=[r20[-2]["cut"], r20[-1]["cut"]] if len(r20) >= 2 else None,
                   change_longest_pair_L0_20=ch)
    return summary


# ------------------------------------------------------------------ TEST C
def analysis_C(cell, M=500, n_traces=10, cuts=CUTS):
    runs = [q for q in cell_runs(cell, M) if not q[2]][:n_traces]
    worst = {2: 0.0, 4: 0.0, 8: 0.0}; changed = {2: 0, 4: 0, 8: 0}; per = []
    for r, path, _ in runs:
        t, x, nup = _load(path)
        dt = (t[-1] - t[0]) / (len(t) - 1)
        for N in cuts:
            n = _prefix(t, nup, N)
            if n is None:
                continue
            n8 = n - n % 8
            P, df = _spectrum(x[:n8], dt); _, nu1, _ = _peak(P, df)
            for s in (2, 4, 8):
                Ps, dfs = _spectrum(x[:n8:s], dt * s); _, nus, _ = _peak(Ps, dfs)
                rel = abs(nus - nu1) / nu1
                worst[s] = max(worst[s], rel); changed[s] += int(rel > 0)
                per.append(dict(run=r, cut=N, stride=s, nu_every_sample=nu1, nu_strided=nus, rel_change=rel))
    samples_per_period = None
    if runs:
        t, x, nup = _load(runs[0][1]); samples_per_period = len(t) / (float(t[-1]) * nup)
    return dict(n_traces=len(runs), worst=worst, changed=changed, rows=per, samples_per_period=samples_per_period,
                passed=bool(runs) and worst[8] < 0.001)


# ------------------------------------------------------------------ TEST B analysis
def analysis_B(LB, out_prefix, cell_dir_fn=None, holds=HOLDS, masses=B_MASSES):
    tasks = []
    for L0s, L0 in L0S:
        for hold in holds:
            for M in masses:
                cell = cell_dir_fn(L0s, L0, hold, M) if cell_dir_fn else os.path.join(TROOT, "B_hold", f"L0_{tag(L0s)}", f"hold_{hold}", f"m_{M}")
                runs = cell_runs(cell, M)
                if runs:
                    tasks.append(((L0, hold, M), L0, M, runs, [LB]))
    results = analyse(tasks)
    byk = {tuple(q["key"]): q["cuts"][LB] for q in results}
    rows, pm = [], []
    for _, L0 in L0S:
        for hold in holds:
            xs, ym, yd, yp = [], [], [], []
            used = disc = 0
            for M in masses:
                c = byk.get((L0, hold, M))
                if c is None:
                    continue
                xs.append(x_of(M, L0)); ym.append(c["nu_mean"]); yd.append(c["nu_median"]); yp.append(c["par_mean"])
                used += c["n_used"]; disc += c["n_disc"]
                pm.append(dict(L0=L0, hold=hold, M=M, n_used=c["n_used"], nu_mean=c["nu_mean"], nu_sd=c["nu_sd"],
                               nu_median=c["nu_median"], drift_fraction_bins_1_3=c["drift_frac"]))
            cm, sm, nm = slope(xs, ym); cd, _, _ = slope(xs, yd); cp, _, _ = slope(xs, yp)
            # paired difference against the current hold, same seeds
            dper = []
            base = holds[0]
            for M in masses:
                a, b = byk.get((L0, base, M)), byk.get((L0, hold, M))
                if a is None or b is None:
                    continue
                ua = {r: nu for r, nu, d in zip(a["r"], a["nu"], a["disc"]) if not d}
                ub = {r: nu for r, nu, d in zip(b["r"], b["nu"], b["disc"]) if not d}
                common = sorted(set(ua) & set(ub))
                if common:
                    dper.append(np.mean([ub[r] - ua[r] for r in common]) / x_of(M, L0))
            dmean = float(np.mean(dper)) if dper else float("nan")
            derr = float(np.std(dper, ddof=1) / math.sqrt(len(dper))) if len(dper) > 1 else float("nan")
            rows.append(dict(L0=L0, hold=hold, record_length=LB, roman_mean=cm, roman_mean_scatter=sm, roman_median=cd,
                             parabolic=cp, n_masses=nm, trajectories_used=used, trajectories_discarded=disc,
                             paired_delta_cs_vs_hold_2000=dmean, paired_delta_err=derr))
    write_csv(f"{out_prefix}_cs_vs_hold.csv", rows)
    write_csv(f"{out_prefix}_per_mass.csv", pm)
    return dict(rows=rows, per_mass=pm, LB=LB)


# ------------------------------------------------------------------ A1 v2 analysis and final figures
def analysis_D(table, TD, cell_fn=None, outdir=PLOTS, stem="260914_A1v2"):
    """A1 v2: c_s per density by every estimator from the same trajectories. The two headline versions are
    written in the column set plot_cs_meeting.py reads and drawn with it in mode v2 (drift-first seeding, T_i = 1,
    error bar = 1 sigma scatter of per-mass c_s):
      roman_mean_of_single_peaks  - peak of each trajectory's spectrum, averaged over seeds (Roman 2002)
      peak_of_averaged_spectrum   - spectra averaged over seeds first, then the largest bin of the average"""
    tasks = []
    for leaf in table:
        L0 = float(leaf["L0"])
        for M in A1_MASSES:
            cell = cell_fn(leaf, M) if cell_fn else os.path.join(DROOT, leaf["leaf"], f"m_{M}")
            runs = cell_runs(cell, M)
            if runs:
                tasks.append(((leaf["eta"], M), L0, M, runs, [TD]))
    results = analyse(tasks)
    byk = {tuple(q["key"]): q["cuts"][TD] for q in results if TD in q["cuts"]}
    rows = []
    for leaf in table:
        L0 = float(leaf["L0"]); eta = leaf["eta"]
        xs, ym, yd, ya, yf, drift = [], [], [], [], [], []
        used = disc = 0
        for M in A1_MASSES:
            c = byk.get((eta, M))
            if c is None:
                continue
            xs.append(x_of(M, L0)); ym.append(c["nu_mean"]); yd.append(c["nu_median"])
            ya.append(c["avg_peak"]); yf.append(c["fit_f0"])
            used += c["n_used"]; disc += c["n_disc"]; drift.append(c["drift_frac"])
        if not xs:
            continue
        cm, sm, nm = slope(xs, ym); cd, sdd, _ = slope(xs, yd); ca, sa, _ = slope(xs, ya); cf, sf, _ = slope(xs, yf)
        rows.append(dict(eta=f"{eta:.6f}", L0=leaf["L0"], roman_mean=cm, roman_scatter=sm, roman_median=cd,
                         median_scatter=sdd, avg_spectrum_peak=ca, avg_spectrum_scatter=sa, spectrum_fit=cf,
                         spectrum_fit_scatter=sf, n_masses=nm, trajectories_used=used, trajectories_discarded=disc,
                         max_drift_fraction=float(np.nanmax(drift)) if drift else float("nan"),
                         cs_KR2006=kr_cs(eta) if eta <= 0.69 else float("nan")))
    write_csv(os.path.join(outdir, f"{stem}_cs_vs_eta_all_estimators.csv"), rows)
    figs, logs = [], []
    plot_script = os.path.join(PLOTS, "plot_cs_meeting.py")
    manifest = os.path.join(PLOTS, "routeA_fit_input_manifest_20260909.csv")
    for key, err, name, label in (
            ("roman_mean", "roman_scatter", "roman_mean_of_single_peaks", "mean of single-trajectory FFT peaks (Román 2002)"),
            ("avg_spectrum_peak", "avg_spectrum_scatter", "peak_of_averaged_spectrum", "peak of the seed-averaged FFT spectrum"),
            ("roman_median", "median_scatter", "median_of_single_peaks", "median of single-trajectory FFT peaks (secondary)"),
            ("spectrum_fit", "spectrum_fit_scatter", "fit_of_averaged_spectrum", "resonance fit to the seed-averaged spectrum (secondary)")):
        csvp = os.path.join(outdir, f"{stem}_{name}_cs_vs_eta.csv")
        write_csv(csvp, [dict(eta=r["eta"], L0=r["L0"], c_s=f"{r[key]:.5f}", c_s_err="0", c_s_scatter_mass=f"{r[err]:.5f}")
                         for r in rows if np.isfinite(r[key]) and np.isfinite(r[err])])
        prefix = os.path.join(outdir, f"{stem}_{name}")
        res = subprocess.run(["/opt/homebrew/bin/python3", plot_script, csvp, manifest, prefix, "v2", label],
                             cwd=PLOTS, capture_output=True, text=True, timeout=900)
        logs.append(res.stdout[-400:] + res.stderr[-400:])
        figs += [os.path.basename(prefix) + "_cs_vs_eta.png", os.path.basename(prefix) + "_cs_idealgas_zoom.png"]
    return dict(rows=rows, TD=TD, figures=figs, plot_logs=logs)


# ------------------------------------------------------------------ ledgers, report
def ledger_stats(test):
    path = os.path.join(STATE, f"ledger_{test}.jsonl")
    latest = {}
    if os.path.exists(path):
        for ln in open(path):
            try:
                rec = json.loads(ln)
            except Exception:  # noqa: BLE001
                continue
            latest[(rec["cell"], rec["r"])] = rec
    recs = list(latest.values())
    out = dict(attempted=len(recs), ok=sum(r["ok"] for r in recs), failed=sum(not r["ok"] for r in recs),
               with_health=sum(1 for r in recs if any(r[c] for c in COUNTERS)),
               seconds=sum(r["seconds"] for r in recs))
    for c in COUNTERS:
        out[c] = sum(r[c] for r in recs)
    return out, recs


def fmt(v, nd=4):
    try:
        return "n/a" if v is None or not np.isfinite(v) else f"{v:.{nd}f}"
    except TypeError:
        return str(v)


def load_json(name):
    p = os.path.join(STATE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def write_report():
    A = load_json("A_analysis.json"); C = load_json("C_analysis.json"); B = load_json("B_analysis.json")
    D = load_json("D_decision.json")
    L = [f"# Overnight tests 2026-09-13: report", "",
         f"Last refreshed {now()}. Generated by `hspist3/validation/tests_20260913.py`; the step-by-step log is "
         f"`260913_tests_STATUS.md`. Estimator exactly as in `260913_method_roman_fft_COWORK.md`. Numbers in this "
         f"file are read from the CSV/JSON outputs listed under each test.", ""]
    # ---- A
    L += ["## TEST A: record length", ""]
    if A:
        L += ["Same 600 trajectories, each cut to the first N predicted oscillations after release. "
              "Román mean is primary; its error is the 1σ scatter of per-mass c_s.", ""]
        for _, L0 in L0S:
            eta = eta_of(L0); kr = kr_cs(eta); rom, re_ = ROMAN_TABLE[L0]
            L += [f"**L0 = {L0:g} (η = {eta:.4f}).** Román 2002 Table I {rom:.2f} ± {re_:.2f}; Kolafa–Rottner {kr:.4f}.", "",
                  "| oscillations | Román mean | ± scatter | median peak | peak of averaged spectrum | 3-bin parabolic | spectrum fit | used / discarded | mean incl. discarded |",
                  "|---|---|---|---|---|---|---|---|---|"]
            for r in sorted([q for q in A["rows"] if q["L0"] == L0], key=lambda q: q["cut"]):
                L.append(f"| {r['cut']} | {fmt(r['roman_mean'])} | {fmt(r['roman_mean_scatter'])} | {fmt(r['roman_median'])} | {fmt(r.get('avg_spectrum_peak'))} | "
                         f"{fmt(r['parabolic'])} | {fmt(r['spectrum_fit'])} | {r['trajectories_used']} / {r['trajectories_discarded']} | "
                         f"{fmt(r['roman_mean_incl_discarded'])} |")
            L += ["", f"Fraction of seeds whose peak sits in bins 1–3 (the drift peak), L0 = {L0:g}:", "",
                  "| M | " + " | ".join(str(c) for c in CUTS) + " |", "|---|" + "---|" * len(CUTS)]
            for M in ROMAN_MASSES:
                cells = []
                for N in CUTS:
                    q = next((p for p in A["per_mass"] if p["L0"] == L0 and p["M"] == M and p["cut"] == N), None)
                    cells.append(fmt(q["drift_fraction_bins_1_3"], 2) if q else "n/a")
                L.append(f"| {M} | " + " | ".join(cells) + " |")
            L.append("")
        L += [f"Record length chosen for TEST B: **{A['LB']} oscillations** (shortest length whose Román mean changes "
              f"by < 0.3 % to the next length; per L0: {A['choice_per_L0']}; the larger is used, "
              f"and 200, the default named in the go, if either L0 never settles).",
              f"Condition (i) for D: Román mean at L0 = 20 changes by {100 * A['change_longest_pair_L0_20']:.3f} % between "
              f"{A['longest_pair']} oscillations (limit 0.5 %).", "",
              "Figure: `260909_plots/260913_testA_cs_vs_record_length.png` (and .pdf). Tables: "
              "`260909_plots/260913_testA_cs_vs_record_length.csv`, `260909_plots/260913_testA_per_mass.csv`.", ""]
    else:
        L += ["Not analysed yet.", ""]
    # ---- B
    L += ["## TEST B: settling time", ""]
    if B:
        L += [f"Record length {B['LB']} oscillations. Seeds are identical at every hold, so the Δ column is a paired "
              "difference against hold = 2000 steps (mean over masses of the per-mass paired Δν/x_M, ± SD/√masses).", "",
              "| L0 | hold steps | Román mean | ± scatter | median peak | 3-bin parabolic | Δ vs 2000 | ± | used / discarded |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in B["rows"]:
            L.append(f"| {r['L0']:g} | {r['hold']} | {fmt(r['roman_mean'])} | {fmt(r['roman_mean_scatter'])} | {fmt(r['roman_median'])} | "
                     f"{fmt(r['parabolic'])} | {fmt(r['paired_delta_cs_vs_hold_2000'])} | {fmt(r['paired_delta_err'])} | "
                     f"{r['trajectories_used']} / {r['trajectories_discarded']} |")
        r75 = [r for r in B["rows"] if r["L0"] == 7.5 and r["hold"] == max(HOLDS)]
        if r75 and np.isfinite(r75[0]["paired_delta_cs_vs_hold_2000"]) and np.isfinite(r75[0]["paired_delta_err"]):
            d, e = r75[0]["paired_delta_cs_vs_hold_2000"], r75[0]["paired_delta_err"]
            moved = abs(d) > 2 * e
            L += ["", f"**Does the L0 = 7.5 value move with hold?** Hold {max(HOLDS)} against 2000: Δc_s = {d:+.4f} ± {e:.4f} "
                  f"({d / e:+.1f}σ). By a 2σ criterion that is **{'a move' if moved else 'no significant move'}**."]
        L += ["", "Tables: `260909_plots/260913_testB_cs_vs_hold.csv`, `260909_plots/260913_testB_per_mass.csv`.", ""]
    else:
        L += ["Not run or not analysed yet.", ""]
    # ---- C
    L += ["## TEST C: stride", ""]
    if C:
        L += [f"{C['n_traces']} TEST A traces at L0 = 20, M = 500, about {fmt(C['samples_per_period'], 0)} samples per oscillation "
              "as recorded. Peak recomputed from every 2nd, 4th and 8th sample at every cut length "
              "(records trimmed to a multiple of 8 samples so the frequency grid is identical).", "",
              "| stride | samples per oscillation | largest relative change in ν | peaks that changed bin |", "|---|---|---|---|"]
        for s in ("2", "4", "8"):
            L.append(f"| {s} | ~{fmt((C['samples_per_period'] or 0) / int(s), 0)} | {100 * C['worst'][s]:.4f} % | {C['changed'][s]} |")
        L += ["", f"Pass criterion < 0.1 % at stride 8: **{'PASSED' if C['passed'] else 'FAILED'}**.", ""]
    else:
        L += ["Not run yet.", ""]
    # ---- D
    L += ["## D: A1 v2 decision", ""]
    if D:
        L += [f"Decided {D['decided']}. Conditions: (i) L0 = 20 longest-pair change {100 * D['change_i']:.3f} % < 0.5 %: "
              f"**{D['cond_i']}**; (ii) TEST C passed: **{D['cond_ii']}**; (iii) free disk after the D estimate "
              f"{D['free_after_gb']:.1f} GB > 20 GB: **{D['cond_iii']}**.",
              f"Target oscillations {D['TD']}; ~{D_SAMPLES_PER_PERIOD} samples per oscillation via an explicit stride per (η, M); "
              f"{D['n_runs']} trajectories over {D['n_eta']} densities; estimated disk {D['disk_gb']:.1f} GB.",
              f"Runtime estimate from A1 step rates: {D['wall_h_a1']:.1f} h wall on {JOBS} slots; calibrated by TEST A's measured "
              f"step rate: {D['wall_h_cal']:.1f} h (used for the 10 h rule). Order: {D['order']}.",
              f"Launched: **{D['launched']}**. {D.get('why', '')}", ""]
        dst, _ = ledger_stats("D")
        if D["launched"]:
            L += [f"Progress at refresh: {dst['ok']}/{D['n_runs']} trajectories complete, {dst['failed']} failed.", ""]
        DA = load_json("D_analysis.json")
        if DA:
            L += ["**A1 v2 result.** c_s per density from the same trajectories, two headline ways: the Román mean of "
                  "single-trajectory peaks, and the peak of the seed-averaged spectrum. Error = 1σ scatter of per-mass c_s.", "",
                  "| η | Román mean | ± | peak of averaged spectrum | ± | median peak | spectrum fit | KR | max drift fraction | used / discarded |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
            for r in DA["rows"]:
                L.append(f"| {r['eta']} | {fmt(r['roman_mean'])} | {fmt(r['roman_scatter'])} | {fmt(r['avg_spectrum_peak'])} | "
                         f"{fmt(r['avg_spectrum_scatter'])} | {fmt(r['roman_median'])} | {fmt(r['spectrum_fit'])} | "
                         f"{fmt(r['cs_KR2006'])} | {fmt(r['max_drift_fraction'], 2)} | {r['trajectories_used']} / {r['trajectories_discarded']} |")
            L += ["", "Final figures: " + ", ".join(f"`260909_plots/{f}`" for f in DA["figures"]), ""]
    else:
        L += ["Not decided yet.", ""]
    # ---- health
    L += ["## Health counts", "", "| test | attempted | complete | failed | with any health event | forced_advance | clamp_repair | overlap_repair | wall_overdue | core-hours |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for t in ("A", "B", "D"):
        s, _ = ledger_stats(t)
        L.append(f"| {t} | {s['attempted']} | {s['ok']} | {s['failed']} | {s['with_health']} | {s['forced_advance']} | "
                 f"{s['clamp_repair']} | {s['overlap_repair']} | {s['wall_overdue']} | {s['seconds'] / 3600:.1f} |")
    L.append("")
    # ---- commands
    L += ["## Commands, verbatim", ""]
    lc = os.path.join(STATE, "launch_command.txt")
    if os.path.exists(lc):
        L += ["Orchestrator:", "", "```", open(lc).read().strip(), "```", ""]
    for t in ("A", "B", "D"):
        cf = os.path.join(STATE, f"command_{t}.txt")
        if os.path.exists(cf):
            L += [f"First {t} trajectory (every run differs only in --lengths, --wall-masses, --wall-hold-steps, "
                  f"--target-oscillations, --speed-sound-log-stride, --seed, --speed-sound-exact-seed and the run dir):",
                  "", "```", open(cf).read().strip(), "```", ""]
    L += ["Seeds: run r of mass index m in a cell uses speed_sound_run_seed(base, 0, m, r) with base "
          f"A = {A_BASE} + L0 index, B = {B_BASE} + L0 index (same at every hold), D = {D_BASE} + density index.", ""]
    # ---- git
    L += ["## git status --short", "", "```"]
    try:
        L.append(subprocess.run(["git", "status", "--short"], cwd=REPO, capture_output=True, text=True, timeout=120).stdout.rstrip())
    except Exception as e:  # noqa: BLE001
        L.append(f"(git status failed: {e})")
    L += ["```", "", "## git log -1", "", "```"]
    try:
        L.append(subprocess.run(["git", "log", "-1"], cwd=REPO, capture_output=True, text=True, timeout=60).stdout.rstrip())
    except Exception as e:  # noqa: BLE001
        L.append(f"(git log failed: {e})")
    L += ["```", ""]
    tmp = REPORT + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\n".join(L))
    os.replace(tmp, REPORT)


def save_json(name, obj):
    with open(os.path.join(STATE, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))


def marker(name):
    return os.path.join(STATE, name)


# ------------------------------------------------------------------ pipeline
def pipeline():
    os.makedirs(STATE, exist_ok=True)
    with open(marker("launch_command.txt"), "w") as fh:
        fh.write("caffeinate -i -s /opt/homebrew/bin/python3 " + " ".join(shlex.quote(a) for a in sys.argv) + "\n")
    ok, why = power_ok()
    status(f"orchestrator start (pid {os.getpid()}); power: {why}; free disk "
           f"{shutil.disk_usage(ROOT).free / 1e9:.0f} GB. The queue paused at 10:22 (top-up part B, alpha = 2, A4) stays paused.")
    table = a1_leaf_table()
    status(f"step rates read from {len(table)} A1 densities")

    # ---------------- TEST A (+ C beside it)
    if not os.path.exists(marker("A.done")):
        jobs = jobs_A(table)
        est = lpt([j["cost"] for j in jobs], JOBS)
        status(f"TEST A launch: {len(jobs)} trajectories (2 L0 × 12 masses × 25 seeds, 1000 oscillations, "
               f"--oscillation-safety=1.0 so the record is exactly 1000 predicted oscillations); estimate {est / 3600:.1f} h wall")

        def c_watch():
            if os.path.exists(marker("C.done")):
                return
            cell = os.path.join(TROOT, "A_length", "L0_20p0", "m_500")
            while len([q for q in cell_runs(cell, 500) if not q[2]]) < 10:
                time.sleep(60)
            try:
                res = analysis_C(cell)
                save_json("C_analysis.json", res)
                write_csv(os.path.join(PLOTS, "260913_testC_stride.csv"), res["rows"])
                open(marker("C.done"), "w").close()
                status(f"TEST C done beside A: largest relative change in ν at stride 2/4/8 = "
                       f"{100 * res['worst'][2]:.4f} / {100 * res['worst'][4]:.4f} / {100 * res['worst'][8]:.4f} % "
                       f"over {res['n_traces']} traces × {len(CUTS)} lengths; {'PASSED' if res['passed'] else 'FAILED'} (< 0.1 % at stride 8)")
            except Exception as e:  # noqa: BLE001
                status(f"TEST C FAILED with an exception: {e!r}")

        threading.Thread(target=c_watch, daemon=True).start()
        done = run_jobs("A", jobs)
        st, _ = ledger_stats("A")
        open(marker("A.done"), "w").close()
        status(f"TEST A runs done: {done}; trajectories with any health event: {st['with_health']}")
    if not os.path.exists(marker("C.done")):
        res = analysis_C(os.path.join(TROOT, "A_length", "L0_20p0", "m_500"))
        save_json("C_analysis.json", res)
        write_csv(os.path.join(PLOTS, "260913_testC_stride.csv"), res["rows"])
        open(marker("C.done"), "w").close()
        status(f"TEST C done after A: worst at stride 8 = {100 * res['worst'][8]:.4f} %; {'PASSED' if res['passed'] else 'FAILED'}")

    # ---------------- A analysis
    A = load_json("A_analysis.json")
    if A is None:
        try:
            A = analysis_A(os.path.join(PLOTS, "260913_testA"))
            save_json("A_analysis.json", A)
            parts = []
            for _, L0 in L0S:
                r = sorted([q for q in A["rows"] if q["L0"] == L0], key=lambda q: q["cut"])
                parts.append(f"L0={L0:g}: " + ", ".join(f"{q['cut']}→{fmt(q['roman_mean'], 3)}" for q in r))
            status("TEST A analysed. Román mean c_s by length: " + "; ".join(parts) +
                   f". Record length for B = {A['LB']}. Figure and tables written to 260909_plots/260913_testA_*")
        except Exception:  # noqa: BLE001
            status("TEST A analysis FAILED:\n\n```\n" + traceback.format_exc() + "```")
            A = None
    write_report()
    LB = A["LB"] if A else 200

    # ---------------- TEST B
    if not os.path.exists(marker("B.done")):
        jobs = jobs_B(table, LB)
        est = lpt([j["cost"] for j in jobs], JOBS)
        status(f"TEST B launch: {len(jobs)} trajectories (2 L0 × 3 holds × 5 masses × 25 seeds, record {LB} oscillations, "
               f"paired seeds across holds); estimate {est / 3600:.1f} h wall")
        done = run_jobs("B", jobs, heartbeat=write_report)
        st, _ = ledger_stats("B")
        open(marker("B.done"), "w").close()
        status(f"TEST B runs done: {done}; trajectories with any health event: {st['with_health']}")
    B = load_json("B_analysis.json")
    if B is None:
        try:
            B = analysis_B(LB, os.path.join(PLOTS, "260913_testB"))
            save_json("B_analysis.json", B)
            parts = [f"L0={r['L0']:g} hold {r['hold']}: {fmt(r['roman_mean'], 4)} (Δ {fmt(r['paired_delta_cs_vs_hold_2000'], 4)} ± {fmt(r['paired_delta_err'], 4)})"
                     for r in B["rows"]]
            status("TEST B analysed. " + "; ".join(parts))
        except Exception:  # noqa: BLE001
            status("TEST B analysis FAILED:\n\n```\n" + traceback.format_exc() + "```")
    write_report()

    # ---------------- D decision
    D = load_json("D_decision.json")
    if D is None:
        C = load_json("C_analysis.json")
        # ##CHRIS 2026-09-13 11:10, instruction from Chris: run A1 v2 tonight so the final c_s(eta) figure exists.
        # Record length fixed at the 200-oscillation cap. The "2 x L_B" rule rests on the Roman-mean convergence
        # criterion, which TEST A showed is confounded by the slow-drift peak. A 200-oscillation record can be cut
        # to 100 in analysis; a 100-oscillation record cannot be extended.
        TD = D_TARGET_CAP
        jobs = jobs_D(table, TD)
        wall_a1 = lpt([j["cost"] for j in jobs], JOBS)
        # calibrate with TEST A's measured step rate at both L0
        cal = []
        _, recs = ledger_stats("A")
        for L0s, L0 in L0S:
            leaf = leaf_for_L0(table, L0)
            meas = [r["seconds"] / (1000 / (nu_pred(leaf["nu50"], r["M"]) * DT_SIGMA) + HOLD_NOW) for r in recs if r["ok"] and r["L0"] == L0s]
            if meas:
                cal.append(float(np.median(meas)) / leaf["s_per_step"])
        factor = max(cal) if cal else 1.0
        wall_cal = wall_a1 * factor
        a_traces = glob.glob(os.path.join(TROOT, "A_length", "*", "*", "wall_x_positions_*.csv"))[:20]
        bpr = 117.0
        if a_traces:
            nb = sum(os.path.getsize(p) for p in a_traces)
            nr = 0
            for p in a_traces:
                with open(p, "rb") as fh:
                    nr += sum(1 for _ in fh) - 1
            bpr = nb / max(nr, 1)
        disk = sum(j["rows"] for j in jobs) * bpr
        free = shutil.disk_usage(ROOT).free
        change_i = A["change_longest_pair_L0_20"] if A else float("nan")
        cond_i = bool(A) and np.isfinite(change_i) and change_i < 0.005
        cond_ii = bool(C) and bool(C.get("passed"))
        cond_iii = (free - disk) > 20e9
        order = "all densities, longest runs first"
        if wall_cal > D_WALL_LIMIT_S:
            etas = sorted({j["eta"] for j in jobs})[:D_FIRST_ETAS]
            first = sorted([j for j in jobs if j["eta"] in etas], key=lambda j: -j["cost"])
            rest = sorted([j for j in jobs if j["eta"] not in etas], key=lambda j: -j["cost"])
            jobs = first + rest
            order = f"estimate above 10 h: the {D_FIRST_ETAS} lowest densities first (η ≤ {max(etas):.6f}), the other {len(set(j['eta'] for j in rest))} queued behind"
        else:
            jobs = sorted(jobs, key=lambda j: -j["cost"])
        # condition (i) is reported but no longer gates the launch (instruction from Chris, 2026-09-13 11:10)
        launched = cond_ii and cond_iii
        why = ("Condition (i) is advisory since 11:10 on Chris's instruction to run A1 v2 tonight"
               + ("" if cond_i else "; it failed, which is the slow-drift peak on long records, not a problem with the new data")
               + ". Record length fixed at the 200-oscillation cap. ")
        if not launched:
            why += "Not launched because " + ", ".join(
                n for n, c in (("(ii) failed", cond_ii), ("(iii) failed", cond_iii)) if not c) + "."
        D = dict(decided=now(), TD=TD, change_i=change_i, cond_i=cond_i, cond_ii=cond_ii, cond_iii=cond_iii,
                 free_after_gb=(free - disk) / 1e9, disk_gb=disk / 1e9, wall_h_a1=wall_a1 / 3600, wall_h_cal=wall_cal / 3600,
                 calibration_factor=factor, n_runs=len(jobs), n_eta=len({j['eta'] for j in jobs}), order=order,
                 launched=launched, why=why)
        save_json("D_decision.json", D)
        status(f"D estimate first: {len(jobs)} trajectories, target {TD} oscillations, ~{D_SAMPLES_PER_PERIOD} samples/oscillation, "
               f"disk {disk / 1e9:.1f} GB, wall {wall_a1 / 3600:.1f} h from A1 rates, {wall_cal / 3600:.1f} h calibrated by TEST A (×{factor:.2f}). "
               f"Conditions (i) {cond_i} (ii) {cond_ii} (iii) {cond_iii}. {'LAUNCHING D. Order: ' + order if launched else why}")
        write_report()
        with open(marker("D_jobs.json"), "w") as fh:
            json.dump([{k: v for k, v in j.items()} for j in jobs], fh)
    if D["launched"] and not os.path.exists(marker("D.done")):
        jobs = json.load(open(marker("D_jobs.json")))
        done = run_jobs("D", jobs, heartbeat=write_report)
        open(marker("D.done"), "w").close()
        st, _ = ledger_stats("D")
        status(f"D runs done: {done}; trajectories with any health event: {st['with_health']}")
    if D["launched"] and os.path.exists(marker("D.done")) and load_json("D_analysis.json") is None:
        try:
            Dres = analysis_D(table, D["TD"])
            save_json("D_analysis.json", Dres)
            status("A1 v2 analysed; final figures written: " + ", ".join(Dres["figures"]))
        except Exception:  # noqa: BLE001
            status("A1 v2 analysis FAILED:\n\n```\n" + traceback.format_exc() + "```")
    write_report()
    status("orchestrator finished; report refreshed")


def report_clock():
    """Refresh the report every 30 minutes whatever step is running, so something current is on disk by 08:00."""
    while True:
        time.sleep(1800)
        try:
            write_report()
        except Exception as e:  # noqa: BLE001
            print(f"report refresh failed: {e!r}", flush=True)


# ------------------------------------------------------------------ self-test on existing A1 traces
def selftest(outdir):
    global QUIET
    QUIET = True
    os.makedirs(outdir, exist_ok=True)
    leaf = {20.0: os.path.join(A1_DIR, "eta_0p196350"), 7.5: os.path.join(A1_DIR, "eta_0p523599")}

    def a_cell(L0s, L0, M):  # A1 leaves are flat: all masses in one dir, health in its run.log
        return leaf[L0]

    t0 = time.time()
    A = analysis_A(os.path.join(outdir, "selftest_A"), cuts=[25, 37], masses=A1_MASSES, cell_dir_fn=a_cell)
    print(f"[selftest] A-analysis on A1 traces ({time.time() - t0:.0f} s); B length would be {A['LB']}")
    for r in A["rows"]:
        print(f"   L0={r['L0']:>4g} cut={r['cut']:>3}  mean {fmt(r['roman_mean'])} ± {fmt(r['roman_mean_scatter'])}  "
              f"median {fmt(r['roman_median'])}  parabolic {fmt(r['parabolic'])}  fit {fmt(r['spectrum_fit'])} "
              f"({r['spectrum_fit_masses']} masses)  used {r['trajectories_used']}")
    C = analysis_C(leaf[20.0], M=500, cuts=[25, 37])
    print(f"[selftest] C on A1 L0=20 M=500: worst rel change stride 2/4/8 = {C['worst']}; changed {C['changed']}; "
          f"samples/period {fmt(C['samples_per_period'], 0)}; passed={C['passed']}")

    def b_cell(L0s, L0, hold, M):  # the same A1 data under every 'hold' label: paired Δ must be exactly 0
        return leaf[L0]

    B = analysis_B(37, os.path.join(outdir, "selftest_B"), cell_dir_fn=b_cell, holds=[2000, 20000], masses=[100, 200, 500, 1000, 2000])
    for r in B["rows"]:
        print(f"   B selftest L0={r['L0']:g} hold={r['hold']}: {fmt(r['roman_mean'])}  paired Δ {fmt(r['paired_delta_cs_vs_hold_2000'])} ± {fmt(r['paired_delta_err'])}")
    table = a1_leaf_table()
    print(f"[selftest] A1 step-rate table: {len(table)} densities; A estimate "
          f"{lpt([j['cost'] for j in jobs_A(table)], JOBS) / 3600:.2f} h; D (target 200) estimate "
          f"{lpt([j['cost'] for j in jobs_D(table, 200)], JOBS) / 3600:.2f} h, rows {sum(j['rows'] for j in jobs_D(table, 200)) / 1e6:.1f} M")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", metavar="OUTDIR")
    ap.add_argument("--report-only", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest(a.selftest)
    elif a.report_only:
        write_report()
    else:
        os.makedirs(STATE, exist_ok=True)
        threading.Thread(target=report_clock, daemon=True).start()
        try:
            pipeline()
        except Exception:  # noqa: BLE001
            status("ORCHESTRATOR CRASHED:\n\n```\n" + traceback.format_exc() + "```")
            try:
                write_report()
            except Exception:  # noqa: BLE001
                pass
            raise
