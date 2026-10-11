#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 28, REGISTRATION C1; 261012 sec. 4.7.57; proposed in sec. 4.7.56): C1 -- gen3's schedule
audit after every executed event on Mac starts of the registered Test G M = 50 gen3 command, each trajectory with an unaudited twin.
Committed before the first C1 run.
THE COMMAND: the B command of cluster/gen3_koa_261009/testG_worker_koa.sh (unchanged since 1cdb7ff6), READ FROM THAT FILE and filled
with the task's fields (the worker's own variables), run as the worker runs it: in a fresh run folder as working directory, with
HD_KE_TRACE=1 (in the command) and HD_CONTACT_AUDIT=1 (exported by the worker). The audited run appends --gen3-audit-every=1; its twin
is the same command without it. The binary: the Mac build of engine-gen3 a643bc7 (build_mac_a643bc7.sh; bin_a643bc7/00ALLINONE,
SHA-256 checked against build_record.txt).
THE TASKS: the Test G KOA task list (tasks_testG_koa_261009.txt, its registered SHA-256 checked), the lines of M = 50, engine gen3,
block main, in task-list order (r = 0..399). The sample is the first n of them; the TIMING trajectory, outside the sample, is the
last of them (r = 399), audited, run alone first. THE n RULE (sec. 4.7.56): n = 200, unless the timing trajectory takes more than
4.5 min (270 s), then n = floor(0.9 x 57,600 s / t_1) (8 processes for 2 h).
THE OUTPUTS: the worker's handling: gen3's psi6(t) file gets its SHA-256 in the run folder's .sha256_uncompressed, then gzip -9 (as
testG_worker_koa.sh does right after every trajectory); the trace stays as written, its SHA-256 recorded too (the twin comparison
below is information). Run folders: runs/<audit|twin>/r<r>/ (data, not committed); the per-run table c1_runs.tsv and the verdict
c1_verdict_output.txt are committed.
THE VERDICT (registered):
  PASS          in all n: every count of the run's [EDMD3-AUDIT] summary line is 0 (missing pair/wall/cross/div/piston, extra
                pair/wall/cross/div/piston, deferred earlier than eligible pair/wall/div/piston, bands missing/extra/short, duplicate
                disagreeing, cells inconsistent), max |dt| <= 1e-9 units, the event hash and the per-class event counts (ev_pair,
                ev_wall, ev_cross, ev_div, ev_piston, ev_band of [EDMD3-HEALTH]) equal to the twin's, and clean=1 in both;
  FAIL          any missing, extra or deferred-early event in any of the n (the first finding line is quoted);
  NOT DECIDED   otherwise (a time disagreement only, a twin that differs, a run that failed or lacks its lines).
usage (from hspist3/):
  python3 experiments_gen3_c1_261011/c1_audit_mac_261011.py --selftest
  python3 experiments_gen3_c1_261011/c1_audit_mac_261011.py timing              # the timing trajectory, alone; fixes n
  python3 experiments_gen3_c1_261011/c1_audit_mac_261011.py sample --jobs 6     # the n sample trajectories and their twins
  python3 experiments_gen3_c1_261011/c1_audit_mac_261011.py verdict             # the table and the verdict (reads the run folders)
"""
import concurrent.futures as cf, gzip, hashlib, math, os, re, shlex, shutil, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
TASKS = os.path.join(HS, "cluster", "gen3_koa_261009", "tasks_testG_koa_261009.txt")
TASKS_SHA = "cff8e6f9ef3d95440f54b4c621fff776b64e48700083e7258e472351fc6412b9"
WORKER = os.path.join(HS, "cluster", "gen3_koa_261009", "testG_worker_koa.sh")
BIN = os.path.join(HERE, "bin_a643bc7", "00ALLINONE")
BUILD_RECORD = os.path.join(HERE, "bin_a643bc7", "build_record.txt")
RUNS = os.path.join(HERE, "runs")
NFILE = os.path.join(HERE, "c1_n.txt")
TSV = os.path.join(HERE, "c1_runs.tsv")
N_REG, T_LIMIT, BUDGET = 200, 270.0, 57600.0
AUDIT_FLAG = "--gen3-audit-every=1"
HEALTH = re.compile(r"EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue")     # the worker's HEALTH pattern
EV = ("ev_pair", "ev_wall", "ev_cross", "ev_div", "ev_piston", "ev_band")
AUD_RE = re.compile(r"^\[EDMD3-AUDIT\] audited states (\d+); missing pair/wall/cross/div/piston (\d+)/(\d+)/(\d+)/(\d+)/(\d+); "
                    r"extra (\d+)/(\d+)/(\d+)/(\d+)/(\d+); deferred earlier than eligible pair/wall/div/piston (\d+)/(\d+)/(\d+)/(\d+); "
                    r"bands missing/extra/short (\d+)/(\d+)/(\d+); duplicate disagreeing (\d+); cells inconsistent (\d+); "
                    r"max \|dt\| (\S+); max \|dt\|/horizon \(\|dt\| > 1e-9 only\) (\S+)$", re.M)
AUD_KEYS = (["states"] + [f"missing_{c}" for c in ("pair", "wall", "cross", "div", "piston")] +
            [f"extra_{c}" for c in ("pair", "wall", "cross", "div", "piston")] +
            [f"deferred_early_{c}" for c in ("pair", "wall", "div", "piston")] +
            ["band_missing", "band_extra", "band_short", "dup_disagree", "cells_inconsistent"])
FINDING_RE = re.compile(r"^\[EDMD3-AUDIT\] (?!audited states)\S+ .*$", re.M)


# ---------------------------------------------------------------- the command, read from the worker
def worker_template(text=None):
    """the worker's B command (the binary's arguments), as written in testG_worker_koa.sh"""
    t = open(WORKER).read() if text is None else text
    m = re.search(r'HD_KE_TRACE=1 "\$HD_BIN" (.*?) > "\$tmp/stdout\.log" 2>&1', t, re.S)
    if not m: sys.exit("STOP: the B command was not found in the worker")
    if 'export HD_CONTACT_AUDIT=1' not in t: sys.exit("STOP: the worker does not export HD_CONTACT_AUDIT=1")
    return m.group(1).replace("\\\n", " ")


def argv_for(task, run_dir, binp=BIN, audited=False, text=None):
    """task = the task line's fields: B rel M r seed L0 H NS stride base eng block"""
    _, rel, M, r, seed, L0, H, NS, stride, base, eng, blk = task
    s = worker_template(text)
    for k, v in (("$((2*NS))", str(2 * int(NS))), ('"$tmp"', run_dir), ("$NS", NS), ("$H", H), ("$L0", L0), ("$M", M),
                 ("$base", base), ("$stride", stride), ("$seed", seed), ("$eng", eng)):
        s = s.replace(k, v)
    if "$" in s: sys.exit(f"STOP: an unfilled variable in the worker's command: {s}")
    return [binp] + shlex.split(s) + ([AUDIT_FLAG] if audited else [])


def tasks():
    raw = open(TASKS, "rb").read()
    if hashlib.sha256(raw).hexdigest() != TASKS_SHA: sys.exit("STOP: the Test G task list is not the registered one")
    T = [l.split() for l in raw.decode().splitlines()]
    T = [p for p in T if p[0] == "B" and p[2] == "50" and p[10] == "gen3" and p[11] == "main"]
    if [int(p[3]) for p in T] != list(range(400)): sys.exit("STOP: the M = 50 gen3 main tasks are not r = 0..399 in order")
    return T


def sha256f(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def recorded_binary_sha(text):
    """the binary's line of build_record.txt ('<sha>  00ALLINONE' at the end of a line; not the source line '00ALLINONE.c')"""
    m = re.findall(r"^\s*([0-9a-f]{64})  00ALLINONE\s*$", text, re.M)
    if len(m) != 1: sys.exit(f"STOP: {len(m)} binary SHA-256 lines in build_record.txt")
    return m[0]


def check_binary():
    want = recorded_binary_sha(open(BUILD_RECORD).read())
    if sha256f(BIN) != want: sys.exit("STOP: bin_a643bc7/00ALLINONE is not the recorded build")
    v = subprocess.run([BIN, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    if v != "00ALLINONE  git a643bc7  target mac-O3-e0pre": sys.exit(f"STOP: version line {v!r}")
    return v


# ---------------------------------------------------------------- one run, and its record
def run_one(task, mode):
    r = int(task[3]); d = os.path.join(RUNS, mode, f"r{r}")
    log = os.path.join(d, "stdout.log")
    if os.path.exists(d):
        if os.path.exists(os.path.join(d, "c1_done.txt")): return d, "done before"
        sys.exit(f"STOP: {d} exists without c1_done.txt (an interrupted run): not overwriting; nothing else started")
    os.makedirs(d)
    env = dict(os.environ, HD_CONTACT_AUDIT="1", HD_KE_TRACE="1")
    t0 = time.monotonic()
    with open(log, "w") as fo:
        rc = subprocess.run(argv_for(task, d, audited=(mode != "twin")), cwd=d, env=env, stdout=fo, stderr=subprocess.STDOUT).returncode
    wall = time.monotonic() - t0
    shas = []
    for f in sorted(os.listdir(d)):
        if f.startswith("wall_x_positions_") and f.endswith(".csv"): shas.append((sha256f(os.path.join(d, f)), f))
        if f.startswith("psi6_t_") and f.endswith(".csv"):
            p = os.path.join(d, f); shas.append((sha256f(p), f))
            with open(p, "rb") as fi, gzip.open(p + ".gz", "wb", compresslevel=9) as fo2: shutil.copyfileobj(fi, fo2)
            os.remove(p)                                      # as the worker's gzip -9: the .gz replaces the file
    with open(os.path.join(d, ".sha256_uncompressed"), "w") as fs:
        for h, f in shas: fs.write(f"{h}  {f}\n")
    with open(os.path.join(d, "c1_done.txt"), "w") as fd: fd.write(f"rc {rc} wall_s {wall:.3f} mode {mode} r {r} seed {task[4]}\n")
    return d, f"rc {rc}, {wall:.1f} s"


def record(d):
    """the run's numbers from its folder"""
    out = {"dir": d}
    done = open(os.path.join(d, "c1_done.txt")).read().split()
    out["rc"], out["wall_s"] = int(done[1]), float(done[3])
    txt = open(os.path.join(d, "stdout.log"), errors="replace").read()
    h = re.findall(r"^\[EDMD3-HEALTH\] .*$", txt, re.M)
    out["n_health"] = len(h); out["n_built"] = len(re.findall(r"^\[EDMD3\] built #", txt, re.M))
    out["other_health"] = sum(1 for l in txt.splitlines() if HEALTH.search(l) and not l.startswith("[EDMD3-HEALTH]"))
    if len(h) == 1:
        kv = dict(x.split("=", 1) for x in h[0].split() if "=" in x)
        out["clean"] = int(kv.get("clean", -1)); out["engine"] = kv.get("engine")
        for k in EV + ("ev_stale",): out[k] = int(kv[k])
        out["hash"] = kv["hash"]
    a = AUD_RE.findall(txt)
    if a:
        v = a[-1]; out.update({k: int(x) for k, x in zip(AUD_KEYS, v[:-2])}); out["max_dt"] = float(v[-2]); out["max_rel"] = float(v[-1])
        out["n_audit_lines"] = len(a)
    out["findings"] = FINDING_RE.findall(txt)
    g = re.search(r"^\[EDMD3-GAP\] contact audit over (\d+) events: max gap \[px\] pair (\S+), wall (\S+), divider (\S+),", txt, re.M)
    if g: out["gap_events"] = int(g.group(1)); out["gap_pair"], out["gap_wall"], out["gap_div"] = (float(g.group(i)) for i in (2, 3, 4))
    s = {}
    p = os.path.join(d, ".sha256_uncompressed")
    if os.path.exists(p):
        for l in open(p):
            hh, f = l.split(); s["trace" if f.startswith("wall_x") else "psi6"] = hh
    out["sha"] = s
    return out


def judge(pairs):
    """pairs: [(audited record, twin record)] -> (verdict, per-pair notes, the first finding or None)"""
    fail_cls = [k for k in AUD_KEYS if k.startswith(("missing_", "extra_", "deferred_early_"))]
    zero_cls = [k for k in AUD_KEYS if k != "states"]
    first, fail, undecided = None, False, []
    for a, t in pairs:
        ok_lines = all(k in a for k in AUD_KEYS) and "hash" in a and "hash" in t
        if not ok_lines: undecided.append((a["dir"], "a run lacks its [EDMD3-AUDIT] or [EDMD3-HEALTH] line")); continue
        if any(a[k] for k in fail_cls):
            fail = True
            if first is None: first = (a["dir"], a["findings"][0] if a["findings"] else "(no finding line printed)")
        why = []
        if any(a[k] for k in zero_cls if k not in fail_cls): why.append("a band, duplicate or cell count above 0")
        if not (a["max_dt"] <= 1e-9 and a["max_rel"] == 0.0): why.append(f"max |dt| {a['max_dt']:.3g}")
        if a["hash"] != t["hash"]: why.append("hash differs from the twin's")
        if any(a[k] != t[k] for k in EV): why.append("per-class event counts differ from the twin's")
        if not (a.get("clean") == 1 and t.get("clean") == 1): why.append("clean != 1")
        if a["rc"] or t["rc"]: why.append("exit code")
        if a["n_audit_lines"] != 1: why.append("not exactly one audit summary line")
        if why: undecided.append((a["dir"], "; ".join(why)))
    verdict = "FAIL" if fail else ("PASS" if not undecided else "NOT DECIDED")
    return verdict, undecided, first


# ---------------------------------------------------------------- the steps
def step_timing():
    T = tasks(); v = check_binary()
    if os.path.exists(NFILE): sys.exit(f"STOP: {NFILE} exists: the timing has been run")
    t = T[-1]
    print(f"C1 timing trajectory (outside the sample): M = 50 gen3 r = {t[3]} seed {t[4]}, audited; {v}", flush=True)
    d, msg = run_one(t, "timing"); rec = record(d); t1 = rec["wall_s"]
    n = N_REG if t1 <= T_LIMIT else int(math.floor(0.9 * BUDGET / t1))
    with open(NFILE, "w") as f:
        f.write(f"timing r {t[3]} seed {t[4]} wall_s {t1:.3f} audited_states {rec.get('states')} clean {rec.get('clean')}\n"
                f"n {n}  (rule: n = {N_REG} if t_1 <= {T_LIMIT:.0f} s, else floor(0.9 x {BUDGET:.0f} s / t_1))\n")
    print(open(NFILE).read(), end="")


def step_sample(jobs):
    T = tasks(); v = check_binary()
    if not os.path.exists(NFILE): sys.exit("STOP: run the timing step first (it fixes n)")
    n = int(re.search(r"^n (\d+)", open(NFILE).read(), re.M).group(1))
    if not 1 <= jobs <= 8: sys.exit("STOP: at most 8 processes")
    S = T[:n]; work = [(t, m) for t in S for m in ("audit", "twin")]
    print(f"C1 sample: n = {n} (r = 0..{n - 1}), audited and twin: {len(work)} runs on {jobs} processes; {v}", flush=True)
    t0 = time.monotonic(); k = 0
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        fut = {ex.submit(run_one, t, m): (t, m) for t, m in work}
        for f in cf.as_completed(fut):
            k += 1; d, msg = f.result()
            if k % 20 == 0 or k == len(work): print(f"  {k}/{len(work)} runs ({(time.monotonic() - t0) / 60:.1f} min)", flush=True)
    print("sample done", flush=True)


def step_verdict():
    T = tasks(); n = int(re.search(r"^n (\d+)", open(NFILE).read(), re.M).group(1))
    pairs = []
    for t in T[:n]:
        r = int(t[3]); a, b = os.path.join(RUNS, "audit", f"r{r}"), os.path.join(RUNS, "twin", f"r{r}")
        if not (os.path.exists(os.path.join(a, "c1_done.txt")) and os.path.exists(os.path.join(b, "c1_done.txt"))):
            sys.exit(f"STOP: r = {r} has no complete audited run and twin")
        pairs.append((record(a), record(b)))
    with open(TSV, "w") as f:
        cols = ["r", "seed", "rc_a", "rc_t", "wall_a_s", "wall_t_s", "clean_a", "clean_t", "hash_a", "hash_t"] + [f"{k}_a" for k in EV] + \
               [f"{k}_t" for k in EV] + AUD_KEYS + ["max_dt", "max_rel", "findings", "trace_sha_equal", "psi6_sha_equal"]
        f.write("\t".join(cols) + "\n")
        for t, (a, b) in zip(T[:n], pairs):
            row = [t[3], t[4], a["rc"], b["rc"], f"{a['wall_s']:.2f}", f"{b['wall_s']:.2f}", a.get("clean"), b.get("clean"), a.get("hash"),
                   b.get("hash")] + [a.get(k) for k in EV] + [b.get(k) for k in EV] + [a.get(k) for k in AUD_KEYS] + \
                  [a.get("max_dt"), a.get("max_rel"), len(a["findings"]), a["sha"].get("trace") == b["sha"].get("trace"),
                   a["sha"].get("psi6") == b["sha"].get("psi6")]
            f.write("\t".join(str(x) for x in row) + "\n")
    verdict, und, first = judge(pairs)
    A = [a for a, _ in pairs]; B = [b for _, b in pairs]
    print(f"# C1 (REGISTERED, decision 28; 261012 sec. 4.7.57): gen3's schedule audit after every event, M = 50, Mac starts\n")
    print(open(NFILE).read())
    print(f"binary: {check_binary()} (bin_a643bc7, SHA-256 as build_record.txt); the command: testG_worker_koa.sh's B command "
          f"(read from the file), audited run + {AUDIT_FLAG}\n")
    print("| quantity | over the n audited runs | twin runs |\n|---|---|---|")
    print(f"| trajectories | {len(A)} | {len(B)} |")
    print(f"| exit code 0 | {sum(a['rc'] == 0 for a in A)} | {sum(b['rc'] == 0 for b in B)} |")
    print(f"| clean=1 in [EDMD3-HEALTH] (exactly one record, one build line, no other health line) | "
          f"{sum(a.get('clean') == 1 and a['n_health'] == 1 and a['n_built'] == 1 and a['other_health'] == 0 for a in A)} | "
          f"{sum(b.get('clean') == 1 and b['n_health'] == 1 and b['n_built'] == 1 and b['other_health'] == 0 for b in B)} |")
    print(f"| [EDMD3-AUDIT] summary lines (one per audited run) | {sum(a.get('n_audit_lines', 0) for a in A)} | "
          f"{sum(b.get('n_audit_lines', 0) for b in B)} |")
    print(f"| audited states, total (min / median / max per run) | {sum(a.get('states', 0) for a in A)} "
          f"({min(a.get('states', 0) for a in A)} / {sorted(a.get('states', 0) for a in A)[len(A) // 2]} / "
          f"{max(a.get('states', 0) for a in A)}) | - |")
    for k in AUD_KEYS[1:]:
        print(f"| {k.replace('_', ' ')}, summed | {sum(a.get(k, 0) for a in A)} | - |")
    print(f"| max \\|dt\\| over all runs [units] | {max(a.get('max_dt', math.nan) for a in A):.3g} | - |")
    print(f"| max \\|dt\\|/horizon (\\|dt\\| > 1e-9 only), over all runs | {max(a.get('max_rel', math.nan) for a in A):.3g} | - |")
    print(f"| finding lines printed | {sum(len(a['findings']) for a in A)} | {sum(len(b['findings']) for b in B)} |")
    print(f"| event hash equal to the twin's | {sum(a.get('hash') == b.get('hash') for a, b in pairs)} of {len(pairs)} | |")
    print(f"| per-class event counts ({', '.join(EV)}) equal to the twin's | "
          f"{sum(all(a.get(k) == b.get(k) for k in EV) for a, b in pairs)} of {len(pairs)} | |")
    print(f"| executed events per run (pair+wall+cross+div), mean | "
          f"{sum(a.get('ev_pair', 0) + a.get('ev_wall', 0) + a.get('ev_cross', 0) + a.get('ev_div', 0) for a in A) / len(A):.0f} | "
          f"{sum(b.get('ev_pair', 0) + b.get('ev_wall', 0) + b.get('ev_cross', 0) + b.get('ev_div', 0) for b in B) / len(B):.0f} |")
    print(f"| wall time per run [s]: median (min-max) | {sorted(a['wall_s'] for a in A)[len(A) // 2]:.1f} "
          f"({min(a['wall_s'] for a in A):.1f}-{max(a['wall_s'] for a in A):.1f}) | {sorted(b['wall_s'] for b in B)[len(B) // 2]:.1f} "
          f"({min(b['wall_s'] for b in B):.1f}-{max(b['wall_s'] for b in B):.1f}) |")
    print(f"| information: trace SHA-256 equal to the twin's | {sum(a['sha'].get('trace') == b['sha'].get('trace') for a, b in pairs)} "
          f"of {len(pairs)} | |")
    print(f"| information: psi6(t) SHA-256 equal to the twin's | {sum(a['sha'].get('psi6') == b['sha'].get('psi6') for a, b in pairs)} "
          f"of {len(pairs)} | |")
    print(f"| information: max contact gap [px] pair / wall / divider | {max(a.get('gap_pair', 0) for a in A):.3g} / "
          f"{max(a.get('gap_wall', 0) for a in A):.3g} / {max(a.get('gap_div', 0) for a in A):.3g} | "
          f"{max(b.get('gap_pair', 0) for b in B):.3g} / {max(b.get('gap_wall', 0) for b in B):.3g} / {max(b.get('gap_div', 0) for b in B):.3g} |")
    if first: print(f"\nthe first finding: {os.path.relpath(first[0], HS)}: {first[1]}")
    if und:
        print(f"\nruns that keep the verdict from PASS ({len(und)}):")
        for d, why in und[:20]: print(f"  {os.path.relpath(d, HS)}: {why}")
    print(f"\nC1: {verdict}")


# ---------------------------------------------------------------- self-test (no simulation)
def selftest():
    ok = {}
    W = open(WORKER).read()
    t = ["B", "testG/gen3/epi8_H_H10_L10/m_50", "50", "7", "123456", "10.000000", "10.000000", "50", "26", "20261220", "gen3", "main"]
    want = ["/b", "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first", "--edmd-acc=0",
            "--particles=100", "--particles-boxes=50,50", "--height=10.000000", "--particle-radius=0.5", "--wall-thickness=0.05",
            "--wall-thickness-vis=0.05", "--lengths=10.000000", "--wall-masses=50", "--repeats=1", "--seed=20261220",
            "--wall-hold-steps=2000", "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", "--speed-sound-log-stride=26",
            "--speed-sound-run-dir=/x/r7", "--speed-sound-exact-seed=123456", "--engine=gen3"]
    ok["the worker's B command, filled from a task line, is the expected argv (hand-written from the file)"] = \
        argv_for(t, "/x/r7", binp="/b", text=W) == want and argv_for(t, "/x/r7", binp="/b", audited=True, text=W) == want + [AUDIT_FLAG]
    ok["the worker exports HD_CONTACT_AUDIT=1 and runs the binary with HD_KE_TRACE=1"] = \
        "export HD_CONTACT_AUDIT=1" in W and 'HD_KE_TRACE=1 "$HD_BIN"' in W
    br = open(BUILD_RECORD).read()
    ok["the build record's binary line is read, not its source line 00ALLINONE.c (the first timing attempt's fault)"] = \
        recorded_binary_sha(br) == sha256f(BIN) and re.search(r"([0-9a-f]{64})  00ALLINONE\.c", br).group(1) != recorded_binary_sha(br)
    ok["the n rule: 200 at t_1 = 200 s and 270 s; floor(0.9 x 57600 / 300) = 172 at 300 s"] = \
        [N_REG if x <= T_LIMIT else int(math.floor(0.9 * BUDGET / x)) for x in (200.0, 270.0, 300.0)] == [200, 200, 172]
    line = ("[EDMD3-AUDIT] audited states 1131871; missing pair/wall/cross/div/piston 0/0/0/0/0; extra 0/0/0/0/0; deferred earlier than "
            "eligible pair/wall/div/piston 0/0/0/0; bands missing/extra/short 0/0/0; duplicate disagreeing 0; cells inconsistent 0; "
            "max |dt| 4e-11; max |dt|/horizon (|dt| > 1e-9 only) 0")
    v = AUD_RE.findall(line)
    ok["the [EDMD3-AUDIT] summary line is parsed (22 numbers)"] = len(v) == 1 and len(v[0]) == 22 and v[0][0] == "1131871" and \
        float(v[0][-2]) == 4e-11
    import tempfile
    def fake(dirp, aud, hash_, ev, clean=1, findings=""):
        os.makedirs(dirp)
        open(os.path.join(dirp, "c1_done.txt"), "w").write("rc 0 wall_s 1.0 mode x r 0 seed 1\n")
        h = (f"[EDMD3-HEALTH] L0=10.0 M=50 run=0 seed=1: clean={clean} ev_pair={ev} ev_wall=2 ev_cross=3 ev_div=4 ev_piston=0 ev_band=0 "
             f"ev_stale=5 hash={hash_} engine=gen3 build=\"x\"\n")
        open(os.path.join(dirp, "stdout.log"), "w").write("[EDMD3] built #1 at t = 0\n" + findings + (aud + "\n" if aud else "") + h)
        return record(dirp)
    with tempfile.TemporaryDirectory() as td:
        bad = line.replace("missing pair/wall/cross/div/piston 0/0/0/0/0", "missing pair/wall/cross/div/piston 1/0/0/0/0")
        fl = "[EDMD3-AUDIT] missing PAIR a=3 b=9 now=12 t_bruteforce=13 t_heap=nan\n"
        a1, t1 = fake(os.path.join(td, "a1"), line, "abc", 10), fake(os.path.join(td, "t1"), "", "abc", 10)
        a2, t2 = fake(os.path.join(td, "a2"), bad, "abc", 10, findings=fl), fake(os.path.join(td, "t2"), "", "abc", 10)
        a3, t3 = fake(os.path.join(td, "a3"), line, "abc", 10), fake(os.path.join(td, "t3"), "", "abd", 10)
        a4, t4 = fake(os.path.join(td, "a4"), line.replace("max |dt| 4e-11", "max |dt| 2e-09").replace("only) 0", "only) 3e-12"), "abc", 10), \
            fake(os.path.join(td, "t4"), "", "abc", 10)
        ok["verdict: all clean -> PASS"] = judge([(a1, t1)])[0] == "PASS"
        v2 = judge([(a1, t1), (a2, t2)])
        ok["verdict: one missing event -> FAIL, and its finding line is quoted"] = v2[0] == "FAIL" and v2[2][1] == fl.strip()
        ok["verdict: a twin with another hash -> NOT DECIDED"] = judge([(a1, t1), (a3, t3)])[0] == "NOT DECIDED"
        ok["verdict: max |dt| above 1e-9 only -> NOT DECIDED"] = judge([(a4, t4)])[0] == "NOT DECIDED"
    for k, x in ok.items(): print(f"| {k} | {'PASS' if x else '**FAIL**'} |")
    print(f"\nSELF-TEST: {'PASS' if all(ok.values()) else 'FAIL'} ({sum(ok.values())} of {len(ok)})")
    return all(ok.values())


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--selftest" in a: sys.exit(0 if selftest() else 1)
    if a[:1] == ["timing"]: step_timing()
    elif a[:1] == ["sample"]: step_sample(int(a[a.index("--jobs") + 1]) if "--jobs" in a else 6)
    elif a[:1] == ["verdict"]: step_verdict()
    else: sys.exit(__doc__)
