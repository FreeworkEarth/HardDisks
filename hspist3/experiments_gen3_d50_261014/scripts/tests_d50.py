#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 50 item 3; 261012 sec. 4.7.113): THE MAC TESTS OF OPTION (ii), Part A's held path on gen3.

  T1    the Part B path unchanged: registered stage 1 task lines (cluster/gen3_stage1_261012/tasks_stage1_koa_261012.txt), each run
        with the worker's command (s1_worker_koa.sh, unchanged) under the reference build OLD and the new build NEW; equal when the
        run record's event hash, t_end and every output file's SHA-256 are equal.
  T1AF  the legacy held path unchanged (no --gen3-exact-box): the worker's AF command (conf_worker.sh) under gen3 at N 400, eta pi/8,
        L0 20, the divider at L0 + sigma_x/2 (Report D49-1's probe Q2), with the event log; OLD and NEW equal as in T1, the log included.
  T2    the exact geometry: --gen3-exact-box on the AF command at all 16 N 400 cells x the 5 stencil positions L0 + j dL (j = -2..2,
        dL = sigma_x / 2 by the registered formulas), and at eta pi/8, j = 0, for N 100, 900 and 1600. Passes when the run exits 0
        with one clean run record and its [EDMD3-HELD-GEOMETRY] line gives |box - 2 L0| <= 1e-9, |x_div - requested| <= 1e-9 and
        N_L = N_R = N_s.
  C2    gate C2, the held start against the free start (step b): a registered Part B line (shortened: the start does not depend on
        the hold or the record) and the held path at j = 0 with --seed = the line's exact seed, both with --gen3-dump-initial; equal
        when the two dumps (the box, the divider, every disk's x, y, vx, vy in %a) are byte-identical and the run headers' lattices
        agree. At pi/8 for N 100, 400, 900, 1600, and at eta 0.05 and 0.69 for N 400.
  R0    the readout changes nothing (step c): held lattice runs with and without --gen3-held-readout=K; equal when the event hash,
        t_end and every output file's SHA-256 are equal.
  R1    the readout against the event log (step c): held lattice runs with --gen3-held-readout=K and HD_PISTON_EVENTS; F_L, F_R, n_L,
        n_R of the [EDMD3-HELD] window line against reduce_AF.py's formulas on the log's D0 rows in [t0, t1) (t0, t1 the readout's;
        the log prints 12 significant digits, so the forces must agree to 1e-9 and the counts exactly);
        also u_wall_max and W_div from the log, v_div, dx_div and dW_div from the readout.
Runs go into experiments_gen3_d50_261014/tests/<label>/<run>/ (a new folder each; an existing one is refused). At most 6 processes (Chris's two viewers run beside them).
usage (from hspist3/): python3 tests_d50.py <label> <T1|T1AF|T2|C2|R0|R1> <NEW binary> [OLD binary] [--seeder=legacy|lattice]
"""
import concurrent.futures as cf, hashlib, math, os, re, subprocess, sys, time, importlib.util

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks_gen3/hspist3"
OUTROOT = os.path.join(HS, "experiments_gen3_d50_261014", "tests")
TASKS = os.path.join(HS, "cluster", "gen3_stage1_261012", "tasks_stage1_koa_261012.txt")
ENV = dict(os.environ, SDL_VIDEODRIVER="dummy", HD_CONTACT_AUDIT="1")
HEALTH = re.compile(r"EDMD-HEALTH|forced_advance|clamp_repair|overlap_repair|wall_overdue")
spec = importlib.util.spec_from_file_location("prereg", os.path.join(HS, "validation", "paper1_confinement_prereg_20261012.py"))
P = importlib.util.module_from_spec(spec); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS); spec.loader.exec_module(P)

# T1: the picks, by (cell folder, mass, replicate): dilute default, the snapshot subset, block S (orientation 0, turned), block D, N 400
T1_PICKS = [("stage1/gen3/B_N0100_eta0.050000/m_50", "0"), ("stage1/gen3/B_N0100_eta0.050000/m_500", "0"),
            ("stage1/gen3/S_N0100_eta0.600000_start0/m_100", "0"), ("stage1/gen3/S_N0100_eta0.690000_turned/m_1000", "1"),
            ("stage1/gen3/D_N0100_eta0.020000/m_50", "0"), ("stage1/gen3/B_N0400_eta0.392699/m_200", "0")]
ETAS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, math.pi / 8, 0.45, 0.50, 0.55, 0.60, 0.63, 0.66, 0.68, 0.69]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def health(log):
    """the worker's gen3 rule: one build line, one run record with clean=1, no other HEALTH line; returns (ok, hash, t_end)"""
    L = open(log, errors="replace").read().splitlines()
    other = [l for l in L if HEALTH.search(l) and not l.startswith("[EDMD3-HEALTH]")]
    nb = sum(l.startswith("[EDMD3] built #") for l in L); rec = [l for l in L if l.startswith("[EDMD3-HEALTH]")]
    ok = not other and nb == 1 and len(rec) == 1 and " clean=1 " in rec[0] and " engine=gen3 " in rec[0]
    hv = re.search(r" hash=([0-9a-f]+)", rec[0]).group(1) if rec else None
    te = re.search(r" t_end=(\S+)", rec[0]).group(1) if rec else None
    return ok, hv, te


def run(d, argv, env, cwd=None):
    os.makedirs(d)   # refuses an existing folder
    open(os.path.join(d, "cmd.txt"), "w").write(" ".join(argv) + "\n")
    t = time.time()
    with open(os.path.join(d, "stdout.log"), "w") as o:
        rc = subprocess.run(argv, stdout=o, stderr=subprocess.STDOUT, env=env, cwd=cwd or d).returncode
    return rc, time.time() - t


def b_cmd(binary, d, f):
    rel, M, r, seed, L0, H, NS, stride, base, hold, dpsi, tsnap, start, blk = f
    extra = ([] if tsnap == "0" else [f"--gen3-snapshots={tsnap}"]) + ([] if start == "default" else [f"--gen3-lattice-orientation={start}"])
    return [binary, "--mode=edmd", "--experiment=speed_of_sound", "--headless", "--kbt1", "--seed-drift-order=drift-first",
            "--edmd-acc=0", f"--particles={2 * int(NS)}", f"--particles-boxes={NS},{NS}", f"--height={H}", "--particle-radius=0.5",
            "--wall-thickness=0.05", "--wall-thickness-vis=0.05", f"--lengths={L0}", f"--wall-masses={M}", "--repeats=1", f"--seed={base}",
            f"--wall-hold-steps={hold}", "--fixed-dt=0.4", "--target-oscillations=200", "--oscillation-safety=1.0",
            "--oscillation-min-steps=10000", "--oscillation-max-steps=400000000", f"--speed-sound-log-stride={stride}",
            f"--speed-sound-run-dir={d}", f"--speed-sound-exact-seed={seed}", "--engine=gen3",
            "--gen3-seeding=lattice", "--gen3-exact-box", "--gen3-psi6-compartments", f"--psi6-every={dpsi}", "--gen3-virial-blocks=50"] + extra


def af_cmd(binary, d, N, H, L0, xw, seed, hold, post, every, exact, seeder):
    NS = N // 2
    a = [binary, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0", "--seed-drift-order=drift-first",
         f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv", f"--trace-every={every}", f"--particles={N}",
         f"--particles-boxes={NS},{NS}", "--particle-radius=0.5", f"--l0={L0:.17g}", f"--height={H:.17g}", "--num-walls=1",
         f"--wall-positions={xw:.17g}", "--wall-mass-factors=1000000000", "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
         "--eff-output=wall-ke", f"--wall-hold-steps={hold}", f"--steps={post}", "--fixed-dt=0.4", "--kbt1", f"--seed={seed}", "--engine=gen3"]
    if exact: a.append("--gen3-exact-box")
    if seeder == "lattice": a.append("--gen3-seeding=lattice")
    return a


def outputs(d):
    """every output file but the run's own command and log, the summary (a time stamp) and the command files (the binary's path)"""
    return {f: sha(os.path.join(d, f)) for f in sorted(os.listdir(d)) if f not in ("stdout.log", "cmd.txt", "summary.csv")
            and not f.endswith(".failures.csv") and not f.endswith(".md") and os.path.isfile(os.path.join(d, f))}


def pair(label, name, mk, NEW, OLD):
    res = {}
    for tag, b in (("old", OLD), ("new", NEW)):
        d = os.path.join(OUTROOT, label, f"{name}_{tag}")
        argv, env, cwd = mk(b, d)
        rc, dt = run(d, argv, env, cwd)
        ok, hv, te = health(os.path.join(d, "stdout.log"))
        res[tag] = dict(rc=rc, ok=ok, hash=hv, t_end=te, dt=dt, files=outputs(d))
    o, n = res["old"], res["new"]
    same_files = o["files"] == n["files"]
    verdict = "EQUAL" if (o["rc"] == n["rc"] == 0 and o["ok"] and n["ok"] and o["hash"] == n["hash"] and o["t_end"] == n["t_end"]
                          and same_files) else "DIFFERENT"
    return (f"| {name} | {o['rc']} / {n['rc']} | {int(o['ok'])} / {int(n['ok'])} | {o['hash']} / {n['hash']} | {o['t_end']} / {n['t_end']} "
            f"| {len(o['files'])} / {len(n['files'])}, {'all equal' if same_files else 'DIFFER: ' + ', '.join(sorted(set(o['files'].items()) ^ set(n['files'].items()), key=str)[:4].__repr__())} "
            f"| {o['dt']:.0f} / {n['dt']:.0f} | {verdict} |")


def t1(label, NEW, OLD):
    lines = [l.split() for l in open(TASKS)]
    picks = []
    for rel, r in T1_PICKS:
        f = next(x for x in lines if x[1] == rel and x[3] == r)
        picks.append((rel.split("/")[-2] + "_" + rel.split("/")[-1] + f"_r{r}", f[1:]))
    hdr = ("| run (cell, mass, replicate) | exit old / new | clean old / new | event hash old / new | t_end old / new | output files, SHA-256 "
           "| seconds old / new | verdict |\n|---|---|---|---|---|---|---|---|")
    mk = lambda f: (lambda b, d: (b_cmd(b, d, f), dict(ENV, HD_KE_TRACE="1"), d))
    with cf.ThreadPoolExecutor(4) as ex:   # 4 pairs at a time, each pair sequential: at most 4 processes
        rows = list(ex.map(lambda p: pair(label, p[0], mk(p[1]), NEW, OLD), picks))
    return hdr + "\n" + "\n".join(rows)


def stencil(N, eta):
    NS, H = N // 2, 10.0 * math.sqrt(N / 100.0)
    L0 = NS * math.pi * P.R * P.R / (H * eta)
    e = P.eta_of(NS, H, L0); sx = math.sqrt(P.leff(L0) ** 2 / (2 * NS * P.cs2(e)))
    return NS, H, L0, sx / 2.0


def t1af(label, NEW, OLD):
    NS, H, L0, dL = stencil(400, math.pi / 8)
    xw = L0 + dL
    mk = lambda b, d: (af_cmd(b, d, 400, H, L0, xw, 7740102, 600, 60, 60, False, "legacy"), dict(ENV, HD_PISTON_EVENTS=os.path.join(d, "ev.csv")), d)
    hdr = ("| run | exit old / new | clean old / new | event hash old / new | t_end old / new | output files (trace, event log), SHA-256 "
           "| seconds old / new | verdict |\n|---|---|---|---|---|---|---|---|")
    return hdr + "\n" + pair(label, f"AF_N400_pi8_L0_{L0:g}_xw_L0+dL", mk, NEW, OLD)


def geo(log):
    for l in open(log, errors="replace"):
        if l.startswith("[EDMD3-HELD-GEOMETRY]"):
            g = lambda k: re.search(k + r"=(\S+)", l).group(1)
            m = re.search(r"\|dW\| (\S+)\) x_div=(\S+) sigma \(requested (\S+), \|dx\| (\S+)\)", l)
            return dict(box=float(re.search(r"box=(\S+) x", l).group(1)), dW=float(m.group(1)), x=float(m.group(2)), dx=float(m.group(4)),
                        NL=int(g("N_L")), NR=int(g("N_R")), LeL=float(g("L_eff_L")), LeR=float(g("L_eff_R")), eta=float(g("eta")))
    return None


def t2(label, NEW, seeder):
    jobs = []
    for i, eta in enumerate(ETAS):
        NS, H, L0, dL = stencil(400, eta)
        for j in (-2, -1, 0, 1, 2): jobs.append((400, eta, j, H, L0, L0 + j * dL, 7750000 + 10 * i + (j + 2)))
    for N in (100, 900, 1600):
        NS, H, L0, dL = stencil(N, math.pi / 8)
        jobs.append((N, math.pi / 8, 0, H, L0, L0, 7760000 + N))

    def one(job):
        N, eta, j, H, L0, xw, seed = job
        d = os.path.join(OUTROOT, label, f"N{N:04d}_eta{eta:.6f}_j{j:+d}")
        rc, dt = run(d, af_cmd(NEW, d, N, H, L0, xw, seed, 120, 12, 12, True, seeder), ENV)
        ok, hv, te = health(os.path.join(d, "stdout.log")); g = geo(os.path.join(d, "stdout.log"))
        stop = next((l.strip() for l in open(os.path.join(d, "stdout.log"), errors="replace") if "ABORTING" in l or "STOP" in l), "")
        p = rc == 0 and ok and g and g["dW"] <= 1e-9 and g["dx"] <= 1e-9 and g["NL"] == g["NR"] == N // 2
        if not g:
            return f"| {N} | {eta:.6f} | {j:+d} | {L0:.12f} | {xw:.12f} | - | - | - | - | {rc} | {int(ok)} | {dt:.1f} | FAIL {stop[:120]} |"
        return (f"| {N} | {eta:.6f} | {j:+d} | {L0:.12f} | {xw:.12f} | {g['dW']:.1e} | {g['dx']:.1e} | {g['NL']} / {g['NR']} | "
                f"{g['LeL']:.9f} / {g['LeR']:.9f} | {rc} | {int(ok)} | {dt:.1f} | {'PASS' if p else 'FAIL ' + stop[:120]} |")
    with cf.ThreadPoolExecutor(6) as ex: rows = list(ex.map(one, jobs))   # Chris's two viewers run: 6 + 2 = 8 processes
    hdr = ("| N | eta | j | L0 [sigma] | x_div requested [sigma] | abs(box - 2 L0) | abs(x_div - requested) | N_L / N_R | L_eff left / right [sigma] "
           "| exit | clean | seconds | verdict |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    return hdr + "\n" + "\n".join(rows) + f"\n\n{sum('PASS' in r for r in rows)} of {len(rows)} PASS"


C2_CELLS = ["B_N0100_eta0.392699", "B_N0400_eta0.392699", "B_N0900_eta0.392699", "B_N1600_eta0.392699",
            "B_N0400_eta0.050000", "B_N0400_eta0.690000"]


def lattice_of(log):
    """the run header's lattice description (left=[...] right=[...]), or None"""
    for l in open(log, errors="replace"):
        if l.startswith("[EDMD3-RUN]"):
            m = re.search(r" (left=\[.*?\]) (right=\[.*?\])", l)
            return (m.group(1) + " " + m.group(2)) if m else None
    return None


def c2(label, NEW):
    lines = [l.split() for l in open(TASKS)]
    rows = []

    def one(cell):
        f = next(x for x in lines if x[1].split("/")[-2] == cell and x[3] == "0")
        rel, M, r, seed, L0, H, NS = f[1:8]
        df = os.path.join(OUTROOT, label, cell + "_free"); dh = os.path.join(OUTROOT, label, cell + "_held")
        fb = list(f[1:]); fb[9] = "60"                                     # the hold: 60 steps (the start comes before it)
        argv = b_cmd(NEW, df, fb)
        argv = [a for a in argv if not a.startswith(("--target-oscillations=", "--oscillation-min-steps=", "--oscillation-max-steps="))]
        argv += ["--target-oscillations=1", "--oscillation-min-steps=60", "--oscillation-max-steps=120", f"--gen3-dump-initial={df}/initial.dump"]
        rcf, _ = run(df, argv, dict(ENV, HD_KE_TRACE="1"), df)
        N = 2 * int(NS)
        ah = af_cmd(NEW, dh, N, float(H), float(L0), float(L0), int(seed), 60, 6, 6, True, "lattice") + [f"--gen3-dump-initial={dh}/initial.dump"]
        rch, _ = run(dh, ah, ENV)
        a, b = os.path.join(df, "initial.dump"), os.path.join(dh, "initial.dump")
        same = os.path.exists(a) and os.path.exists(b) and open(a, "rb").read() == open(b, "rb").read()
        la, lb = lattice_of(os.path.join(df, "stdout.log")), lattice_of(os.path.join(dh, "stdout.log"))
        okf, _, _ = health(os.path.join(df, "stdout.log")); okh, _, _ = health(os.path.join(dh, "stdout.log"))
        n = sum(1 for _ in open(a)) - 6 if os.path.exists(a) else -1
        return (f"| {cell} | {seed} | {L0} | {N} | {rcf} / {rch} | {int(okf)} / {int(okh)} | {n} | "
                f"{sha(a)[:16] if os.path.exists(a) else '-'} / {sha(b)[:16] if os.path.exists(b) else '-'} | "
                f"{'equal' if la and la == lb else 'DIFFERENT'} | {'EQUAL' if same and la == lb and rcf == rch == 0 and okf and okh else 'DIFFERENT'} |")
    with cf.ThreadPoolExecutor(3) as ex: rows = list(ex.map(one, C2_CELLS))   # each cell runs its two starts in turn: at most 3 processes
    hdr = ("| cell | exact seed | L0 [sigma] | N | exit free / held | clean free / held | disks in the dump | dump SHA-256 free / held "
           "| run header's lattice | verdict |\n|---|---|---|---|---|---|---|---|---|---|")
    return hdr + "\n" + "\n".join(rows) + f"\n\n{sum('| EQUAL |' in r for r in rows)} of {len(rows)} EQUAL"


R_CELLS = [(400, 0.05, 1), (400, math.pi / 8, 1), (400, 0.69, 1), (100, math.pi / 8, -2), (1600, 0.60, 2)]
R_HOLD, R_K, R_POST = 3000, 600, 60       # 50 sigma-time held, T_eq = 10 sigma-time, a 1 sigma-time tail


def held_window(log):
    for l in open(log, errors="replace"):
        if l.startswith("[EDMD3-HELD] window"):
            return {k: float(v) for k, v in re.findall(r" (\w+)=(\S+)", l)}
    return None


def held_t(log, which):
    for l in open(log, errors="replace"):
        if l.startswith(f"[EDMD3-HELD] {which} "):
            return {k: v for k, v in re.findall(r" (\w+)=(\S+)", l.split(" exact ")[0])}
    return None


def r_jobs():
    out = []
    for i, (N, eta, j) in enumerate(R_CELLS):
        NS, H, L0, dL = stencil(N, eta)
        out.append((N, eta, j, H, L0, L0 + j * dL, 7770000 + 10 * i))
    return out


def r0(label, NEW):
    rows = []
    for N, eta, j, H, L0, xw, seed in r_jobs():
        res = {}
        for tag, extra in (("off", []), ("on", [f"--gen3-held-readout={R_K}"])):
            d = os.path.join(OUTROOT, label, f"N{N:04d}_eta{eta:.6f}_j{j:+d}_{tag}")
            rc, dt = run(d, af_cmd(NEW, d, N, H, L0, xw, seed, R_HOLD, R_POST, 6, True, "lattice") + extra, ENV)
            ok, hv, te = health(os.path.join(d, "stdout.log"))
            res[tag] = dict(rc=rc, ok=ok, hash=hv, t_end=te, files=outputs(d), win=held_window(os.path.join(d, "stdout.log")))
        a, b = res["off"], res["on"]
        eq = a["rc"] == b["rc"] == 0 and a["ok"] and b["ok"] and a["hash"] == b["hash"] and a["t_end"] == b["t_end"] and a["files"] == b["files"]
        rows.append(f"| {N} | {eta:.6f} | {j:+d} | {a['rc']} / {b['rc']} | {int(a['ok'])} / {int(b['ok'])} | {a['hash']} / {b['hash']} | "
                    f"{a['t_end']} / {b['t_end']} | {len(a['files'])} / {len(b['files'])}, {'all equal' if a['files'] == b['files'] else 'DIFFER'} | "
                    f"{'none' if a['win'] is None else 'PRESENT'} / {'present' if b['win'] else 'MISSING'} | {'EQUAL' if eq and b['win'] else 'DIFFERENT'} |")
    hdr = ("| N | eta | j | exit off / on | clean off / on | event hash off / on | t_end off / on | output files, SHA-256 | window line off / on | verdict |"
           "\n|---|---|---|---|---|---|---|---|---|---|")
    return hdr + "\n" + "\n".join(rows) + f"\n\n{sum('| EQUAL |' in r for r in rows)} of {len(rows)} EQUAL"


def r1(label, NEW):
    import numpy as np, pandas as pd
    rows = []
    for N, eta, j, H, L0, xw, seed in r_jobs():
        d = os.path.join(OUTROOT, label, f"N{N:04d}_eta{eta:.6f}_j{j:+d}")
        rc, dt = run(d, af_cmd(NEW, d, N, H, L0, xw, seed, R_HOLD, R_POST, 6, True, "lattice") + [f"--gen3-held-readout={R_K}"],
                     dict(ENV, HD_PISTON_EVENTS=os.path.join(d, "ev.csv")))
        ok, hv, te = health(os.path.join(d, "stdout.log"))
        w = held_window(os.path.join(d, "stdout.log")); a, b = held_t(os.path.join(d, "stdout.log"), "t0"), held_t(os.path.join(d, "stdout.log"), "t1")
        t0, t1 = float(a["t_sigma"]), float(b["t_sigma"])
        e = pd.read_csv(os.path.join(d, "ev.csv"), usecols=["t_sigma", "kind", "u_wall", "dE", "dp"])
        dd = e[(e["kind"] == "D0") & (e["t_sigma"] >= t0) & (e["t_sigma"] < t1)]; dp = dd["dp"].to_numpy(); Tw = t1 - t0
        FL, FR = -dp[dp < 0].sum() / Tw, dp[dp > 0].sum() / Tw
        nL, nR = int((dp < 0).sum()), int((dp > 0).sum())
        uw, Wd = float(np.abs(dd["u_wall"]).max()) if len(dd) else 0.0, float(dd["dE"].sum())
        rl, rr = abs(w["F_L"] / FL - 1), abs(w["F_R"] / FR - 1)
        good = (rc == 0 and ok and int(w["n_L"]) == nL and int(w["n_R"]) == nR and rl < 1e-9 and rr < 1e-9 and uw == 0.0 and Wd == 0.0
                and w["v_div_t0"] == 0.0 and w["v_div_t1"] == 0.0 and w["dx_div"] == 0.0 and w["dW_div"] == 0.0 and abs(w["T_w"] / Tw - 1) < 1e-15)
        rows.append(f"| {N} | {eta:.6f} | {j:+d} | {t0:.6f}, {t1:.6f} | {w['F_L']:.12f} / {FL:.12f} | {rl:.1e} | {w['F_R']:.12f} / {FR:.12f} | {rr:.1e} | "
                    f"{int(w['n_L'])} / {nL} | {int(w['n_R'])} / {nR} | {uw:g}, {Wd:g} | {w['v_div_t0']:g}, {w['v_div_t1']:g}, {w['dx_div']:g}, {w['dW_div']:g} | "
                    f"{rc} | {int(ok)} | {'PASS' if good else 'FAIL'} |")
    hdr = ("| N | eta | j | window (t0, t1) [sigma-time] | F_L readout / log | rel. diff | F_R readout / log | rel. diff | n_L readout / log | "
           "n_R readout / log | log: u_wall_max, W_div | readout: v_div t0, t1, dx_div, dW_div | exit | clean | verdict |\n"
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    return hdr + "\n" + "\n".join(rows) + f"\n\n{sum('| PASS |' in r for r in rows)} of {len(rows)} PASS"


def main():
    label, test, NEW = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    OLD = os.path.abspath(sys.argv[4]) if len(sys.argv) > 4 and not sys.argv[4].startswith("--") else None
    seeder = next((a.split("=", 1)[1] for a in sys.argv[4:] if a.startswith("--seeder=")), "legacy")
    d = os.path.join(OUTROOT, label)
    if os.path.exists(d): sys.exit(f"REFUSED: {d} exists")
    os.makedirs(d)
    head = (f"# {test} ({label}), {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\nNEW {os.path.relpath(NEW, HS)} SHA-256 {sha(NEW)}\n"
            + (f"OLD {os.path.relpath(OLD, HS)} SHA-256 {sha(OLD)}\n" if OLD else "") + (f"seeder {seeder}\n" if test == "T2" else "") + "\n")
    body = (t1(label, NEW, OLD) if test == "T1" else t1af(label, NEW, OLD) if test == "T1AF" else c2(label, NEW) if test == "C2"
            else r0(label, NEW) if test == "R0" else r1(label, NEW) if test == "R1" else t2(label, NEW, seeder))
    out = head + body + "\n"
    open(os.path.join(d, "RESULTS.md"), "w").write(out)
    print(out)


if __name__ == "__main__":
    main()
