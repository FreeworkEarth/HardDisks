#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7, the plan author's CC task of 2026-10-09): every number of the GENERATION-3 ENGINE DESIGN
NOTE. Notes and numbers only -- no engine code. Sections:
  1  baseline: collisions per second per core on KOA for the two engines (279282b = the legacy path of 7b08827, byte-identical:
     261012 sec. 4.4.8, E0) at pi/8 (held and free divider) and near eta 0.70 (dense), N = 100 and 400.
     The KOA wall times are the profile job 15008378 (times.tsv; its legacy column reproduces the 279282b job 14983181 of sec. 4.3).
     The profile logs count divider and outer-wall events but NOT disk-disk collisions, so the collision counts come from the
     same six commands re-run on the Mac with a build of the same sources (git 7b08827, -O3 -ffp-contract=off) and the read-only
     contact audit (HD_CONTACT_AUDIT=1: '[EDMD-CONTACT] executed events N' counts every executed event) plus the event log
     (HD_PISTON_EVENTS: D and W rows, as the profile script counts them). A different binary gives a different trajectory of
     the same cell; event RATES are statistical properties -- the Mac divider/wall counts are printed against KOA's as a check.
     Then the collision-rate model: pair collisions per disk per sigma-time nu = 4 (Z - 1)/sqrt(pi) (kT = m = sigma = 1; the
     virial theorem with the flux-weighted normal speed sqrt(pi kT/m) per collision), checked against the measured counts;
     outer-wall and divider events from the wall theorem, n Z sqrt(kT/(2 pi m)) per unit length of the centre-accessible
     boundary. Z: KR (the module) for the fluid, Engel's plateau P* = 9.17 at 0.70, and the Alder-Hoover-Young high-density form
     Z = 2/a + 1.90 + 0.67 a, a = eta_cp/eta - 1, for the solid (recalled, not checked against a PDF: OPEN).
     The cost model: the time of a 2e4 sigma-time trajectory at N = 100/400/900/1600 (design geometry H = 10 sqrt(N/100), L0 from
     eta) for (a) 279282b (two-term model fitted to the measured legacy times: alpha N per non-divider event + beta N^2 per divider
     event), (b) 7b08827 (cost per event ~ N^q with q measured), (c) a constant-work engine at 2e5 and 5e5 collisions per second.
  3  exactness: double-precision ulp at the engine's absolute times (internal time unit = sigma-time / 24, positions in px, 24 px per
     sigma); the cancellation of the current quadratic root against the stable form, by exact rational arithmetic; every bare
     tolerance literal of edmd.c at 7b08827 with its line; long double on this machine; the double-precision budget of heavy dividers.
  4  feature numbers: commensurate triangular boxes, storage of psi6(t) and position snapshots.
usage (from hspist3/):
  python3 validation/gen3_design_numbers_261009.py --measure <mac binary>   (re-runs the six cells on the Mac, ~1 min; writes
                                                                             experiments_gen3_design_261009/eventcounts/)
  python3 validation/gen3_design_numbers_261009.py                           (prints every table of sec. 4.7)
"""
import csv, math, os, re, subprocess, sys, time
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import tests_20260913 as T

PROFILE = os.path.join(HS, "experiments_resched_gate2_261007", "profile_edmd_15008378", "times.tsv")
COUNTS = os.path.join(HS, "experiments_gen3_design_261009", "eventcounts")
CELLS = [("held", "H10", 10.0, 10.0, 50), ("held", "H40", 10.0, 40.0, 200), ("free", "H10", 10.0, 10.0, 50),
         ("free", "H40", 10.0, 40.0, 200), ("dense", "D10", 5.604167, 10.0, 50), ("dense", "D40", 5.604167, 40.0, 200)]
STEPS = {"held": (42000, 60), "free": (12000, 30000), "dense": (6000, 6000)}
DT_SIGMA = 0.4 / 24.0                      # sigma-time per fixed step (tests_20260913.DT_SIGMA)
ETA_CP = math.pi / (2 * math.sqrt(3))
T_TRAJ = 2.0e4
NS_DESIGN = (100, 400, 900, 1600)
ETAS = (("pi/8", math.pi / 8), ("0.70", 0.70), ("0.85", 0.85), ("0.90", 0.90))
SECTION_REF_279282B = {("held", 100): 2.3, ("held", 400): 56.4}     # 261012 sec. 4.3 table (job 14983181), quoted for comparison


# ------------------------------------------------------------------------------------------------------------- measuring
def measure(binp, acc=False):
    """acc=True: the same six commands with --edmd-acc=1 (the accelerated backend, edmd_accelerated.c), timing only (it has no contact
    audit); written to eventcounts/acc_<cell>."""
    os.makedirs(COUNTS, exist_ok=True)
    for kind, tag, L0, H, NS in CELLS:
        d = os.path.join(COUNTS, f"{'acc_' if acc else ''}{kind}_{tag}")
        if os.path.exists(d): sys.exit(f"STOP: {d} exists -- not overwriting")
        os.makedirs(d); hold, post = STEPS[kind]; MF = 1000000000 if kind == "held" else 5 * 2 * NS
        cmd = [binp, "--mode=edmd", "--experiment=energy_transfer", "--headless", "--quiet", "--edmd-acc=0", "--seed-drift-order=drift-first",
               f"--energy-transfer-summary={d}/summary.csv", f"--energy-transfer-trace={d}/tr.csv", "--trace-every=600",
               f"--particles={2 * NS}", f"--particles-boxes={NS},{NS}", "--particle-radius=0.5", f"--l0={L0:.6f}", f"--height={H:.6f}",
               "--num-walls=1", f"--wall-positions={L0:.6f}", f"--wall-mass-factors={MF}", "--wall-thickness=0.05", "--wall-thickness-vis=0.05",
               "--eff-output=wall-ke", f"--wall-hold-steps={hold}", f"--steps={post}", "--fixed-dt=0.4", "--kbt1", "--seed=9700"]
        if acc: cmd[cmd.index("--edmd-acc=0")] = "--edmd-acc=1"
        env = dict(os.environ, HD_CONTACT_AUDIT="1", HD_PISTON_EVENTS=f"{d}/ev.csv")
        t0 = time.time(); rc = subprocess.run(cmd, stdout=open(f"{d}/run.log", "w"), stderr=subprocess.STDOUT, env=env, cwd=d).returncode
        w = time.time() - t0
        nd = nw = 0
        with open(f"{d}/ev.csv") as fh:
            next(fh)
            for l in fh:
                k = l.split(",")[1]
                nd += k.startswith("D"); nw += k.startswith("W")
        m = re.search(r"executed events (\d+)", open(f"{d}/run.log").read())
        if acc and not re.search(r"EDMD backend: accelerated", open(f"{d}/run.log").read()):
            print(f"note: {d}/run.log does not print the backend line in --quiet mode; the flag --edmd-acc=1 was passed")
        ver = subprocess.run([binp, "--version"], capture_output=True, text=True).stdout.split("\n")[0]
        open(f"{d}/counts.tsv", "w").write(f"executed\tdivider\touter_wall\tmac_wall_s\texit\tbinary\n{m.group(1) if m else -1}\t{nd}\t{nw}\t{w:.2f}\t{rc}\t{ver}\n")
        os.remove(f"{d}/ev.csv"); [os.remove(f"{d}/{f}") for f in ("tr.csv",) if os.path.exists(f"{d}/{f}")]
        print(f"{kind} {tag}: executed {m.group(1) if m else '?'}, divider {nd}, outer wall {nw}, {w:.1f} s, exit {rc}")


# ------------------------------------------------------------------------------------------------------------- physics
def z_of(eta):
    """Compressibility factor used for the event-rate model, with its source label."""
    if eta <= 0.69:
        a = np.array([eta]); import plot_speed_of_sound_edmd as S
        return float(S.Z_kolafa_rottner_2006(a)[0]), "KR (module)"
    if eta < 0.72:
        return 9.17 * math.pi / (4 * eta), "Engel plateau P* = 9.17"
    a = ETA_CP / eta - 1
    return 2 / a + 1.90 + 0.67 * a, "Alder-Hoover-Young high-density form (OPEN: not checked against the paper)"


def nu_pair(Z):
    """Pair collisions each disk undergoes per sigma-time (kT = m = sigma = 1)."""
    return 4 * (Z - 1) / math.sqrt(math.pi)


def wall_rate_per_length(eta, Z):
    n = 4 * eta / math.pi                  # number density (sigma = diameter = 1)
    return n * Z / math.sqrt(2 * math.pi)


def geometry(N, eta):
    Ns = N // 2; H = 10 * math.sqrt(N / 100); L0 = Ns * math.pi * 0.25 / (eta * H)
    return Ns, H, L0


def events_per_sigma(N, eta):
    Z, _ = z_of(eta); Ns, H, L0 = geometry(N, eta)
    pair = N * nu_pair(Z) / 2; wl = wall_rate_per_length(eta, Z)
    outer = wl * (4 * L0 + 2 * H); div = wl * 2 * H
    return pair, outer, div


# ------------------------------------------------------------------------------------------------------------- tables
def read_profile():
    t = {}
    for r in csv.DictReader(open(PROFILE), delimiter="\t"):
        t[(r["kind"], int(r["N"]), r["policy"])] = dict(wall=float(r["wall_s"]), div=int(r["div_events"]), wall_ev=int(r["wall_events"]))
    return t


def read_counts():
    c = {}
    for kind, tag, L0, H, NS in CELLS:
        f = os.path.join(COUNTS, f"{kind}_{tag}", "counts.tsv")
        r = list(csv.DictReader(open(f), delimiter="\t"))[0]
        c[(kind, 2 * NS)] = dict(executed=int(r["executed"]), div=int(r["divider"]), wall=int(r["outer_wall"]), mac_s=float(r["mac_wall_s"]),
                                 binary=r["binary"], L0=L0, H=H, Ns=NS)
    return c


def section1():
    prof, cnt = read_profile(), read_counts()
    print("## 1. Baseline: collisions per second per core on KOA (Xeon E5-2680 v2), measured\n")
    print(f"event counts: {sorted({v['binary'] for v in cnt.values()})} on this Mac (same commands as cluster/profile_edmd_koa.sh, seed 9700, "
          "minimal policy, HD_CONTACT_AUDIT=1); KOA wall times: profile job 15008378 (git 7b08827 target koa)\n")
    print("| kind | N | eta | sigma-time | executed events (Mac) | of which pair | divider (Mac / KOA) | outer wall (Mac / KOA) | "
          "pair collisions per disk per sigma-time | KOA s, 7b08827 minimal | KOA s, legacy = 279282b | sec. 4.3 (279282b) s | "
          "events/s 7b08827 | events/s 279282b | pair events/s 7b08827 | pair events/s 279282b |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = {}
    for kind, tag, L0, H, NS in CELLS:
        N = 2 * NS; c = cnt[(kind, N)]; hold, post = STEPS[kind]; Tsig = (hold + post) * DT_SIGMA
        eta = NS * math.pi * 0.25 / (H * L0); pair = c["executed"] - c["div"] - c["wall"]
        pm, pl = prof[(kind, N, "minimal")], prof[(kind, N, "legacy")]
        ref = SECTION_REF_279282B.get((kind, N)); refs = f"{ref:.1f}" if ref else "-"
        rows[(kind, N)] = dict(T=Tsig, eta=eta, ex=c["executed"], pair=pair, div=c["div"], wall=c["wall"], tm=pm["wall"], tl=pl["wall"])
        print(f"| {kind} | {N} | {eta:.4f} | {Tsig:.0f} | {c['executed']} | {pair} | {c['div']} / {pm['div']} | {c['wall']} / {pm['wall_ev']} | "
              f"{2 * pair / (N * Tsig):.3f} | {pm['wall']:.2f} | {pl['wall']:.2f} | {refs} | {c['executed'] / pm['wall']:.3g} | {c['executed'] / pl['wall']:.3g} | "
              f"{pair / pm['wall']:.3g} | {pair / pl['wall']:.3g} |")
    print("\n(collisions each disk undergoes = 2 x pair events / (N x sigma-time). The target of the decision: >= 2e5 events/s per core on KOA; "
          "Engel et al. Table II (Isobe's EDMD): 1.7e9 collisions/h = 4.7e5/s at N = 512^2.)")
    print("\n### The same six cells on this Mac: the default backend (7b08827, minimal policy) against the accelerated backend "
          "(edmd_accelerated.c: 9-cell neighbour predictions and cell-crossing events, but still the global position jump and grid_build per "
          "event and a full reschedule after every divider event; never validated for production)\n")
    print("| kind | N | default backend [s] | accelerated backend [s] | accelerated exit | first failure of the accelerated run (summary.failures.csv) |"
          "\n|---|---|---|---|---|---|")
    for kind, tag, L0, H, NS in CELLS:
        N = 2 * NS; c = read_counts()[(kind, N)]
        fa = os.path.join(COUNTS, f"acc_{kind}_{tag}", "counts.tsv")
        if not os.path.exists(fa): print(f"| {kind} | {N} | {c['mac_s']:.2f} | (not measured) | | |"); continue
        ra = list(csv.DictReader(open(fa), delimiter="\t"))[0]
        ff = os.path.join(COUNTS, f"acc_{kind}_{tag}", "summary.failures.csv"); fail = "none"
        if os.path.exists(ff):
            q = list(csv.DictReader(open(ff)))[0]
            fail = f"t = {float(q['time_sigma']):.2f} sigma-time ({q['phase']}): {q['reason']}, {float(q['value']):.3g} px (tolerance {float(q['limit']):.1e})"
        print(f"| {kind} | {N} | {c['mac_s']:.2f} | {float(ra['mac_wall_s']):.2f}{' (stopped)' if ra['exit'] != '0' else ''} | {ra['exit']} | {fail} |")
    print("\n### The collision-rate model against the measured counts\n")
    print("| kind | N | eta | Z used (source) | model: pair collisions per disk per sigma-time | measured | measured/model | model outer+divider "
          "events per sigma-time | measured | measured/model |\n|---|---|---|---|---|---|---|---|---|---|")
    for (kind, N), r in rows.items():
        Z, src = z_of(r["eta"]); L0 = [c[2] for c in CELLS if c[0] == kind and 2 * c[4] == N][0]; H = [c[3] for c in CELLS if c[0] == kind and 2 * c[4] == N][0]
        nm = nu_pair(Z); nmeas = 2 * r["pair"] / (N * r["T"])
        wl = wall_rate_per_length(r["eta"], Z); wm = wl * (4 * L0 + 2 * H + 2 * H); wmeas = (r["div"] + r["wall"]) / r["T"]
        print(f"| {kind} | {N} | {r['eta']:.4f} | {Z:.3f} ({src}) | {nm:.3f} | {nmeas:.3f} | {nmeas / nm:.3f} | {wm:.1f} | {wmeas:.1f} | {wmeas / wm:.3f} |")
    return rows


def fit_costs(rows):
    """(a) 279282b: t = alpha N E_nondiv + beta N^2 E_div (held cells, legacy times); (b) 7b08827: t = c N^q E_all (held cells, minimal)."""
    h1, h4 = rows[("held", 100)], rows[("held", 400)]
    A = np.array([[100 * (h1["ex"] - h1["div"]), 100 ** 2 * h1["div"]], [400 * (h4["ex"] - h4["div"]), 400 ** 2 * h4["div"]]], float)
    alpha, beta = np.linalg.solve(A, np.array([h1["tl"], h4["tl"]]))
    q = math.log((h4["tm"] / h4["ex"]) / (h1["tm"] / h1["ex"])) / math.log(4)
    c = (h1["tm"] / h1["ex"]) / 100 ** q
    return alpha, beta, q, c


def cost_table(rows):
    alpha, beta, q, c = fit_costs(rows)
    print("\n### Cost models fitted to the held-divider cells, checked on the free and dense cells\n")
    print(f"(a) 279282b: t = alpha N E_non-divider + beta N^2 E_divider: alpha = {alpha:.3e} s, beta = {beta:.3e} s")
    print(f"(b) 7b08827: t = c N^q E_all: q = {q:.3f} (cost per event grows as N^q), c = {c:.3e} s\n")
    print("| kind | N | KOA 279282b [s] | model (a) [s] | KOA 7b08827 [s] | model (b) [s] |\n|---|---|---|---|---|---|")
    for (kind, N), r in rows.items():
        ta = alpha * N * (r["ex"] - r["div"]) + beta * N * N * r["div"]; tb = c * N ** q * r["ex"]
        print(f"| {kind} | {N} | {r['tl']:.2f} | {ta:.2f} | {r['tm']:.2f} | {tb:.2f} |")
    print(f"\n### Time of one {T_TRAJ:.0e} sigma-time trajectory, design geometry (H = 10 sqrt(N/100), L0 from eta, free divider), one KOA core\n")
    print("| eta | Z (source) | N | pair events/sigma | outer-wall/sigma | divider/sigma | events per trajectory | (a) 279282b | (b) 7b08827 | "
          "(c) 2e5/s | (c) 5e5/s |\n|---|---|---|---|---|---|---|---|---|---|---|")
    def fmt(s):
        return f"{s:.0f} s" if s < 120 else (f"{s / 60:.1f} min" if s < 7200 else (f"{s / 3600:.1f} h" if s < 172800 else f"{s / 86400:.1f} d"))
    for lab, eta in ETAS:
        Z, src = z_of(eta)
        for N in NS_DESIGN:
            pe, oe, de = events_per_sigma(N, eta); E = (pe + oe + de) * T_TRAJ
            ta = alpha * N * (pe + oe) * T_TRAJ + beta * N * N * de * T_TRAJ; tb = c * N ** q * E
            print(f"| {lab} | {Z:.2f} ({src.split(' (')[0]}) | {N} | {pe:.3g} | {oe:.3g} | {de:.3g} | {E:.3g} | {fmt(ta)} | {fmt(tb)} | {fmt(E / 2e5)} | {fmt(E / 5e5)} |")


def section3():
    print("\n## 3. Exactness numbers\n")
    print("### Double-precision spacing (ulp) at the engine's times; internal time unit = sigma-time / 24 (dt = 0.4 per step = 1/60 sigma-time)\n")
    print("| absolute time [internal] | = sigma-time | ulp [internal] | position error of one ulp at thermal speed ~1 px per internal unit [px] | [sigma] |\n|---|---|---|---|---|")
    for t in (1e3, 1e4, 4.8e5, 1e6, 1e8, 9e8):
        u = float(np.spacing(t)); print(f"| {t:.1e} | {t / 24:.3g} | {u:.2e} | {u:.2e} | {u / 24:.2e} |")
    print("(4.8e5 = the end of a 2e4 sigma-time trajectory; 1e6-9e8 = the prediction horizons the schedule audit met, sec. 4.4.10 D1. "
          "With a floating origin reset every 1e4 internal units, every stored time stays below ~2e4: ulp <= 3.6e-12)")
    print("\n### The quadratic root near contact: current form t = (-b - sqrt(disc))/vv against the textbook stable form c/(-b + sqrt(disc))\n")
    print("Random oblique approaches (2000 per gap): disk i at a random absolute position in [0, 960] px (the box), j at distance sigma + gap "
          "in a random direction, relative velocity of thermal size with r.v < 0 and impact parameter < sigma; rx = xj - xi etc. formed in "
          "double as the engine does. Reference: the exact root of the SAME double inputs at 50 digits (decimal). Errors in internal time units.\n")
    print("| gap [px] | median exact root | current: median abs. error | current: max abs. error | stable: median abs. error | stable: max abs. error | "
          "floor: ulp(960 px) / closing speed, median |\n|---|---|---|---|---|---|---|")
    import random
    from decimal import Decimal, getcontext
    getcontext().prec = 50
    rng = random.Random(20261009); sig = 24.0
    for gap in (1e-1, 1e-4, 1e-7, 1e-10):
        ec, es, ex, fl = [], [], [], []
        while len(ec) < 2000:
            xi, yi = rng.uniform(0, 960), rng.uniform(0, 960); th = rng.uniform(0, 2 * math.pi)
            xj = xi + (sig + gap) * math.cos(th); yj = yi + (sig + gap) * math.sin(th)
            vx, vy = rng.gauss(0, 1.4), rng.gauss(0, 1.4)
            rx, ry = xj - xi, yj - yi
            b = rx * vx + ry * vy; vv = vx * vx + vy * vy; c = (rx * rx + ry * ry) - sig * sig; disc = b * b - vv * c
            if b >= 0 or disc <= 0 or c <= 0: continue
            cur = (-b - math.sqrt(disc)) / vv; stab = c / (-b + math.sqrt(disc))
            R = [Decimal(v) for v in (rx, ry, vx, vy)]; Bd = R[0] * R[2] + R[1] * R[3]; Vd = R[2] ** 2 + R[3] ** 2
            Cd = R[0] ** 2 + R[1] ** 2 - Decimal(sig) ** 2; Dd = Bd * Bd - Vd * Cd
            if Dd <= 0: continue
            tx = float((-Bd - Dd.sqrt()) / Vd)
            ec.append(abs(cur - tx)); es.append(abs(stab - tx)); ex.append(tx); fl.append(float(np.spacing(960.0)) / (-b / sig))
        print(f"| {gap:.0e} | {np.median(ex):.3e} | {np.median(ec):.1e} | {max(ec):.1e} | {np.median(es):.1e} | {max(es):.1e} | {np.median(fl):.1e} |")
    print("\n### Every bare tolerance literal in edmd_core/edmd.c at 7b08827\n")
    raw = subprocess.run(["git", "show", "7b08827:hspist3/edmd_core/edmd.c"], capture_output=True, text=True, cwd=HS).stdout
    code = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), raw, flags=re.S)       # drop block comments, keep line numbers
    code = re.sub(r"//[^\n]*", "", code); code = re.sub(r'"(?:\\.|[^"\\])*"', '""', code)        # line comments, string literals
    pat = re.compile(r"(?<![\w.])(1e-[0-9]+|1\.0e-[0-9]+)")
    print("(code only: comments and string literals removed before matching)\n")
    print("| line | literal(s) | code |\n|---|---|---|")
    src = raw.split("\n")
    for i, l in enumerate(code.split("\n"), 1):
        m = pat.findall(l)
        if m: print(f"| {i} | {', '.join(sorted(set(m)))} | `{src[i - 1].strip()[:110]}` |")
    fi = np.finfo(np.longdouble)
    print(f"\nlong double on this machine ({os.uname().machine}): {np.dtype(np.longdouble).itemsize * 8} bits stored, eps = {fi.eps:.3e} "
          f"({'= double: no extended precision' if fi.eps == np.finfo(np.float64).eps else 'extended'}); x86-64 (KOA) long double = 80-bit x87, eps = 1.08e-19")
    print("\n### Heavy dividers in double precision (kT = m = 1, positions in px: sigma = 24 px)\n")
    print("| divider mass M | thermal speed sqrt(kT/M) [sigma per sigma-time] | velocity kick per collision ~2 v_gas/M | kick / thermal speed | "
          "ulp of that speed | divider period at N = 100, eta 0.70 [sigma-time] (cot K = alpha K, c_s 15) | at N = 1600 |\n|---|---|---|---|---|---|---|")
    for M in (1e2, 1e4, 1e6, 4e7, 1e8):
        vth = 1 / math.sqrt(M); kick = 2 / M
        per = []
        for N in (100, 1600):
            Ns, H, L0 = geometry(N, 0.70); al = M / (2 * Ns); K = T.k_root(al); Le = L0 - 1.025
            per.append(2 * math.pi * Le / (15.0 * K))
        print(f"| {M:.0e} | {vth:.2e} | {kick:.1e} | {kick / vth:.1e} | {float(np.spacing(vth)):.1e} | {per[0]:.3g} | {per[1]:.3g} |")


def section4():
    print("\n## 4. Feature numbers\n")
    print("### Commensurate triangular crystals in one compartment (hard walls; lattice constant a = sqrt(eta_cp/eta_lattice) sigma)\n")
    print("The outermost rows sit at a distance r + g from each wall and from the divider face (r = 1/2, g = surface gap, here g = 0: "
          "rows touching). PARALLEL = columns of disks along the divider (y), spaced a sqrt(3)/2 in x, alternate columns shifted by a/2: "
          "L0 - t/2 = (n_x - 1) a sqrt(3)/2 + 1 + 2g, H = (n_y - 1/2) a + 1 + 2g. PERPENDICULAR = rows along x, spaced a sqrt(3)/2 in y: "
          "L0 - t/2 = (n_x - 1/2) a + 1 + 2g, H = (n_y - 1) a sqrt(3)/2 + 1 + 2g. N_s = n_x n_y (no vacancy). The nominal eta "
          "(N_s pi/4 / (H L0), the project's definition) is below the lattice's own eta because of the wall layers. The box is also "
          "rounded to the 1/24 sigma pixel grid (methods sec. 14); the column 'strain' is the lattice strain that rounding leaves.\n")
    print("| eta_lattice | orientation | n_x x n_y | N_s | a | L0 exact | H exact | nominal eta | L0 on the 1/24 grid | strain from the grid [%] |"
          "\n|---|---|---|---|---|---|---|---|---|---|")
    tw, g = 0.05, 0.0
    for etal in (0.72, 0.80, 0.85, 0.88, 0.90):
        a = math.sqrt(ETA_CP / etal)
        for orient in ("parallel", "perpendicular"):
            for nx, ny in ((5, 10), (10, 20), (15, 30), (20, 40)):
                if orient == "parallel": Lin = (nx - 1) * a * math.sqrt(3) / 2; H = (ny - 0.5) * a + 1 + 2 * g
                else: Lin = (nx - 0.5) * a; H = (ny - 1) * a * math.sqrt(3) / 2 + 1 + 2 * g
                L0 = Lin + 1 + 2 * g + tw / 2; L0g = math.floor(2 * L0 * 24) / 48
                etan = nx * ny * math.pi / 4 / (H * L0)
                print(f"| {etal:.2f} | {orient} | {nx} x {ny} | {nx * ny} | {a:.5f} | {L0:.6f} | {H:.6f} | {etan:.4f} | {L0g:.6f} | "
                      f"{100 * (L0g - L0) / Lin:+.3f} |")
    print("\n(the grid rounding is the binary's floor of 2 L0 x 24 px (box_delta, methods sec. 14); an exact box length removes it. "
          "Vacancies: N_s = n_x n_y - n_vac, positions removed by a stated rule)")
    print("\n### Storage of the structural clock and snapshots per trajectory (text, ~18 bytes per number incl. separator)\n")
    print("| N | psi6(t) every 0.25 sigma-time, 2e4 sigma-time | positions every 10 sigma-time | positions every 1 sigma-time |\n|---|---|---|---|")
    for N in NS_DESIGN:
        print(f"| {N} | {4 * T_TRAJ * 2 * 18 / 1e6:.2f} MB | {T_TRAJ / 10 * N * 2 * 18 / 1e6:.1f} MB | {T_TRAJ * N * 2 * 18 / 1e6:.0f} MB |")


def main():
    if "--measure" in sys.argv:
        return measure(os.path.abspath(sys.argv[sys.argv.index("--measure") + 1]))
    if "--measure-acc" in sys.argv:
        return measure(os.path.abspath(sys.argv[sys.argv.index("--measure-acc") + 1]), acc=True)
    print("# Generation 3 -- design numbers (261012 sec. 4.7), printed by validation/gen3_design_numbers_261009.py\n")
    rows = section1(); cost_table(rows); section3(); section4()


if __name__ == "__main__":
    main()
