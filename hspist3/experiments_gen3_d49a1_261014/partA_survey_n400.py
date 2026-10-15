#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 49 addendum 1, items 1-3 and 5; 261012 sec. 4.7.111): THE PART A SURVEY AT N = 400, read-only.
Prints (A) the code facts of the held path (file:line at f72d764 and at the engine-gen3 head, by git show); (B) the stencil per N = 400 cell
by the registered formulas, which are the functions of validation/paper1_confinement_prereg_20261012.py (Table A of notes sec. 1.4b):
sigma_x = L_eff / sqrt(2 N_s c_KR^2), dL = sigma_x / 2, the old grid rule dL = max(1/24, grid-rounded sigma_x / 2), the span in sigma_x and
the KR truncation bias of the five-point stencil, flagged above 0.1 %; also L0 rounded to the 1/48 grid (what the held path needs) and its eta;
(C) the event log per (seed, position) at N = 400 (bulk estimate x the factor measured by partA_probes_n400.py); (D) the cost in E5
core-hours: the stage 1 generator's model (gen_stage1_koa_261012: events per sigma-time / rate_old) x the median measured ratio of the
cell's E5-2670 v3 chunks (interim_261013/timing/timing_<group>.tsv); (E) k_S^dyn per N = 400 cell from B400's INTERIM tables.
Nothing is run; nothing of Part A is planned here beyond the registered formulas. usage (from hspist3/):
  python3 experiments_gen3_d49a1_261014/partA_survey_n400.py
"""
import csv, importlib.util, math, os, re, statistics, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); REPO = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS); sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage1_261012"))
spec = importlib.util.spec_from_file_location("prereg", os.path.join(HS, "validation", "paper1_confinement_prereg_20261012.py"))
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
import gen_stage1_koa_261012 as GEN

N, NS, H = 400, 200, 20.0
EPS0 = 0.2710            # notes sec. 1.4b: sigma_F/F = 0.2710 (100/N_s)^(1/2)/sqrt(T), Level 3 c0, eta 0.10005 [PLANNING stand-in for eps0+]
T_SEED, MIN_SEEDS, MAX_SEEDS, NOISE = 5000.0, 8, 400, 0.009
TAIL = 20.0              # the released tail [sigma-time]
LOG_FACTOR = 1.12        # rows measured / bulk estimate at N = 400 (partA_probes_n400.py: 1.10 at eta 0.05, 1.11 at pi/8, 1.12 at 0.6905)
BYTES_ROW = 65.3
FALLBACK = []           # cells whose ratio had to use every CPU model


def grep(rev, path, pat):
    try:
        txt = subprocess.run(["git", "-C", REPO, "show", f"{rev}:{path}"], capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return [f"{rev} {path}: absent"]
    return [f"{rev[:8]} {path}:{i}: {l.strip()[:150]}" for i, l in enumerate(txt.splitlines(), 1) if re.search(pat, l)]


def facts():
    head = subprocess.run(["git", "-C", REPO, "rev-parse", "--short=8", "HEAD"], capture_output=True, text=True).stdout.strip()
    print(f"## (A) The held path's code at f72d764 and at the engine-gen3 head {head}\n")
    items = [
        ("(a) the start and the exact box: refused outside the speed-of-sound experiment (the AF command is --experiment=energy_transfer)",
         "hspist3/00ALLINONE.c", r"else if \(\(cli_gen3_seeding \|\| cli_gen3_exact_box|wired to the speed-of-sound experiment only\";$"),
        ("(b) the box on the 1/48 grid, the divider as a float fraction of the box width, the position check at 1e-5 sigma",
         "hspist3/00ALLINONE.c", r"not on the 1/48 sigma grid|float frac = pos_sigma / fmaxf|all_wall_positions\[w\] = XW1 \+ frac \* width_px;|"
                                 r"position_tolerance_sigma = 1e-5;|initial_wall_position_mismatch"),
        ("(c) the event log under gen3: outer-wall rows and divider/piston rows; u = the body's velocity, dE = the particle's KE change if M = 0",
         "hspist3/edmd_core/edmd_gen3.c", r"evlog_row\(S, WALL_NAME|evlog_row\(S, kb|static void evlog_row"),
        ("(c), (e) the log's switch: HD_PISTON_EVENTS, energy-transfer mode only", "hspist3/00ALLINONE.c", r"getenv\(\"HD_PISTON_EVENTS\"\)|energy-transfer mode only, off unless"),
        ("(c), (d), (e) the AF command: hold (equilibration 200 + record 5000 sigma-time), the released tail, t0 = 200 passed to reduce_AF.py",
         "hspist3/cluster/confinement_20261013/conf_worker.sh", r"AF <rel_pos_dir>|hold = 12000|--wall-hold-steps=\$hold|reduce_AF.py\" .* 200 "),
        ("(c), (e) reduce_AF.py: the whole log read, the D0 rows of [t0, t1), F_L and F_R, u_wall_max and W_div",
         "hspist3/cluster/confinement_20261013/reduce_AF.py", r"pd.read_csv\(ev|kind\"\] == \"D0\"|FL = |u_wall_max="),
        ("(c) the drift f from the temperatures (C4 formula; |1 - f| < 0.002)", "hspist3/validation/paper1_confinement_heldwall_posthoc_261004.py", r"^def drift"),
        ("(e) a streaming reader of the log (gate G6 tool)", "hspist3/cluster/gen3_gate_261009/af_stream_reduce.py", r"data = fh.read\(\)|^def "),
    ]
    for what, path, pat in items:
        print(f"- {what}:")
        for rev in ("f72d764", "HEAD"):
            for l in grep(rev, path, pat): print(f"    {l.replace('HEAD    ', head + ' ')}")
    print()


def stencil_table():
    print("## (B) The stencil per N = 400 cell (H = 20, N_s = 200 per side, t = 0.05; L_eff = L0 - 2r - t/2; KR: Z, Z', c_s^2 = Z + eta Z' + Z^2)\n")
    print("| eta | L0 (registered) | L0 on the 1/48 grid? | L0 rounded to 1/48 | its eta | d eta / eta | sigma_x | dL = sigma_x/2 | bias, exact dL | "
          "grid rule dL | span +-2dL [sigma_x] | bias, grid dL | flag (> 0.1 %) |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for e in GEN.ETA_B:
        L0 = N * math.pi / (8.0 * H * e); og = abs(L0 * 48 - round(L0 * 48)) < 1e-9; L0g = round(L0 * 48) / 48; eg = P.eta_of(NS, H, L0g)
        sx = math.sqrt(P.leff(L0) ** 2 / (2 * NS * P.cs2(e))); dL = sx / 2; dLg = max(1 / 24, round(sx / 2 * 24) / 24)

        def bias(d):
            F = lambda L: NS * P.Z(P.eta_of(NS, H, L)) / L
            ex = NS * (P.Z(e) + e * P.dZ(e)) / L0 ** 2
            st5 = -(F(L0 - 2 * d) - 8 * F(L0 - d) + 8 * F(L0 + d) - F(L0 + 2 * d)) / (12 * d)
            return st5 / ex - 1
        b, bg = bias(dL), bias(dLg)
        flag = ("grid" if abs(bg) > 0.001 else "") + (" exact" if abs(b) > 0.001 else "")
        print(f"| {e:.4g} | {L0:.6f} | {'yes' if og else 'no'} | {L0g:.6f} | {eg:.6f} | {100 * (eg / e - 1):+.3f} % | {sx:.4f} | {dL:.4f} | {b:+.1e} | "
              f"{dLg:.4f} ({round(dLg * 24)}/24) | {2 * dLg / sx:.2f} | {bg:+.1e} | {flag.strip() or '-'} |")
        rows.append(dict(eta=e, L0=L0, sx=sx, dL=dL, dLg=dLg))
    print()
    return rows


def log_table():
    print(f"## (C) The event log per (seed, position) at N = 400: held T_eq + 5000 sigma-time + the {TAIL:g} sigma-time tail\n")
    print(f"rows per sigma-time = {LOG_FACTOR} x the bulk estimate rho Z_KR / sqrt(2 pi) x (4 H + 4 L0); the divider's share 2 H / (4 H + 4 L0); "
          f"{BYTES_ROW} bytes per row (partA_probes_n400.py)\n")
    print("| eta | L0 | T_eq | rows per sigma-time | rows | log size | divider share | a divider-only log would be |\n|---|---|---|---|---|---|---|---|")
    out = {}
    for e in GEN.ETA_B:
        L0 = N * math.pi / (8.0 * H * e); Teq = 1000.0 if e <= 0.60 + 1e-12 else 5000.0
        rate = LOG_FACTOR * (NS / (L0 * H)) * P.Z(e) / math.sqrt(2 * math.pi) * (4 * H + 4 * L0); T = Teq + T_SEED + TAIL
        rows = rate * T; mb = rows * BYTES_ROW / 1e6; sh = 2 * H / (4 * H + 4 * L0)
        out[e] = mb
        print(f"| {e:.4g} | {L0:.4f} | {Teq:g} | {rate:.1f} | {rows:,.0f} | {mb:.0f} MB | {100 * sh:.1f} % | {mb * sh:.0f} MB |")
    print()
    return out


def measured_ratio(group, cid):
    p = os.path.join(HS, "experiments_gen3_s1_261012", "interim_261013", "timing", f"timing_{group}.tsv")
    rows = [x for x in csv.DictReader(open(p), delimiter="\t") if x["cell"] == cid and x["ratio"] not in ("", "nan")]
    r = [float(x["ratio"]) for x in rows if "E5-2670 v3" in x["cpu"]]
    if not r: FALLBACK.append(cid); r = [float(x["ratio"]) for x in rows]   # no E5-2670 v3 chunk: every CPU model (printed)
    return statistics.median(r) if r else math.nan


def cost_tables(st):
    print("## (D) The cost in E5 core-hours (the generator's model x the median measured ratio of the cell's E5-2670 v3 chunks)\n")
    print("| eta | T_eq | model s per sigma-time (N 400) | measured ratio (B400) | E5 s per sigma-time | per (seed, position) [core-h] | "
          "the N 400 pilot, 4 seeds x 5 positions [core-h] | seeds per position (eps0 0.2710, PLANNING) | full N 400 Part A [core-h] |\n"
          "|---|---|---|---|---|---|---|---|---|")
    tot_pilot = tot_full = 0.0; per_cell = {}
    for row in st:
        e = row["eta"]; Teq = 1000.0 if e <= 0.60 + 1e-12 else 5000.0
        model = GEN.events_per_sigma(N, e) / GEN.rate_old(N, e); ratio = measured_ratio("B400", GEN.cell_id("B", N, e, "default"))
        sps = model * ratio; ch = sps * (Teq + T_SEED + TAIL) / 3600
        L0 = row["L0"]; dL = row["dL"]
        F = lambda L: NS * P.Z(P.eta_of(NS, H, L)) / L
        exact = NS * (P.Z(e) + e * P.dZ(e)) / L0 ** 2
        eps = EPS0 * math.sqrt(100 / NS); Tpos = ((math.sqrt(130) / 12) * eps * (F(L0) / exact) / (dL * NOISE)) ** 2 / 2
        seeds = min(MAX_SEEDS, max(MIN_SEEDS, math.ceil(Tpos / T_SEED)))
        pilot = 4 * 5 * ch; full = seeds * 5 * ch; tot_pilot += pilot; tot_full += full; per_cell[e] = (ch, seeds)
        print(f"| {e:.4g} | {Teq:g} | {model:.4f} | {ratio:.3f} | {sps:.4f} | {ch:.3f} | {pilot:.1f} | {seeds} (T_pos {Tpos:.3g}) | {full:.1f} |")
    print(f"\nthe N 400 pilot: 16 cells x 5 positions x 4 seeds = 320 trajectories, {tot_pilot:.0f} E5 core-hours; "
          f"the full N 400 Part A at the planning eps0: {tot_full:.0f} E5 core-hours\n")
    # the registered pilot: N 100 and N 1600
    for n, grp, proxy in ((100, "B100", None), (1600, "B900", "B900's measured ratios as a stand-in (B1600 not fetched)")):
        Hn = 10.0 * math.sqrt(n / 100.0); tot = 0.0
        for e in GEN.ETA_B:
            Teq = 1000.0 if e <= 0.60 + 1e-12 else 5000.0
            model = GEN.events_per_sigma(n, e) / GEN.rate_old(n, e)
            ratio = measured_ratio(grp, GEN.cell_id("B", 900 if grp == "B900" else n, e, "default"))
            tot += 4 * 5 * model * ratio * (Teq + T_SEED + TAIL) / 3600
        print(f"the registered pilot at N {n}: 16 x 5 x 4 = 320 trajectories, {tot:.0f} E5 core-hours" + (f" ({proxy})" if proxy else ""))
    if FALLBACK: print("(cells without an E5-2670 v3 chunk, their ratio from every CPU model: " + ", ".join(FALLBACK) + ")")
    print()
    return per_cell


def ksdyn_table():
    print("## (E) k_S^dyn per N = 400 cell (sec. 4.8.6.5; B400's INTERIM tables, interim_261013/B400/a_interim_tables.txt)\n")
    p = os.path.join(HS, "experiments_gen3_s1_261012", "interim_261013", "B400", "a_interim_tables.txt")
    print("| eta | k_S^dyn +- SE [kT/sigma^2] | relative SE | N_s c_KR^2 / L_eff^2 (KR, for scale) | k_S^dyn / that |\n|---|---|---|---|---|")
    for l in open(p):
        if not l.startswith("| B_N0400"): continue
        c = [x.strip() for x in l.split("|")[1:-1]]
        e = float(re.search(r"eta([\d.]+)", c[0]).group(1)); ks = c[12]
        m = re.match(r"([\d.e+-]+) \+- ([\d.e+-]+)", ks)
        L0 = N * math.pi / (8.0 * H * e); kr = NS * P.cs2(e) / P.leff(L0) ** 2
        if m:
            v, s = float(m.group(1)), float(m.group(2))
            print(f"| {e:.4g} | {v:.6g} +- {s:.2g} | {100 * s / v:.3f} % | {kr:.6g} | {v / kr:.4f} |")
        else:
            print(f"| {e:.4g} | {ks} | - | {kr:.6g} | - |")
    print()


def main():
    print("# Part A survey at N = 400 (partA_survey_n400.py; decision 49 addendum 1; read-only)\n")
    facts(); st = stencil_table(); sizes = log_table(); per_cell = cost_tables(st)
    sh = {e: 2 * H / (4 * H + 4 * N * math.pi / (8.0 * H * e)) for e in GEN.ETA_B}
    pil = sum(20 * sizes[e] for e in GEN.ETA_B); full = sum(5 * per_cell[e][1] * sizes[e] for e in GEN.ETA_B)
    print("## (C') The event-log volume at N = 400, written by the AF command as it is (one full log per (seed, position))\n")
    print(f"- the N 400 pilot (16 cells x 5 positions x 4 seeds): {pil / 1e3:.1f} GB; a divider-only log: {sum(20 * sizes[e] * sh[e] for e in GEN.ETA_B) / 1e3:.1f} GB")
    print(f"- the full N 400 Part A at the planning seeds of (D): {full / 1e6:.2f} TB; a divider-only log: "
          f"{sum(5 * per_cell[e][1] * sizes[e] * sh[e] for e in GEN.ETA_B) / 1e6:.2f} TB")
    print(f"- the largest single log: {max(sizes.values()):.0f} MB (eta 0.69); reduce_AF.py reads it whole\n")
    ksdyn_table()


if __name__ == "__main__":
    main()
