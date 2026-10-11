#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task M1; 261012 sec. 4.7.59): the TABLES behind the four figures for the meeting of
Wednesday 2026-10-14. Existing data only, no runs. m1_figures_261011.py draws the figures from the CSV files written here, and from
nothing else.
SOURCES. Every number comes from a printed table that is committed. It is read with git from the commit named in SRC below, so neither
an edit of the working tree nor the branch that is checked out changes it. The one exception is Figure 1's nu_d: it needs the fetched
per-trajectory ACFs (acf_runs.npz, 1 MB each, not committed). Their SHA-256 are first checked against the committed fetch manifest.
nu_d is then computed by the registered Test G verdict's own procedure (information_nud of validation/testG_verdict_261009.py:
resched_testT_followup_261007.fit_acf, started at the mass's pooled argmax mean of both engines). The per-trajectory table is written
too, so every Figure 1 number can be recomputed from committed tables.
  Fig. 1  Test G on KOA (registered; ADDITION: nu_d, information): per engine and mass, the mean nu by argmax (red_nu.csv) and by nu_d,
          with SE = SD/sqrt(n). gen3 - gen2 per mass with the verdict's own stats() (imported unchanged): relative [%], its 95 % interval
          100 (d -+ 1.96 s) / mean_gen2, and z. Masses 50, 300 and 1500 (the verdict's, 400 per engine), and 2000 (information, 100 per
          engine).
  Fig. 2  Report G2's printed delta_M table (POST HOC, INFORMATION): estimator, engine, c_ref, delta_M [%] with SE; plus methods sec.
          16's reference, as G2 prints it.
  Fig. 3  P6's printed per-seed tables (PILOT): per compartment Z_c and, where measured, |Psi_c| and the orientation. At N = 1600 the
          compressed starts' compartments are put into the two Z groups by the cut that P6 printed (its largest gap). P6 printed no
          split at N = 400, so none is made there.
  Fig. 4  K2-2 part (c): the gen3 rate table on KOA's Gold node, as printed. Decision 19 item 2's window table (the Mac), with gen3 /
          gen2 per N and time scale.
CHECKS. Each one is printed, and the script stops on any failure:
  - the SHA-256 of the fetched ACFs equal the committed manifest's;
  - Fig. 1 reproduces, at the printed precision, every cell of Report G-KOA's registered table, its nu_d information table and its
    M = 2000 row, and G2's eight printed mean nu_d;
  - Fig. 3 reproduces P6's per-start mean (SE) and the two N = 1600 groups' counts and mean (SE), within the rounding of the printed
    per-seed values (1e-4);
  - Fig. 4: 12 rate rows, every run clean; 8 window rows, gen3 clean.
usage (from hspist3/ of main, where the fetched ACFs are): python3 experiments_m1_figures_261011/m1_tables_261011.py [--out DIR]
"""
import argparse, hashlib, io, math, os, re, subprocess, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); ROOT = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
KOA = "hspist3/experiments_gen3_gate_koa_261009"
SRC = {   # name: (commit, path) -- every committed table this script reads
    "verdict": ("2eff7c0d5c42008f84cd2066c2040e723d54c5f0", f"{KOA}/REPORT_G_KOA_261011/verdict_testG_koa_output.txt"),
    "fetched": ("2eff7c0d5c42008f84cd2066c2040e723d54c5f0", f"{KOA}/REPORT_G_KOA_261011/FETCHED_SHA256_FILES.txt"),
    "rates": ("2eff7c0d5c42008f84cd2066c2040e723d54c5f0", f"{KOA}/REPORT_G_KOA_261011/koa_k2_m5_rates_output.txt"),
    "g2": ("76f1f1157117b88894cc6671e1ce15cbced2c392", f"{KOA}/REPORT_G2_261011/testG_koa_posthoc_output.txt"),
    "p6": ("c540b7460b3b3d43d25cd9f1702e9eb9664b706e", "hspist3/experiments_gen3_p6_261010/p6_tables_output.txt"),
    "d19": ("cb2d32cef319704f58942251ba29de8413d53336", "hspist3/experiments_gen3_d19_261010/d19_tables_output.txt"),
}
RED_COMMIT = "2eff7c0d5c42008f84cd2066c2040e723d54c5f0"           # the fetched red_nu.csv files, committed with Report G-KOA
CELL, ENG, MASSES, VERDICT_MASSES = "epi8_H_H10_L10", ("gen2", "gen3"), (50, 300, 1500, 2000), (50, 300, 1500)
Z95 = 1.96
FAILS = []


def stop(msg):
    print(f"STOP: {msg}"); sys.exit(2)


def check(label, ok, detail=""):
    print(f"- {label}: {'PASS' if ok else 'FAIL'}{(' -- ' + detail) if detail else ''}")
    if not ok: FAILS.append(label)


def git(*a):
    return subprocess.run(["git", "-C", ROOT, *a], check=True, capture_output=True, text=True).stdout


def committed(name):
    c, p = SRC[name]
    return git("show", f"{c}:{p}")


def md_rows(text, header):
    """the rows of the one markdown table in `text` whose header line starts with `header`"""
    lines = text.splitlines()
    at = [i for i, l in enumerate(lines) if l.startswith(header)]
    if len(at) != 1: stop(f"table '{header}': {len(at)} found, 1 expected")
    rows = []
    for l in lines[at[0] + 2:]:
        if not l.startswith("|"): break
        rows.append([c.strip() for c in l.strip().strip("|").split("|")])
    return rows


def pm(s):
    """'-0.689 +- 0.146' -> (-0.689, 0.146)"""
    m = re.fullmatch(r"([-+]?[\d.]+) \+- ([\d.]+)", s.strip())
    if not m: stop(f"not 'value +- SE': {s!r}")
    return float(m.group(1)), float(m.group(2))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- Fig. 1: Test G on KOA
def fig1(out):
    import testG_verdict_261009 as V                             # the registered rule's stats() and outcome(), unchanged
    import resched_testT_followup_261007 as F                    # fit_acf: sec. 4.4.13 item 4a, as information_nud uses it
    import warnings; warnings.filterwarnings("ignore")
    man = {}
    for l in committed("fetched").splitlines():
        p = l.split()
        if len(p) == 3: man[p[2]] = p[0]
    print("\n## Fig. 1: Test G on KOA\n")
    print(f"red_nu.csv read from commit {RED_COMMIT[:8]}; acf_runs.npz from the fetched tree, SHA-256 against {SRC['fetched'][1]} at "
          f"{SRC['fetched'][0][:8]}\n")
    traj = []
    for M in MASSES:
        red = {e: pd.read_csv(io.StringIO(git("show", f"{RED_COMMIT}:{KOA}/testG/{e}/{CELL}/m_{M}/red_nu.csv"))).sort_values("run")
               .reset_index(drop=True) for e in ENG}
        if list(red["gen2"]["run"]) != list(red["gen3"]["run"]): stop(f"M = {M}: the two engines' run lists differ")
        ref = 0.5 * (red["gen3"]["nu"].mean() + red["gen2"]["nu"].mean())          # information_nud's start value
        for e in ENG:
            p = os.path.join(HS, "experiments_gen3_gate_koa_261009", "testG", e, CELL, f"m_{M}", "acf_runs.npz")
            key = f"./testG/{e}/{CELL}/m_{M}/acf_runs.npz"
            if not os.path.exists(p): stop(f"{p} is missing (the fetched Test G tree)")
            got = sha256(p); check(f"{e} M = {M}: acf_runs.npz SHA-256 = the manifest's", man.get(key) == got, got[:16])
            z = np.load(p); dt = float(red[e]["dt"].iloc[0])
            for r, s, nu in zip(red[e]["run"], red[e]["seed"], red[e]["nu"]):
                try: nd = F.fit_acf(z[f"run{r}"].astype(float), dt, ref)[4] / (2 * math.pi)
                except Exception: nd = math.nan
                traj.append(dict(engine=e, M=M, run=int(r), seed=int(s), nu_arg=float(nu), nu_d=nd))
    T = pd.DataFrame(traj)
    T.to_csv(os.path.join(out, "fig1_testG_koa_per_trajectory.csv"), index=False)
    means, diffs = [], []
    for M in MASSES:
        for est, col in (("argmax", "nu_arg"), ("nu_d", "nu_d")):
            a = T[(T.engine == "gen3") & (T.M == M)][col].to_numpy(float); b = T[(T.engine == "gen2") & (T.M == M)][col].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]                # as information_nud: pairs by run
            for e, v in (("gen3", a), ("gen2", b)):
                se = v.std(ddof=1) / math.sqrt(len(v))
                means.append(dict(engine=e, M=M, estimator=est, n=len(v), mean=v.mean(), se=se, lo95=v.mean() - Z95 * se,
                                  hi95=v.mean() + Z95 * se, role="registered" if (est == "argmax" and M in VERDICT_MASSES) else "information",
                                  fit_failures=int((~ok).sum()) if est == "nu_d" else 0))
            r = V.stats(a, b)
            diffs.append(dict(M=M, estimator=est, n_gen3=r["na"], n_gen2=r["nb"], mean_gen3=r["ma"], se_gen3=r["sa"], mean_gen2=r["mb"],
                              se_gen2=r["sb"], difference=r["d"], rel_pct=r["rel"], lo95_pct=r["lo"], hi95_pct=r["hi"], z=r["z"],
                              corr_same_seed=np.corrcoef(a, b)[0, 1],
                              outcome=V.outcome(r["z"]) if (est == "argmax" and M in VERDICT_MASSES) else "",
                              role="registered" if (est == "argmax" and M in VERDICT_MASSES) else "information"))
    Mn, D = pd.DataFrame(means), pd.DataFrame(diffs)
    Mn.to_csv(os.path.join(out, "fig1_testG_koa_means.csv"), index=False)
    D.to_csv(os.path.join(out, "fig1_testG_koa_difference.csv"), index=False)
    print("\n| M | estimator | gen3 mean (SE) | gen2 mean (SE) | gen3 - gen2 [%] | 95 % interval [%] | z | role |\n|---|---|---|---|---|---|---|---|")
    for _, d in D.iterrows():
        print(f"| {d.M} | {d.estimator} | {d.mean_gen3:.7f} ({d.se_gen3:.2e}) | {d.mean_gen2:.7f} ({d.se_gen2:.2e}) | {d.rel_pct:+.3f} | "
              f"[{d.lo95_pct:+.3f}, {d.hi95_pct:+.3f}] | {d.z:+.2f} | {d.role}{(': ' + d.outcome) if d.outcome else ''} |")
    print(f"\nnu_d fit failures: {int(T.nu_d.isna().sum())} of {len(T)}\n")
    # the checks against the printed reports, at the printed precision
    vt = committed("verdict")
    for row in md_rows(vt, "| mass | judged on (per engine) |"):
        M = int(row[0]); d = D[(D.M == M) & (D.estimator == "argmax")].iloc[0]
        mine = [str(M), "400", str(d.n_gen3), str(d.n_gen2), f"{d.mean_gen3:.7f}", f"{d.se_gen3:.2e}", f"{d.mean_gen2:.7f}", f"{d.se_gen2:.2e}",
                f"{d.difference:+.3e}", f"{d.rel_pct:+.3f}", f"[{d.lo95_pct:+.3f}, {d.hi95_pct:+.3f}]", f"{d.z:+.2f}", f"**{d.outcome}**"]
        check(f"Report G-KOA's registered row M = {M} reproduced", row == mine, "" if row == mine else f"printed {row} / here {mine}")
    for row in md_rows(vt, "| M | estimator | gen3 mean | gen2 mean |"):
        M = int(row[0]); est = "argmax" if row[1].startswith("argmax") else "nu_d"; d = D[(D.M == M) & (D.estimator == est)].iloc[0]
        mine = [f"{d.mean_gen3:.7f}", f"{d.mean_gen2:.7f}", f"{d.rel_pct:+.3f}", f"{d.z:+.2f}", f"{d.corr_same_seed:+.3f}"]
        check(f"Report G-KOA's information row M = {M}, {est} reproduced", row[2:] == mine, "" if row[2:] == mine else f"printed {row[2:]} / here {mine}")
    row = [r for r in md_rows(vt, "| row | n gen3 | n gen2 |") if r[0] == "nu, M = 2000"][0]
    d = D[(D.M == 2000) & (D.estimator == "argmax")].iloc[0]
    mine = [str(d.n_gen3), str(d.n_gen2), f"{d.mean_gen3:.7f}", f"{d.se_gen3:.2e}", f"{d.mean_gen2:.7f}", f"{d.se_gen2:.2e}", f"{d.rel_pct:+.3f}",
            f"[{d.lo95_pct:+.3f}, {d.hi95_pct:+.3f}]", f"{d.z:+.2f}"]
    check("Report G-KOA's information row 'nu, M = 2000' reproduced", row[1:10] == mine, "" if row[1:10] == mine else f"printed {row[1:10]} / here {mine}")
    for row in md_rows(committed("g2"), "| engine | M | trajectories | nu_d fit failures |"):
        e, M = row[0], int(row[1]); printed = row[6].split(" +- ")[0]
        here = f"{Mn[(Mn.engine == e) & (Mn.M == M) & (Mn.estimator == 'nu_d')].iloc[0]['mean']:.7f}"
        check(f"G2's mean nu_d, {e} M = {M}, reproduced", printed == here, f"{printed} / {here}")


# ---------------------------------------------------------------- Fig. 2: delta_M per engine (Report G2)
def fig2(out):
    g2 = committed("g2"); rows = []
    print("\n## Fig. 2: delta_M per engine, Report G2's printed table (POST HOC, INFORMATION)\n")
    for r in md_rows(g2, "| estimator | engine | c_ref | delta_50 [%] |"):
        for M, cell in zip(MASSES, r[3:7]):
            v, s = pm(cell)
            rows.append(dict(estimator=r[0], engine=r[1], c_ref=float(r[2]), M=M, delta_pct=v, se_pct=s, lo95_pct=v - Z95 * s, hi95_pct=v + Z95 * s,
                             in_reference=M >= 300, source="Report G2 (76f1f115), testG_koa_posthoc_output.txt"))
    m = re.search(r"methods sec\. 16: ([-+]?[\d.]+) \+- ([\d.]+) % with nu_0 \(([^)]*)\)", g2)
    if not m: stop("G2: the methods sec. 16 line is missing")
    v, s = float(m.group(1)), float(m.group(2))
    rows.append(dict(estimator="nu_0", engine="gen2 (methods sec. 16)", c_ref=math.nan, M=50, delta_pct=v, se_pct=s, lo95_pct=v - Z95 * s,
                     hi95_pct=v + Z95 * s, in_reference=False, source=f"methods sec. 16 as Report G2 prints it: {m.group(3)}"))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(out, "fig2_deltaM_G2.csv"), index=False)
    check("G2's delta table: 3 estimators x 2 engines x 4 masses, plus methods sec. 16", len(D) == 25, f"{len(D)} rows")
    print("\n| estimator | engine | c_ref | " + " | ".join(f"delta_{M} [%]" for M in MASSES) + " |\n|---|---|---|---|---|---|---|")
    for (est, e), g in D[D.engine.isin(ENG)].groupby(["estimator", "engine"], sort=False):
        print(f"| {est} | {e} | {g.c_ref.iloc[0]:.5f} | " + " | ".join(f"{x:+.3f} +- {s:.3f}" for x, s in zip(g.delta_pct, g.se_pct)) + " |")
    print(f"\nmethods sec. 16 (gen2): {v:+.3f} +- {s:.3f} % ({m.group(3)})")


# ---------------------------------------------------------------- Fig. 3: the solid, Z per compartment (P6's tables)
STARTS = (  # P6's row label -> key, build, window, kind
    ("P6 turned (ordered)", "P6 turned", "f72d764", "Q4", "turned start"),
    ("P6 0-deg (ADDITION)", "P6 0-deg", "f72d764", "Q4", "0-deg start"),
    ("P4 A, the present family (0 deg)", "P4 A", "06e7248", "Q4", "present M4 family"),
    ("P4 A, the present family (turned)", "P4 A", "06e7248", "Q4", "present M4 family"),
    ("P4a A, first quarter [0, 5000) (all at 0 deg)", "P4a A [0, 5000)", "06e7248", "[0, 5000)", "present M4 family"),
    ("P4 B, compressed u = 0.01", "P4 B u 0.01", "06e7248", "Q4", "compressed"),
    ("P5 B, compressed u = 0.001", "P5 B u 0.001", "06e7248", "Q4", "compressed"),
    ("P5 B, compressed u = 0.0001", "P5 B u 0.0001", "06e7248", "Q4", "compressed"),
)
CELLS = {"N400_eta0.740": (400, 0.740, 20.0), "N1600_eta0.740": (1600, 0.740, 40.0), "N1600_eta0.780": (1600, 0.780, 40.0)}


def mean_se(v):
    v = np.asarray(v, float); return v.mean(), v.std(ddof=1) / math.sqrt(len(v)), len(v)


def fig3(out):
    p6 = committed("p6"); st = {s[0]: s for s in STARTS}; rows, summ = [], []
    print("\n## Fig. 3: the solid, Z per compartment (P6's printed tables; PILOT)\n")
    parts = re.split(r"^## (N\d+_eta[\d.]+)\s*$", p6, flags=re.M)
    secs = dict(zip(parts[1::2], parts[2::2]))
    if set(secs) != set(CELLS): stop(f"P6: cells {sorted(secs)}")
    for name, (N, eta, H) in CELLS.items():
        sec = secs[name]; L0 = N * math.pi / (8 * H * eta)
        cut = re.search(r"split at their largest gap ([\d.]+) \(cut ([\d.]+)\): the lower group (\d+) at ([\d.]+) \(([\d.]+)\), "
                        r"the higher group (\d+) at ([\d.]+) \(([\d.]+)\)", sec)
        for r in md_rows(sec, "| start | seed |"):
            if r[0] not in st: stop(f"P6 {name}: unknown start {r[0]!r}")
            lab, key, build, win, kind = st[r[0]]
            for c in (0, 1):
                Z, ps, o = float(r[2 + c]), r[4 + 2 * c], r[5 + 2 * c]
                orient = math.nan if o == "-" else float(o)
                ocls = "not measured" if o == "-" else ("0 deg" if abs(orient) < 15 else "turned")
                grp = ""
                if cut and kind == "compressed": grp = "lower" if Z < float(cut.group(2)) else "higher"
                rows.append(dict(cell=name, N=N, eta=eta, H=H, L0=L0, radius=0.5, start=key, start_label=lab, kind=kind, build=build, window=win,
                                 seed=int(r[1]), compartment=c, Z=Z, abs_psi=(math.nan if ps == "-" else float(ps)), orientation_deg=orient,
                                 orientation=ocls, group=grp))
        R = pd.DataFrame([x for x in rows if x["cell"] == name])
        for r in md_rows(sec, "| start | compartments |"):
            lab, key = st[r[0]][0], st[r[0]][1]
            m, s, n = mean_se(R[R.start == key].Z)
            pm_ = re.fullmatch(r"([\d.]+) \(([\d.]+)\)", r[2])
            ok = n == int(r[1]) and abs(m - float(pm_.group(1))) <= 1e-4 and abs(s - float(pm_.group(2))) <= 1e-4
            check(f"{name} {key}: mean Z_c (SE) {m:.4f} ({s:.4f}), n {n} = P6's {r[2]}, n {r[1]}", ok)
            summ.append(dict(cell=name, row="start", start=key, start_label=lab, n=n, mean=m, se=s, printed=r[2], cut=math.nan))
        if cut:
            cv = float(cut.group(2)); C = R[R.kind == "compressed"]
            for g, (gn, gm, gs) in (("lower", cut.group(3, 4, 5)), ("higher", cut.group(6, 7, 8))):
                m, s, n = mean_se(C[C.group == g].Z)
                ok = n == int(gn) and abs(m - float(gm)) <= 1e-4 and abs(s - float(gs)) <= 1e-4
                check(f"{name} compressed {g} group (cut {cv}): {n} at {m:.4f} ({s:.4f}) = P6's {gn} at {gm} ({gs})", ok)
                summ.append(dict(cell=name, row=f"compressed {g} group", start="compressed", start_label="all compressed starts", n=n, mean=m,
                                 se=s, printed=f"{gn} at {gm} ({gs}); gap {cut.group(1)}", cut=cv))
        else:
            print(f"- {name}: P6 printed no split of the compressed compartments; none is made here")
    pd.DataFrame(rows).to_csv(os.path.join(out, "fig3_solid_Z_per_compartment.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(out, "fig3_solid_groups.csv"), index=False)
    print(f"\ncompartments: {len(rows)}; orientation measured for {sum(1 for x in rows if x['orientation'] != 'not measured')}")


# ---------------------------------------------------------------- Fig. 4: the KOA rate table and the window speed-up
def fig4(out):
    rt = committed("rates"); rows = []
    print("\n## Fig. 4: the KOA rate table (Gold node) and the window speed-up\n")
    hdr = re.search(r"^# K2-2 \(job (\d+); ([^,]+), ([^)]*?) per decision 22\)", rt, re.M)
    if not hdr: stop("rates: the header line")
    etas = {"pi/8": math.pi / 8, "0.70": 0.70, "0.78": 0.78}
    for r in md_rows(rt, "| N | eta | events | run [s] |"):
        N, eta = int(r[0]), etas[r[1]]; H = 10 * math.sqrt(N / 100.0)
        rows.append(dict(N=N, eta_label=r[1], eta=eta, H=H, L0=N * math.pi / (8 * H * eta), radius=0.5, M=300, events=int(r[2]),
                         run_s=float(r[3]), engine_s=float(r[4]), driver_share=float(r[5]), events_per_s_run=float(r[6]),
                         events_per_s_engine=float(r[7]), clean=int(r[8]), job=hdr.group(1), node=hdr.group(2), cpu=hdr.group(3),
                         build="a643bc7 (target koa)", engine="gen3"))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(out, "fig4_koa_rates_gold.csv"), index=False)
    check("rate table: 12 rows (N 100, 400, 900, 1600 x 3 densities), every run clean", len(R) == 12 and (R.clean == 1).all())
    m = re.search(r"E5/Gold = ([\d.]+) \(Gold ([\d.]+)x faster\)", rt)
    print(f"node {hdr.group(2)} ({hdr.group(3)}), job {hdr.group(1)}; information printed with it: E5/Gold = {m.group(1) if m else '?'} at N 400, eta 0.70")
    w = []
    for r in md_rows(committed("d19"), "| N | time scale | engine | frames per second |"):
        N = int(r[0]); f = math.sqrt(N / 400.0)
        w.append(dict(N=N, time_scale=int(r[1]), engine=r[2], fps=float(r[3]), sigma_time_per_wall_s=float(r[4]),
                      events_per_wall_s=(math.nan if r[5] == "-" else float(r[5])), frames=int(r[6]), steps=int(r[7]),
                      clean=(math.nan if r[8] == "-" else int(r[8])), eta=0.31, H=20.0 * f, L0=25.34 * f, radius=0.5, wall_mass_factor=50,
                      build="7f6adca (target mac-O3-e0pre)"))
    W = pd.DataFrame(w)
    check("window table: 8 rows, gen3 clean", len(W) == 8 and (W[W.engine == "gen3"].clean == 1).all())
    sp = []
    for (N, ts), g in W.groupby(["N", "time_scale"]):
        a, b = g[g.engine == "gen3"].iloc[0], g[g.engine == "gen2"].iloc[0]
        sp.append(dict(N=N, time_scale=ts, speedup_sigma_time=a.sigma_time_per_wall_s / b.sigma_time_per_wall_s, speedup_fps=a.fps / b.fps,
                       gen3_sigma_time_per_wall_s=a.sigma_time_per_wall_s, gen2_sigma_time_per_wall_s=b.sigma_time_per_wall_s,
                       gen3_frames=a.frames, gen2_frames=b.frames, eta=0.31, H=a.H, L0=a.L0, radius=0.5, wall_mass_factor=50, build=a.build))
    S = pd.DataFrame(sp); S.to_csv(os.path.join(out, "fig4_window_speedup.csv"), index=False)
    W.to_csv(os.path.join(out, "fig4_window_table.csv"), index=False)
    print("\n| N | time scale | gen3 / gen2, sigma-time per wall second | gen3 / gen2, frames per second | frames gen3 / gen2 |\n|---|---|---|---|---|")
    for _, s in S.iterrows():
        print(f"| {s.N} | {s.time_scale} | {s.speedup_sigma_time:.3g} | {s.speedup_fps:.3g} | {s.gen3_frames} / {s.gen2_frames} |")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.join(HERE, "tables")); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    print("# M1 (decision 29): the tables of the four meeting figures, printed by experiments_m1_figures_261011/m1_tables_261011.py; "
          "existing data only, no runs\n\n## The committed sources (read with git show)\n\n| name | commit | path | blob |\n|---|---|---|---|")
    for k, (c, p) in SRC.items(): print(f"| {k} | {c[:8]} | {p} | {git('rev-parse', f'{c}:{p}').strip()[:12]} |")
    fig1(a.out); fig2(a.out); fig3(a.out); fig4(a.out)
    print("\n## The tables written\n\n| file | rows | SHA-256 |\n|---|---|---|")
    for f in sorted(os.listdir(a.out)):
        if f.endswith(".csv"): print(f"| {f} | {len(pd.read_csv(os.path.join(a.out, f)))} | {sha256(os.path.join(a.out, f))} |")
    print(f"\nCHECKS: {'all PASS' if not FAILS else 'FAIL: ' + '; '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
