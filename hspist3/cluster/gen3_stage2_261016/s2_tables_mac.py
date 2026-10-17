#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 58; registration S2, notes 261012 sec. 4.9): THE S2 TABLES ON THE MAC -- the registration's
estimators, all pre-specified there, from the fetched tree (hspist3/experiments_gen3_stage2_koa_261016/S2; fetch_s2_koa.sh). Reads only the
merged tables (merged/stage2/gen3/<cell id>/m_<M>/: red_nu.csv, nud_runs.csv, runs_health.csv, acf_runs.npz from the stage 1 reducers;
p2_runs.csv, p2_frames.csv from p2_reduce.py); writes into a NEW folder --out. The cell id S2_N<N>_eta<eta>_k<k>_pin<0|turned> (gen_s2).
  1. E1 per (cell, H, mass): stage 1's per-mass E1 (s1_verdict_261012.per_mass's steps, as p2_tables_mac.per_mass: the used trajectories
     gen3 clean and fit_ok, nu_d-bar, tau_r of the item-4a fit of their seed-averaged ACF, nu_0, Q, the delete-one jackknife SE, c = nu_0 / x
     with x = K(alpha) / (2 pi (L0 - 1.025)), alpha = M / N); the group rule: more than 2 of the 10 trajectories lost -> the mass is out; the
     low-Q flag: Q < 5.
  2. E1 per (cell, H) over the masses: s1_verdict_261012.e1 (the inverse-variance mean of 5 unflagged masses with an SE; the SE scaled by
     sqrt(chi2/dof) when that exceeds 1; no value with fewer than 5).
  3. The relaxed strain e per (cell, H): per run the mean over its record's snapshots (p2_frames.csv, phase 1) of both compartments' strain
     (p2_lattice.config_strain, as p2_reduce.py), then the mean over the (cell, H)'s runs, with its SE (sd / sqrt(n)).
  4. c0 per cell (N, eta, pin): a weighted straight line ln c = ln c0 + b e through its three (cell, H) points (c the (cell, H) E1, e its
     relaxed strain; the effective variance (SE_c / c)^2 + b^2 SE_e^2, iterated); c0 = exp(ln c0) with SE c0 x SE(ln c0); the slope
     b = d ln c / d e with its SE; both SE multiplied by sqrt(chi2/dof) when that exceeds 1 (dof = points - 2).
  5. The orientation: a (cell, H) whose end orientation (the line direction of the last snapshot, config_strain: x for pin 0, y for turned)
     differs from the pin in any compartment of any run is REPORTED, with its count, and KEPT.
  6. The compartment flag: a compartment whose record-mean strain lies more than 1 % (0.01) from the median of its (cell, H)'s compartments
     is FLAGGED; E1 is given with and without its run (per mass, per (cell, H), and c0 per cell).
  7. The outcomes: c0(eta) for N 900 and 1600 with the slopes (table, figure); c0(1600) - c0(900) per eta with z = difference / sqrt(SE^2 + SE^2);
     VERDICT I (isotropy) at eta 0.76 and 0.78 (N 900): z = (c0 pin 0 - c0 turned) / sqrt(SE_0^2 + SE_turned^2); both |z| < 3.5 ->
     ISOTROPY-CONSISTENT, else INCONSISTENT (a missing c0: NO VERDICT); the diagnostics per cell: psi6, Z, sd/sd_E1, nu_d/nu_pred. The
     comparison with the literature elastic constants is descriptive and not computed here (the registration: no verdict).
##CHRIS 2026-10-16 (the independent review before Report D58-1; notes 261012 sec. 4.7.143). Conventions the registration does not fix,
  each printed in s2_tables.md, for the plan author to keep or change:
  - "lost" (the group rule) counts the trajectories of the 10-seed design that are not usable (missing, not gen3 clean, fit failed); a run
    left out for its compartment flag is not lost, so the group rule decides the same in both versions (with and without the flagged runs);
    the 'without' E1 is stage 1's per-mass E1 on the folder less that run (nu_ref of the seed-averaged fit from the remaining runs).
  - c0 "through the three heights": a cell with another number of heights with an E1 (one mass out removes a height: stage 1's e1
    needs 5 masses, S2 has 5) gets its line and c0 from those heights, marked NOT AS REGISTERED. Verdict I: an eta whose two c0 are
    both as registered and give |z| >= 3.5 makes it INCONSISTENT (the registration's "else"), whatever the other eta shows; otherwise
    a missing c0 (or SE) gives NO VERDICT, a c0 not as registered NO VERDICT AS REGISTERED (its z shown for information), and both
    registered |z| < 3.5 ISOTROPY-CONSISTENT.
  - e per (cell, H): equal weight per run (each run's record mean over its snapshots; a snapshot whose compartment has no strain,
    fewer than 10 ordered disks, is skipped in that mean and counted); a run with a compartment that has no strain in any snapshot
    has no e; those compartments and runs are counted.
  - the diagnostics per cell are over all its runs (no exclusion); Verdict I is also shown without the flagged runs (information).
  - H and L0 of a (cell, H) are taken from any of its runs (one value per (cell, H)), else from the S2 task list.
  - a folder whose stage 1 tables are missing or unreadable does not stop the script: it is listed and its mass is out (lost 10); a
    folder without p2_runs.csv or p2_frames.csv is listed: its runs are absent from the run tables, or have no strain, while its
    stage 1 tables still give its per-mass E1. --out is made only after every table has been read.
Writes s2_masses.csv, s2_heights.csv, s2_cells.csv, s2_runs.csv, s2_flags.csv, s2_tables.md, fig_s2_c_vs_e.png, fig_s2_c0_eta.png.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/s2_tables_mac.py --data <fetched S2 folder> --out <new folder>
"""
import glob, math, os, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "cluster", "gen3_testR_261011"), os.path.join(HS, "validation"), HS):
    sys.path.insert(0, p)
import tests_20260913 as T
import s1_verdict_261012 as V
from p2_tables_mac import seedavg
warnings.filterwarnings("ignore")
NSEED, FLAG_E, ZCRIT, NH = 10, 0.01, 3.5, 3
NEED = ("red_nu.csv", "nud_runs.csv", "runs_health.csv", "acf_runs.npz")
C_PIN = {"0": "#2a6fbb", "turned": "#d9822b"}; C_N = {900: "#1f3b5a", 1600: "#c0392b"}


def parse(cid):
    f = cid.split("_"); return dict(cid=cid, N=int(f[1][1:]), eta=float(f[2][3:]), k=int(f[3][1:]), pin=f[4][3:])


def per_mass(d, exclude=()):
    """stage 1's per-mass E1 (s1_verdict_261012.per_mass's steps) on the folder's runs less `exclude`; a missing or unreadable table: no runs"""
    miss = [f for f in NEED if not os.path.exists(os.path.join(d, f))]
    if miss: return dict(present=0, used=0, fit_failures=0, missing=" ".join(miss))
    try: return per_mass_read(d, exclude)
    except Exception as ex: return dict(present=0, used=0, fit_failures=0, missing=f"unreadable ({type(ex).__name__}: {str(ex)[:80]})")


def per_mass_read(d, exclude):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run"); red = red[~red["run"].isin(exclude)]
    nud = pd.read_csv(os.path.join(d, "nud_runs.csv")); hl = pd.read_csv(os.path.join(d, "runs_health.csv")); z = np.load(os.path.join(d, "acf_runs.npz"))
    clean = set(hl[(hl["gen3_clean"] == 1) & (hl["gen3_records"] == 1)]["run"]); ok = set(nud[nud["fit_ok"] == 1]["run"])
    used = sorted(set(red["run"]) & clean & ok); out = dict(present=len(red), used=len(used), fit_failures=len((set(red["run"]) & clean) - ok), missing="")
    if len(used) < 3: return out
    dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    nd = nud.set_index("run").loc[used, "nu_d"].to_numpy(float); acfs = [z[f"run{r}"] for r in used]; p = seedavg(acfs, dt, nu_ref)
    if p is None: return out
    nu0 = math.sqrt(nd.mean() ** 2 + (1 / (2 * math.pi * p[3])) ** 2); jk = []
    for i in range(len(used)):
        keep = [k for k in range(len(used)) if k != i]; fi = seedavg([acfs[k] for k in keep], dt, nu_ref)
        if fi is not None: jk.append(math.sqrt(nd[keep].mean() ** 2 + (1 / (2 * math.pi * fi[3])) ** 2))
    n = len(jk); se = math.sqrt((n - 1) / n * sum((v - np.mean(jk)) ** 2 for v in jk)) if n >= 3 else math.nan
    out.update(nud=float(nd.mean()), tau_r=float(p[3]), Q=math.pi * nd.mean() * p[3], nu0=nu0, nu0_se=se)
    return out


def line_fit(e, se, c, sc):
    """ln c = a + b e, weighted, the effective variance iterated; returns a, b, their SE (scaled when chi2/dof > 1), chi2, dof"""
    e, se, c, sc = (np.asarray(v, float) for v in (e, se, c, sc)); y, sy = np.log(c), sc / c; n = len(e)
    if n < 2: return dict(a=math.nan, b=math.nan, sa=math.nan, sb=math.nan, chi2=math.nan, dof=n - 2, points=n)
    b = np.polyfit(e, y, 1)[0]
    for _ in range(100):
        w = 1 / (sy ** 2 + b ** 2 * se ** 2); S, Sx, Sy, Sxx, Sxy = w.sum(), (w * e).sum(), (w * y).sum(), (w * e * e).sum(), (w * e * y).sum()
        D = S * Sxx - Sx ** 2; bn = (S * Sxy - Sx * Sy) / D
        if abs(bn - b) < 1e-13: b = bn; break
        b = bn
    w = 1 / (sy ** 2 + b ** 2 * se ** 2); S, Sx, Sy, Sxx, Sxy = w.sum(), (w * e).sum(), (w * y).sum(), (w * e * e).sum(), (w * e * y).sum()
    D = S * Sxx - Sx ** 2; a = (Sxx * Sy - Sx * Sxy) / D; chi2 = float((w * (y - a - b * e) ** 2).sum()); dof = n - 2
    k = max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else 1.0
    return dict(a=a, b=b, sa=math.sqrt(Sxx / D) * k, sb=math.sqrt(S / D) * k, chi2=chi2, dof=dof, points=n)


def fx(v, d=3, sign=False): return "-" if v is None or not np.isfinite(v) else (f"{v:+.{d}f}" if sign else f"{v:.{d}f}")
def pm(c, s, d=3): return "-" if not np.isfinite(c) else (f"{c:.{d}f} ± {s:.{d}f}" if np.isfinite(s) else f"{c:.{d}f} (no SE)")
def hts(n): return f"{n}" if n == NH else f"{n} (NOT AS REGISTERED)"


def main():
    av = sys.argv; data = os.path.abspath(av[av.index("--data") + 1]); out = os.path.abspath(av[av.index("--out") + 1])
    if os.path.exists(out): sys.exit(f"REFUSED: {out} exists (a new folder only)")
    folders = sorted(glob.glob(os.path.join(data, "merged", "stage2", "gen3", "S2_*", "m_*")))
    if not folders: sys.exit(f"STOP: no merged/stage2/gen3/S2_*/m_<M> under {data}")
    # the runs: the relaxed strain per run and compartment, the end orientation, the diagnostics
    runs, missing = [], []
    for d in folders:
        c = parse(d.split(os.sep)[-2]); M = int(d.split("_")[-1]); rel = os.path.relpath(d, data)
        if not os.path.exists(os.path.join(d, "p2_runs.csv")): missing.append((rel, "p2_runs.csv")); continue
        R = pd.read_csv(os.path.join(d, "p2_runs.csv"))
        if os.path.exists(os.path.join(d, "p2_frames.csv")): F = pd.read_csv(os.path.join(d, "p2_frames.csv")); F = F[F["phase"] == 1]
        else: missing.append((rel, "p2_frames.csv")); F = pd.DataFrame(columns=["run", "phase", "e_0", "e_1"])
        for r in R.to_dict("records"):
            f = F[F["run"] == r["run"]]; sd = r.get("seed", math.nan)
            runs.append(dict(c, M=M, alpha=M / c["N"], run=int(r["run"]), seed=int(sd) if pd.notna(sd) else -1, H=r["H"], L0=r["L0"],
                             e_rec_0=float(f["e_0"].astype(float).mean()) if len(f) else math.nan, e_rec_1=float(f["e_1"].astype(float).mean()) if len(f) else math.nan,
                             snaps_without_strain=int(f["e_0"].isna().sum() + f["e_1"].isna().sum()) if len(f) else 0,
                             end_line_0=r.get("end_line_0", math.nan), end_line_1=r.get("end_line_1", math.nan), psi6=(r.get("psi6_rec_0", math.nan) + r.get("psi6_rec_1", math.nan)) / 2,
                             Z=r.get("Z", math.nan), x_sd=r.get("x_sd", math.nan), nu_d=r.get("nu_d", math.nan), nu_pred=r.get("nu_pred", math.nan),
                             error=str(r.get("error", "")) if pd.notna(r.get("error", "")) else ""))
    if not runs: sys.exit(f"STOP: no p2_runs.csv row in any m_<M> folder under {data}")
    Rn = pd.DataFrame(runs); Rn["e_rec"] = (Rn["e_rec_0"] + Rn["e_rec_1"]) / 2
    Rn["hkey"] = Rn["cid"]; Rn["ckey"] = Rn.apply(lambda r: f"N{r['N']}_eta{r['eta']:.2f}_pin{r['pin']}", axis=1)
    # H and L0 per (cell, H): from any of its runs (one value per (cell, H)), else from the S2 task list
    TL, tp = {}, os.path.join(HERE, "tasks_s2_koa_261017.txt")
    if os.path.exists(tp):
        for l in open(tp):
            f = l.split()
            if len(f) >= 7: TL[f[1].split("/")[2]] = (float(f[6]), float(f[5]))
    first = lambda s: float(s.dropna().iloc[0]) if s.notna().any() else math.nan
    geo = {cid: (first(g["H"]), first(g["L0"])) for cid, g in Rn.groupby("cid")}
    def hl0(cid):
        h, l = geo.get(cid, (math.nan, math.nan))
        return (h, l) if np.isfinite(h) and np.isfinite(l) else TL.get(cid, (h, l))
    # 6. the compartment flags: more than 1 % from the median of the (cell, H)'s compartments
    fl = []
    for hk, g in Rn.groupby("hkey"):
        allc = np.concatenate([g["e_rec_0"].to_numpy(float), g["e_rec_1"].to_numpy(float)]); med = float(np.nanmedian(allc)) if np.isfinite(allc).any() else math.nan
        for r in g.to_dict("records"):
            for k in (0, 1):
                if np.isfinite(r[f"e_rec_{k}"]) and abs(r[f"e_rec_{k}"] - med) > FLAG_E:
                    fl.append(dict(cid=hk, M=r["M"], run=r["run"], seed=r["seed"], compartment="left" if k == 0 else "right", e_rec=r[f"e_rec_{k}"], median=med))
    FL = pd.DataFrame(fl, columns=["cid", "M", "run", "seed", "compartment", "e_rec", "median"])
    flagged_runs = {(r["cid"], int(r["M"]), int(r["run"])) for r in FL.to_dict("records")}
    # 1. E1 per (cell, H, mass), with and without the flagged runs; "lost" is the same in both versions (a run left out for its flag is not lost)
    ms = []
    for d in folders:
        c = parse(d.split(os.sep)[-2]); M = int(d.split("_")[-1]); L0 = hl0(c["cid"])[1]
        a = M / c["N"]; x = T.k_root(a) / (2 * math.pi * (L0 - 1.025)); ex = tuple(sorted(r for (cid, m, r) in flagged_runs if cid == c["cid"] and m == M))
        lost = None
        for tag, excl in (("all", ()), ("without flagged", ex)):
            if tag != "all" and not ex: continue
            pmm = per_mass(d, excl)
            if tag == "all":
                lost = NSEED - pmm["used"]
                if pmm["missing"]: missing.append((os.path.relpath(d, data), pmm["missing"]))
            ms.append(dict(c, M=M, alpha=a, L0=L0, version=tag, excluded=len(excl), **pmm, lost=lost, flagged=lost > 2, x=x,
                           c=pmm.get("nu0", math.nan) / x, sig_c=pmm.get("nu0_se", math.nan) / x, lowQ=int(pmm.get("Q", math.inf) < 5)))
    MS = pd.DataFrame(ms)
    os.makedirs(out)   # only now: every table has been read
    FL.to_csv(os.path.join(out, "s2_flags.csv"), index=False, float_format="%.8g"); MS.to_csv(os.path.join(out, "s2_masses.csv"), index=False, float_format="%.8g")
    # 2-3. E1 per (cell, H) over the masses; the relaxed e; the orientation
    hs = []
    for hk, g in Rn.groupby("hkey", sort=False):
        c = parse(hk); want = 0.0 if c["pin"] == "0" else 90.0
        lines = g[["end_line_0", "end_line_1"]].to_numpy(float).ravel(); off = int(np.sum(np.isfinite(lines) & (lines != want)))
        n = int(g["e_rec"].notna().sum()); comp = g[["e_rec_0", "e_rec_1"]].to_numpy(float)
        base = dict(c, H=hl0(hk)[0], L0=hl0(hk)[1], runs=len(g), runs_error=int((g["error"] != "").sum()), runs_without_e=int(g["e_rec"].isna().sum()),
                    compartments_without_strain=int(np.isnan(comp).sum()), snapshots_without_strain=int(g["snaps_without_strain"].sum()), e=float(g["e_rec"].mean()), se_e=float(g["e_rec"].std(ddof=1) / math.sqrt(n)) if n > 1 else math.nan,
                    end_off_pin=off, end_compartments=int(np.isfinite(lines).sum()), flagged_compartments=int((FL["cid"] == hk).sum()) if len(FL) else 0,
                    psi6=float(g["psi6"].mean()), Z=float(g["Z"].mean()), nud_over_nup=float((g["nu_d"] / g["nu_pred"]).median()))
        for tag in ("all", "without flagged"):
            m = MS[(MS["cid"] == hk) & (MS["version"] == tag)] if tag == "all" or ((MS["cid"] == hk) & (MS["version"] == tag)).any() else MS[(MS["cid"] == hk) & (MS["version"] == "all")]
            if tag != "all":   # the masses without a flagged run keep their 'all' row
                m = pd.concat([m, MS[(MS["cid"] == hk) & (MS["version"] == "all") & ~MS["M"].isin(m["M"])]])
                gx = g[~g.apply(lambda r: (r["cid"], int(r["M"]), int(r["run"])) in flagged_runs, axis=1)]
                e_x, se_x = float(gx["e_rec"].mean()), float(gx["e_rec"].std(ddof=1) / math.sqrt(max(1, gx["e_rec"].notna().sum())))
            else: e_x, se_x = base["e"], base["se_e"]
            rows = [dict(alpha=r["alpha"], flagged=bool(r["flagged"]), c=r["c"], sig_c=r["sig_c"], Q=r.get("Q", math.nan), M=int(r["M"])) for r in m.to_dict("records")]
            e = V.e1(rows) if rows else None
            usable = sum(1 for r in rows if r["alpha"] >= 1 and not r["flagged"] and np.isfinite(r["c"]) and np.isfinite(r["sig_c"]) and r["sig_c"] > 0)
            hs.append(dict(base, version=tag, e=e_x, se_e=se_x, E1_c=e["cs"] if e else math.nan, E1_sig=e["sig"] if e else math.nan,
                           chi2=e["chi2"] if e else math.nan, masses=e["n"] if e else 0, masses_usable=usable, lowQ=int(e["lowQ"]) if e else 0))
    HSd = pd.DataFrame(hs); HSd.to_csv(os.path.join(out, "s2_heights.csv"), index=False, float_format="%.8g")
    # sd / sd_E1 per run (sd_E1 = sqrt(L_eff / (2 H rho c^2)), c the (cell, H, mass) E1)
    cm = {(r["cid"], int(r["M"])): r["c"] for r in MS[MS["version"] == "all"].to_dict("records")}
    Rn["sd_E1"] = [math.sqrt((r["L0"] - 1.025) / (2 * r["H"] * (4 * r["eta"] / math.pi) * cm[(r["cid"], int(r["M"]))] ** 2))
                   if np.isfinite(cm.get((r["cid"], int(r["M"])), math.nan)) and cm[(r["cid"], int(r["M"]))] > 0 else math.nan for r in Rn.to_dict("records")]
    Rn["sd_over_E1"] = Rn["x_sd"] / Rn["sd_E1"]; Rn.to_csv(os.path.join(out, "s2_runs.csv"), index=False, float_format="%.8g")
    # 4. c0 per cell (as registered only with all three heights)
    cs = []
    for ck, g in HSd.groupby(HSd.apply(lambda r: (r["N"], r["eta"], r["pin"]), axis=1), sort=False):
        N, eta, pin = ck
        for tag in ("all", "without flagged"):
            h = g[g["version"] == tag] if (g["version"] == tag).any() else g[g["version"] == "all"]
            h = h[np.isfinite(h["E1_c"]) & np.isfinite(h["E1_sig"]) & np.isfinite(h["e"]) & np.isfinite(h["se_e"])]
            f = line_fit(h["e"], h["se_e"], h["E1_c"], h["E1_sig"]); c0 = math.exp(f["a"]) if np.isfinite(f["a"]) else math.nan
            rr = Rn[Rn["ckey"] == f"N{N}_eta{eta:.2f}_pin{pin}"]; ga = g[g["version"] == "all"]
            cs.append(dict(N=N, eta=eta, pin=pin, version=tag, points=f["points"], as_registered=int(f["points"] == NH), c0=c0, se_c0=c0 * f["sa"] if np.isfinite(c0) else math.nan,
                           slope=f["b"], se_slope=f["sb"], chi2=f["chi2"], dof=f["dof"], e_min=float(h["e"].min()) if len(h) else math.nan,
                           e_max=float(h["e"].max()) if len(h) else math.nan, psi6=float(rr["psi6"].mean()), Z=float(rr["Z"].mean()),
                           sd_over_E1=float(rr["sd_over_E1"].median()), nud_over_nup=float((rr["nu_d"] / rr["nu_pred"]).median()),
                           lowQ_heights=int(ga["lowQ"].sum()), end_off_pin=int(ga["end_off_pin"].sum()), flagged=int(ga["flagged_compartments"].sum())))
    CS = pd.DataFrame(cs); CS.to_csv(os.path.join(out, "s2_cells.csv"), index=False, float_format="%.8g")
    # 7. the outcomes
    A = CS[CS["version"] == "all"]; W = CS[CS["version"] == "without flagged"]; HA = HSd[HSd["version"] == "all"]
    L = ["# The S2 tables (s2_tables_mac.py; registration sec. 4.9)", "", f"- data: {data}",
         f"- runs: {len(Rn)} (with a p2_reduce error: {int((Rn['error'] != '').sum())}); (cell, H): {Rn['hkey'].nunique()}; cells: {A.shape[0]}; m_<M> folders: {len(folders)}, "
         f"with a missing table: {len({m[0] for m in missing})}",
         f"- flagged compartments: {len(FL)}; (cell, H) with an end orientation off the pin: {int((HA['end_off_pin'] > 0).sum())}; compartments with no strain "
         f"in any snapshot: {int(HA['compartments_without_strain'].sum())}; runs without e: {int(HA['runs_without_e'].sum())}; snapshot compartments without a strain "
         f"(skipped in the record means): {int(HA['snapshots_without_strain'].sum())}",
         "- conventions (not fixed by the registration; the docstring): a run left out for its flag is not 'lost'; c0 as registered needs exactly three heights "
         "(else NOT AS REGISTERED); Verdict I is INCONSISTENT when an eta with both c0 as registered gives |z| >= 3.5, else NO VERDICT (a c0 missing) or NO VERDICT "
         "AS REGISTERED (a c0 not as registered) before CONSISTENT; e per (cell, H) weighs each run equally (its record mean over the snapshots, both compartments; "
         "a snapshot without a strain is skipped; a run with a compartment with no strain in any snapshot has no e); the diagnostics are over all runs of the cell.", ""]
    L += ["## c0 per cell (the straight line of ln c on e through the heights; c0 at e = 0)", "",
          "| N | eta | pin | heights | c0 ± SE | d ln c / d e ± SE | chi2 / dof | e range [%] | psi6 | Z | sd / sd_E1 | nu_d / nu_pred | low-Q heights | end off the pin | flagged compartments |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in A.sort_values(["N", "pin", "eta"]).to_dict("records"):
        L.append(f"| {r['N']} | {r['eta']:.2f} | {r['pin']} | {hts(r['points'])} | {pm(r['c0'], r['se_c0'])} | {pm(r['slope'], r['se_slope'], 2)} | {fx(r['chi2'], 2)} / {r['dof']} | "
                 f"{fx(100 * r['e_min'], 3, True)} to {fx(100 * r['e_max'], 3, True)} | {fx(r['psi6'])} | {fx(r['Z'], 2)} | {fx(r['sd_over_E1'], 2)} | {fx(r['nud_over_nup'])} | "
                 f"{r['lowQ_heights']} | {r['end_off_pin']} | {r['flagged']} |")
    if len(FL):
        L += ["", "c0 without the flagged runs (the cells that have one):", "", "| N | eta | pin | c0 ± SE (all) | heights | c0 ± SE (without) | heights |", "|---|---|---|---|---|---|---|"]
        for r in A.to_dict("records"):
            w = W[(W["N"] == r["N"]) & np.isclose(W["eta"], r["eta"]) & (W["pin"] == r["pin"])]
            if r["flagged"] and len(w): L.append(f"| {r['N']} | {r['eta']:.2f} | {r['pin']} | {pm(r['c0'], r['se_c0'])} | {hts(r['points'])} | {pm(w['c0'].iloc[0], w['se_c0'].iloc[0])} | {hts(w['points'].iloc[0])} |")
    L += ["", "## c0(1600) - c0(900) per eta (turned)", "", "| eta | c0 N 900 | c0 N 1600 | difference ± SE | z |", "|---|---|---|---|---|"]
    for eta in sorted(A["eta"].unique()):
        a9 = A[(A["N"] == 900) & np.isclose(A["eta"], eta) & (A["pin"] == "turned")]; a16 = A[(A["N"] == 1600) & np.isclose(A["eta"], eta) & (A["pin"] == "turned")]
        if not (len(a9) and len(a16)): continue
        d = a16["c0"].iloc[0] - a9["c0"].iloc[0]; s = math.hypot(a16["se_c0"].iloc[0], a9["se_c0"].iloc[0])
        mark = "" if a9["points"].iloc[0] == NH and a16["points"].iloc[0] == NH else " (NOT AS REGISTERED: a c0 not from exactly three heights)"
        L.append(f"| {eta:.2f} | {pm(a9['c0'].iloc[0], a9['se_c0'].iloc[0])} | {pm(a16['c0'].iloc[0], a16['se_c0'].iloc[0])} | {pm(d, s)} | {fx(d / s if s > 0 else math.nan, 2)}{mark} |")

    def verdict_i(X):
        """the z at 0.76 and 0.78 and the verdict, from the cell table X"""
        rows, zs, regs = [], [], []
        for eta in (0.76, 0.78):
            p0 = X[(X["N"] == 900) & np.isclose(X["eta"], eta) & (X["pin"] == "0")]; pt = X[(X["N"] == 900) & np.isclose(X["eta"], eta) & (X["pin"] == "turned")]
            if not (len(p0) and len(pt)) or not (np.isfinite(p0["c0"].iloc[0]) and np.isfinite(pt["c0"].iloc[0])):
                zs.append(math.nan); regs.append(False); rows.append(f"| {eta:.2f} | - | - | - | - |"); continue
            d = p0["c0"].iloc[0] - pt["c0"].iloc[0]; s = math.hypot(p0["se_c0"].iloc[0], pt["se_c0"].iloc[0]); z = d / s if s > 0 else math.nan; zs.append(z)
            reg = p0["points"].iloc[0] == NH and pt["points"].iloc[0] == NH; regs.append(reg)
            rows.append(f"| {eta:.2f} | {pm(p0['c0'].iloc[0], p0['se_c0'].iloc[0])} [{hts(p0['points'].iloc[0])}] | {pm(pt['c0'].iloc[0], pt['se_c0'].iloc[0])} [{hts(pt['points'].iloc[0])}] | "
                        f"{pm(d, s)} | {fx(z, 2)} |")
        # the registration's "else INCONSISTENT" first: one eta with both c0 as registered and |z| >= 3.5 decides it
        v = ("INCONSISTENT" if any(r and np.isfinite(z) and abs(z) >= ZCRIT for z, r in zip(zs, regs)) else
             "NO VERDICT (a c0 or its SE is missing)" if not all(np.isfinite(zs)) else
             "NO VERDICT AS REGISTERED (a c0 not from exactly three heights; the z is shown for information)" if not all(regs) else
             "ISOTROPY-CONSISTENT")
        return rows, v
    rows, verdict = verdict_i(A)
    L += ["", "## VERDICT I (isotropy), N 900: z = (c0 pin 0 - c0 turned) / sqrt(SE_0^2 + SE_turned^2); both |z| < 3.5 -> ISOTROPY-CONSISTENT", "",
          "| eta | c0 pin 0 [heights] | c0 turned [heights] | difference ± SE | z |", "|---|---|---|---|---|"] + rows + ["", f"**VERDICT I: {verdict}**", ""]
    fc = A[(A["N"] == 900) & A["eta"].apply(lambda v: any(np.isclose(v, (0.76, 0.78))))]
    if int(fc["flagged"].sum()):
        rows_w, verdict_w = verdict_i(W)
        L += ["Verdict I without the flagged runs (information; a 0.76 or 0.78 cell has a flagged compartment):", "",
              "| eta | c0 pin 0 [heights] | c0 turned [heights] | difference ± SE | z |", "|---|---|---|---|---|"] + rows_w + ["", f"without the flagged runs: {verdict_w}", ""]
    L += ["## The literature elastic constants of the hard-disk crystal (descriptive, no verdict)", "",
          "Not computed by this script: the comparison is written in the S2 report from the literature values.", ""]
    L += ["## Per (cell, H): E1 over the masses and the relaxed strain", "",
          "e per (cell, H): the mean over its runs of each run's record mean (both compartments; equal weight per run; a snapshot without a strain skipped).", "",
          "| cell id | runs | runs with an error | masses in E1 | E1: c ± SE | chi2 | low Q | e [%] ± SE | runs without e | compartments with no strain in any snapshot | snapshot compartments without a strain | end off the pin (compartments) | flagged compartments |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in HA.sort_values(["N", "eta", "pin", "k"]).to_dict("records"):
        mi = f"{r['masses']}" if r["masses"] else f"0 ({r['masses_usable']} usable, fewer than 5: no E1)"
        L.append(f"| {r['cid']} | {r['runs']} | {r['runs_error']} | {mi} | {pm(r['E1_c'], r['E1_sig'])} | {fx(r['chi2'], 2)} | {r['lowQ']} | {pm(100 * r['e'], 100 * r['se_e'])} | "
                 f"{r['runs_without_e']} | {r['compartments_without_strain']} | {r['snapshots_without_strain']} | {r['end_off_pin']} of {r['end_compartments']} | {r['flagged_compartments']} |")
    L += ["", "## Per (cell, H, mass): E1", "", "| cell id | alpha | used / present | lost | out (> 2 lost) | Q | low Q | c ± SE |", "|---|---|---|---|---|---|---|---|"]
    for r in MS[MS["version"] == "all"].sort_values(["N", "eta", "pin", "k", "alpha"]).to_dict("records"):
        L.append(f"| {r['cid']} | {r['alpha']:g} | {r['used']} / {r['present']} | {r['lost']} | {'OUT' if r['flagged'] else '-'} | {fx(r.get('Q', math.nan), 1)} | {r['lowQ']} | {pm(r['c'], r['sig_c'])} |")
    if len(FL):
        L += ["", "## The flagged compartments (more than 1 % from their (cell, H) median)", "", "| cell id | M | run | seed | compartment | e [%] | median [%] |", "|---|---|---|---|---|---|---|"]
        for r in FL.to_dict("records"):
            L.append(f"| {r['cid']} | {r['M']} | {r['run']} | {r['seed']} | {r['compartment']} | {100 * r['e_rec']:+.3f} | {100 * r['median']:+.3f} |")
        L += ["", "E1 with and without the flagged runs (the (cell, H) and the masses that have one; 'lost' is the same in both versions):", "",
              "| cell id | over | all: c ± SE (used) | without: c ± SE (used) |", "|---|---|---|---|"]
        for hk in sorted({r["cid"] for r in FL.to_dict("records")}):
            for M in sorted({int(r["M"]) for r in FL.to_dict("records") if r["cid"] == hk}):
                ra = MS[(MS["cid"] == hk) & (MS["M"] == M) & (MS["version"] == "all")]; rw = MS[(MS["cid"] == hk) & (MS["M"] == M) & (MS["version"] == "without flagged")]
                if len(ra) and len(rw):
                    L.append(f"| {hk} | M {M} | {pm(ra['c'].iloc[0], ra['sig_c'].iloc[0])} ({ra['used'].iloc[0]}{', OUT' if ra['flagged'].iloc[0] else ''}) | "
                             f"{pm(rw['c'].iloc[0], rw['sig_c'].iloc[0])} ({rw['used'].iloc[0]}{', OUT' if rw['flagged'].iloc[0] else ''}) |")
            ha = HSd[(HSd["cid"] == hk) & (HSd["version"] == "all")]; hw = HSd[(HSd["cid"] == hk) & (HSd["version"] == "without flagged")]
            if len(ha) and len(hw):
                L.append(f"| {hk} | the masses (E1), e [%] | {pm(ha['E1_c'].iloc[0], ha['E1_sig'].iloc[0])} ({ha['masses'].iloc[0]} masses), e {pm(100 * ha['e'].iloc[0], 100 * ha['se_e'].iloc[0])} | "
                         f"{pm(hw['E1_c'].iloc[0], hw['E1_sig'].iloc[0])} ({hw['masses'].iloc[0]} masses), e {pm(100 * hw['e'].iloc[0], 100 * hw['se_e'].iloc[0])} |")
    if missing:
        L += ["", "## Folders with a missing or unreadable table (stage 1 tables: the mass is out; p2_runs.csv: its runs are absent from the run tables; "
              "p2_frames.csv: its runs have no strain)", "", "| folder | missing |", "|---|---|"] + [f"| {f} | {m} |" for f, m in missing]
    open(os.path.join(out, "s2_tables.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    # the figures
    cells = list(A.sort_values(["N", "pin", "eta"])[["N", "eta", "pin"]].itertuples(index=False, name=None))
    nc = 4; nr = max(1, math.ceil(len(cells) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.4 * nr), squeeze=False)
    for i, (N, eta, pin) in enumerate(cells):
        ax = axs[i // nc, i % nc]; h = HSd[(HSd["N"] == N) & np.isclose(HSd["eta"], eta) & (HSd["pin"] == pin) & (HSd["version"] == "all")]
        ax.errorbar(100 * h["e"], h["E1_c"], xerr=100 * h["se_e"], yerr=h["E1_sig"], fmt="o", color=C_PIN[pin], capsize=3, ms=5)
        r = A[(A["N"] == N) & np.isclose(A["eta"], eta) & (A["pin"] == pin)].iloc[0]
        if np.isfinite(r["c0"]):
            xs = np.linspace(min(-0.8, 100 * r["e_min"] - 0.2), max(0.8, 100 * r["e_max"] + 0.2), 30)
            ax.plot(xs, r["c0"] * np.exp(r["slope"] * xs / 100), color="0.4", lw=0.9, ls="--"); ax.errorbar([0], [r["c0"]], yerr=[r["se_c0"]], fmt="s", color="k", ms=4, capsize=3)
        ax.axvline(0, color="0.7", lw=0.6, ls=":"); ax.grid(alpha=0.25); ax.tick_params(labelsize=7)
        ax.set_title(f"N {N}, eta {eta:.2f}, pin {pin}: c0 {fx(r['c0'])} ± {fx(r['se_c0'])}" + ("" if r["points"] == NH else f" ({r['points']} heights)"), fontsize=8)
        ax.set_xlabel("relaxed e [%]", fontsize=7); ax.set_ylabel("c (E1)", fontsize=7)
    for j in range(len(cells), nr * nc): axs[j // nc, j % nc].axis("off")
    fig.suptitle("S2: c (E1 over the masses) against the relaxed strain per (cell, H); the line in ln c and c0 at e = 0 (black square)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig_s2_c_vs_e.png"), dpi=110); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5.2))
    for N in (900, 1600):
        a = A[(A["N"] == N) & (A["pin"] == "turned")].sort_values("eta")
        if len(a): ax.errorbar(a["eta"] + (0.002 if N == 1600 else 0), a["c0"], yerr=a["se_c0"], fmt="o-", color=C_N[N], capsize=3, ms=5, lw=1, label=f"N {N}, turned")
    a = A[(A["N"] == 900) & (A["pin"] == "0")].sort_values("eta")
    if len(a): ax.errorbar(a["eta"] - 0.002, a["c0"], yerr=a["se_c0"], fmt="s", color=C_PIN["0"], mfc="white", capsize=3, ms=6, label="N 900, pin 0 (isotropy check)")
    nr_ = A[(A["points"] != NH) & np.isfinite(A["c0"])]   # a c0 not from exactly three heights: marked
    if len(nr_): ax.plot(nr_["eta"] + np.where(nr_["N"] == 1600, 0.002, np.where(nr_["pin"] == "0", -0.002, 0.0)), nr_["c0"], "x", color="0.3", ms=11, mew=1.5, label="NOT AS REGISTERED (not three heights)")
    ax.set_xlabel("eta", fontsize=9); ax.set_ylabel("c0: c at zero relaxed strain", fontsize=9); ax.grid(alpha=0.25); ax.legend(fontsize=8)
    ax.set_title(f"S2: c0(eta) for N 900 and 1600; VERDICT I: {verdict}", fontsize=10); fig.savefig(os.path.join(out, "fig_s2_c0_eta.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
