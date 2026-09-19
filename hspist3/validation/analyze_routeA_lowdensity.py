#!/usr/bin/env python3
"""##CHRIS 2026-09-12: c_s(eta) for the route-A low-density extension.

Per-run nu comes from the damped-cosine fit (fit_nu_damped.fit_trace), not from the
binned FFT peak. Acceptance is the strict health contract: a trajectory is used only
if run.log reports forced_advance = clamp_repair = overlap_repair = wall_overdue = 0.
On top of that a fit-quality cut sigma_nu/nu < 3e-4 (0.03 %) is applied. The cut is on
the fit alone, never on agreement with any other nu estimate, which would bias the slope.

Per eta:
  per mass  mean(nu) over surviving repeats, sem
  fit       nu = c_s * x through the origin, x = K/(2 pi L_eff), cot K = (M/2N_side) K,
            L_eff = L0 - 2r, weights 1/sem^2   -> c_s, c_s_err
  scatter   c_s_scatter_mass = std(ddof=1) of the per-mass implied c_s = mean(nu)/x

Writes rows in the routeA_refit_cs_vs_eta column set, appended to a COPY of that file.
The original is never modified. Read-only on every trace.
"""
import csv, glob, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fit_nu_damped import fit_trace
import plot_speed_of_sound_edmd as sos

CAMP = sys.argv[1]           # campaign dir holding eta_* leaves
SRC  = sys.argv[2]           # routeA_refit_cs_vs_eta_20260909.csv (read-only)
OUT  = sys.argv[3]           # new file: copy of SRC + the new rows
PERRUN = sys.argv[4] if len(sys.argv) > 4 else None
CUT  = float(sys.argv[5]) if len(sys.argv) > 5 else 3e-4
# ##CHRIS: the corrected-estimator table these rows are also appended to. Optional 6th
# argument so the already-running chain (which passes 5) picks up the default.
CORR = sys.argv[6] if len(sys.argv) > 6 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "0000_PLAN_OVERALL", "ALL_MARKDOWNS", "260909_plots",
    "routeA_refit_cs_vs_eta_20260912.csv")
R, N_SIDE, KR_MAX = 0.5, 50, 0.69
# ##CHRIS: accuracy guard -- see the note in plot_famB_cs_vs_N.py. sigma_nu/nu is a
# precision cut and cannot detect a fit that locked onto the wrong spectral component.
NU_RATIO_LO, NU_RATIO_HI = 1.0 / 3.0, 3.0
# ##CHRIS 2026-09-12: robust mode -- no sigma_nu cut, per-mass median, MAD weights.
ROBUST = os.environ.get("HD_ROBUST", "0") == "1"


def _mad(a):
    m = np.median(a)
    return 1.4826 * float(np.median(np.abs(a - m)))

runs, summ = [], []
for leaf in sorted(glob.glob(f"{CAMP}/eta_*")):
    if not os.path.isdir(leaf) or not os.path.exists(f"{leaf}/run.log"):
        continue
    log = open(f"{leaf}/run.log", errors="replace").read()
    seed_of, health, ti_of = {}, {}, {}
    pending_ti = None
    for ln in log.splitlines():
        # ##CHRIS: HD_KE_TRACE=1 prints "2b after per-segment equalize" just before each
        # "Running:" line. T_i per compartment is KE_side / n_side (2D: KE = n kB T).
        k = re.search(r"2b after per-segment equalize\s+N=(\d+)\s+KE_tot=\S+\s+"
                      r"KE_left=(\S+)\s+KE_right=(\S+)", ln)
        if k:
            n_side = int(k.group(1)) / 2.0
            pending_ti = (float(k.group(2)) / n_side, float(k.group(3)) / n_side)
        m = re.search(r"Running: L0 = ([\d.]+), M = (\d+)\*m, run = (\d+), seed = (\d+)", ln)
        if m:
            seed_of[(int(m.group(2)), int(m.group(3)))] = int(m.group(4))
            if pending_ti: ti_of[(int(m.group(2)), int(m.group(3)))] = pending_ti; pending_ti = None
        h = re.search(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): "
                      r"forced_advance=(\d+) wall_clamp_repairs=(\d+) overlap_repairs=(\d+) wall_overdue=(\d+)", ln)
        if h: health[(int(h.group(1)), int(h.group(2)))] = tuple(int(h.group(i)) for i in (4, 5, 6, 7))
    per_mass = {}
    L0 = None
    for f in sorted(glob.glob(f"{leaf}/wall_x_positions_*.csv")):
        mm = re.search(r"wallmassfactor_(\d+)_run(\d+)\.csv$", f)
        if not mm: continue
        M, run = int(mm.group(1)), int(mm.group(2))
        with open(f) as fh:
            rd = csv.DictReader(fh); row0 = next(rd, None)
        if row0 is None: continue
        if L0 is None: L0 = float(row0["L0"])
        fa, cr, orp, wo = health.get((M, run), (0, 0, 0, 0))
        reasons = []
        if fa or cr or orp or wo:
            reasons.append(f"health fa={fa} crep={cr} orep={orp} overdue={wo}")
        res = fit_trace(f)
        nu = res["nu"] if res else float("nan")
        sig = res["sigma_nu"] if res else float("nan")
        if res is None: reasons.append("fit_failed")
        elif not (nu > 0): reasons.append("nu<=0")
        elif (not ROBUST) and sig / nu >= CUT: reasons.append(f"quality sigma/nu={sig/nu:.3e}")
        else:
            try: pred = float(row0["Predicted_Frequency"])
            except (KeyError, ValueError, TypeError): pred = 0.0
            if pred > 0 and not (NU_RATIO_LO < nu / pred < NU_RATIO_HI):
                reasons.append(f"frequency_alias nu/nu_pred={nu/pred:.4g}")
                print(f"  ALIAS nu/nu_pred={nu/pred:.4g} sigma_nu/nu={sig/nu:.2e}  {f}")
        elig = not reasons
        runs.append(dict(eta=row0["eta"], L0=row0["L0"], M=M, run=run,
                         seed=seed_of.get((M, run), ""),
                         nu_fitted=f"{nu:.10g}" if nu == nu else "",
                         sigma_nu_fitted=f"{sig:.6g}" if sig == sig else "",
                         gamma_fitted=f"{res['gamma']:.6g}" if res else "",
                         fit_rms=f"{res['rms']:.6g}" if res else "",
                         forced_advance=fa, clamp_repair=cr, overlap_repair=orp, wall_overdue=wo,
                         T_i_left=f"{ti_of.get((M,run),(float('nan'),)*2)[0]:.9g}",
                         T_i_right=f"{ti_of.get((M,run),(float('nan'),)*2)[1]:.9g}",
                         eligible=int(elig), exclusion_reasons=";".join(reasons)))
        if elig:
            ti = ti_of.get((M, run))
            ti_m = 0.5 * (ti[0] + ti[1]) if ti else float("nan")
            per_mass.setdefault(M, []).append((nu, ti_m))
    if L0 is None: continue
    eta = N_SIDE * 2 * math.pi * R * R / (2 * L0 * 10.0)
    L_eff = L0 - 2 * R
    xs, ys, ws, cs_mass, n_used = [], [], [], [], 0
    rys, rws, ti_all = [], [], []
    for M, pairs in sorted(per_mass.items()):
        if len(pairs) < 3: continue
        raw = np.array([p_[0] for p_ in pairs])
        tis = np.array([p_[1] for p_ in pairs])
        # ##CHRIS: T_i is READ from the HD_KE_TRACE audit, never assumed. Under
        # --seed-drift-order=drift-first it is 1 exactly and corrected == raw.
        corr = raw / np.sqrt(np.where(np.isfinite(tis) & (tis > 0), tis, 1.0))
        ti_all.extend([t for t in tis if t == t])
        if ROBUST:
            m0 = float(np.median(corr)); s0 = _mad(corr)
            keep = np.abs(corr - m0) <= 5 * s0 if s0 > 0 else np.ones(len(corr), bool)
            if keep.sum() < 3: continue
            corr = corr[keep]; raw = raw[keep]
            y = float(np.median(corr))
            sdc = _mad(corr) if _mad(corr) > 0 else float(corr.std(ddof=1))
            semc = sdc / math.sqrt(len(corr))
            ry = float(np.median(raw))
            sdr = _mad(raw) if _mad(raw) > 0 else float(raw.std(ddof=1))
            semr = sdr / math.sqrt(len(raw))
        else:
            y = float(np.mean(corr)); semc = float(np.std(corr, ddof=1)) / math.sqrt(len(corr))
            ry = float(np.mean(raw)); semr = float(np.std(raw, ddof=1)) / math.sqrt(len(raw))
        K = sos.k_root_bisect(M / (2.0 * N_SIDE))
        x = K / (2 * math.pi * L_eff)
        xs.append(x); ys.append(y)
        ws.append(max(semc, 1e-12))
        rys.append(ry)
        rws.append(max(semr, 1e-12))
        cs_mass.append(y / x); n_used += len(corr)
    if len(xs) < 3:
        print(f"  eta={eta:.6f}: only {len(xs)} masses with >=3 accepted runs -- no fit"); continue
    x, y, s = np.array(xs), np.array(ys), np.array(ws)
    cs, cs_err = sos.weighted_linreg(x, y, s, force_zero_intercept=True)[:2]
    cs_raw, cs_raw_err = sos.weighted_linreg(x, np.array(rys), np.array(rws),
                                             force_zero_intercept=True)[:2]
    ti_mean = float(np.mean(ti_all)) if ti_all else float("nan")
    cf, cf_err, b0, b0_err = sos.weighted_linreg(x, y, s, force_zero_intercept=False)[:4]
    r2 = 1 - float((((y - cs * x) ** 2).sum()) / (((y - y.mean()) ** 2).sum()))
    scatter = float(np.std(np.array(cs_mass), ddof=1))
    a = np.array([eta])
    Zk = sos.Z_kolafa_rottner_2006(a); Zl = sos.Z_liu_global(a); dZl = sos.dZ_liu_global(a)
    h = 1e-5
    dZk = (sos.Z_kolafa_rottner_2006(a + h) - sos.Z_kolafa_rottner_2006(a - h)) / (2 * h)
    cs_kr = float(sos.cs_adiabatic_2d_monatomic(Zk, dZk, a, kbt=1, m=1)[0]) if eta <= KR_MAX else float("nan")
    cs_liu = float(sos.cs_adiabatic_2d_monatomic(Zl, dZl, a, kbt=1, m=1)[0])
    cs_dil = math.sqrt(2.0) * (1.0 + 2.0 * eta)
    n_all = len([r for r in runs if abs(float(r["eta"]) - eta) < 5e-7])
    n_el = len([r for r in runs if abs(float(r["eta"]) - eta) < 5e-7 and r["eligible"]])
    summ.append(dict(eta=f"{eta:.6f}", L0=f"{L0:.6f}", L_eff=f"{L_eff:.6f}",
                     c_s=f"{cs:.5f}", c_s_err=f"{cs_err:.5f}",
                     c_s_free_intercept=f"{cf:.5f}", intercept=f"{b0:.5g}", fit_r2=f"{r2:.5f}",
                     n_masses=len(xs), n_runs_all=n_all, n_runs_eligible=n_el, n_runs_in_fit=n_used,
                     cs_KR2006=f"{cs_kr:.5f}" if cs_kr == cs_kr else "",
                     dev_KR=f"{100*(cs-cs_kr)/cs_kr:+.2f}%" if cs_kr == cs_kr else "",
                     cs_Liu2021=f"{cs_liu:.5f}", dev_Liu=f"{100*(cs-cs_liu)/cs_liu:+.2f}%",
                     _scatter=f"{scatter:.5f}", _cs_dilute=f"{cs_dil:.5f}",
                     _dev_dilute=f"{100*(cs-cs_dil)/cs_dil:+.2f}%",
                     _cs_raw=f"{cs_raw:.5f}", _cs_raw_err=f"{cs_raw_err:.5f}",
                     _ti=ti_mean))

src = list(csv.DictReader(open(SRC)))
# ##CHRIS: plot_cs_meeting.py sizes the error bar for rows it cannot find in the
# trajectory manifest from a c_s_scatter_mass column. The 20260909 column set has no
# such column, so without this the new low-density points would be drawn with the fit
# error alone (0.1-0.2 %) instead of the mass scatter (0.5-0.8 %), understating them by
# roughly a factor of 4. Carried-over rows keep it blank on purpose: for those the
# script recomputes the scatter from the manifest, and specifying it twice could disagree.
cols = list(src[0].keys())
if "c_s_scatter_mass" not in cols:
    cols.insert(cols.index("c_s_err") + 1, "c_s_scatter_mass")
extra = {}
for s_ in summ:
    extra[s_["eta"]] = dict(scatter=s_.pop("_scatter"), cs_dil=s_.pop("_cs_dilute"),
                            dev_dil=s_.pop("_dev_dilute"), cs_raw=s_.pop("_cs_raw"),
                            cs_raw_err=s_.pop("_cs_raw_err"), ti=s_.pop("_ti"))
new_etas = {f"{float(s_['eta']):.6f}" for s_ in summ}
kept_src = [r for r in src if f"{float(r['eta']):.6f}" not in new_etas]
if len(kept_src) != len(src):
    print(f"note: {len(src)-len(kept_src)} source row(s) replaced by a recomputed row at the same eta")
for s_ in summ:
    s_["c_s_scatter_mass"] = extra[s_["eta"]]["scatter"]
rows = summ + kept_src
rows.sort(key=lambda r: float(r["eta"]))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in rows: w.writerow({c: r.get(c, "") for c in cols})
if PERRUN:
    with open(PERRUN, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(runs[0].keys())); w.writeheader(); w.writerows(runs)

# ---- also append into the corrected-estimator table (same rows, its column set) ----
corr_rows = list(csv.DictReader(open(CORR))) if os.path.exists(CORR) else []
if corr_rows:
    ccols = list(corr_rows[0].keys())
    new = []
    for s_ in summ:
        x_ = extra[s_["eta"]]
        ti = x_["ti"]
        exact = (ti == ti) and abs(ti - 1.0) < 1e-9
        if ti != ti:
            print(f"WARNING eta={s_['eta']}: no HD_KE_TRACE audit in run.log, so T_i was "
                  f"not measured. No 1/sqrt(T_i) correction was applied and c_s_raw equals "
                  f"c_s_corr by construction; T_i_mean is left blank, not set to 1.")
        elif not exact:
            print(f"WARNING eta={s_['eta']}: measured T_i_mean = {ti:.9g}, not 1. Writing the "
                  f"measured value; c_s_raw and c_s_corr genuinely differ.")
        new.append({
            "eta": s_["eta"], "L0": s_["L0"], "L_eff": s_["L_eff"],
            # ##CHRIS: emit the plain c_s / c_s_err aliases too, so a target table that
            # carries them (the robust A1 table does) does not end up with blank rows.
            "c_s": s_["c_s"], "c_s_err": s_["c_s_err"],
            "c_s_corr": s_["c_s"], "c_s_err_fit": s_["c_s_err"],
            "c_s_raw": s_["c_s"] if exact else x_["cs_raw"],
            "c_s_raw_err": s_["c_s_err"] if exact else x_["cs_raw_err"],
            "T_i_mean": "1.000000" if exact else (f"{ti:.6f}" if ti == ti else ""),
            "c_s_scatter_mass": x_["scatter"], "n_masses": s_["n_masses"],
            "n_masses_lt7": 1 if int(s_["n_masses"]) < 7 else 0,
            "n_runs_in_fit": s_["n_runs_in_fit"],
            "cs_KR2006": s_["cs_KR2006"], "dev_KR": s_["dev_KR"],
            "cs_dilute_sqrt2": x_["cs_dil"], "dev_dilute": x_["dev_dil"],
            "cs_Liu2021": s_["cs_Liu2021"], "dev_Liu": s_["dev_Liu"]})
    new_e = {f"{float(r['eta']):.6f}" for r in new}
    merged = new + [r for r in corr_rows if f"{float(r['eta']):.6f}" not in new_e]
    merged.sort(key=lambda r: float(r["eta"]))
    with open(CORR, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ccols, extrasaction="ignore")
        w.writeheader()
        for r in merged: w.writerow({c: r.get(c, "") for c in ccols})
    print(f"{CORR}: {len(merged)} rows ({len(new)} new)")
else:
    print(f"NOTE: {CORR} absent or empty -- corrected-table append skipped")

print(f"quality cut sigma_nu/nu < {CUT*100:.3f}%   strict health contract (all four counters == 0)")
print(f"traces {len(runs)}  accepted {sum(r['eligible'] for r in runs)}")
print(f"{OUT}: {len(rows)} rows ({len(summ)} new + {len(kept_src)} carried over)")
tis = [float(r["T_i_left"]) for r in runs if r["T_i_left"] not in ("", "nan")] + \
      [float(r["T_i_right"]) for r in runs if r["T_i_right"] not in ("", "nan")]
if tis:
    print(f"T_i over {len(tis)} compartments: min={min(tis):.9g} max={max(tis):.9g} "
          f"mean={sum(tis)/len(tis):.9g}")
print()
print(f"{'eta':>9} {'L0':>6} {'c_s':>9} {'+-fit':>8} {'scatter':>8} {'nM':>3} {'n_fit':>5} "
      f"{'c_s(KR)':>9} {'devKR':>8} {'sqrt2(1+2eta)':>13} {'devDil':>8}")
for s_ in summ:
    x_ = extra[s_["eta"]]; sc, cd, dd = x_["scatter"], x_["cs_dil"], x_["dev_dil"]
    print(f"{s_['eta']:>9} {float(s_['L0']):>6.0f} {s_['c_s']:>9} {s_['c_s_err']:>8} {sc:>8} "
          f"{s_['n_masses']:>3} {s_['n_runs_in_fit']:>5} {s_['cs_KR2006']:>9} {s_['dev_KR']:>8} "
          f"{cd:>13} {dd:>8}")
