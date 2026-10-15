#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.5-4.8.7; gate G5): THE STAGE 1 VERDICT SCRIPT, Part B
(blocks B, D, S). It reads the fetched groups (fetch_s1_koa.sh: <root>/<group>/merged/stage1/gen3/<cell>/m_<M>/ with red_nu.csv,
nud_runs.csv, acf_runs.npz, runs_health.csv, failed_runs.csv, s1_runs.csv, snap_runs.csv) and prints, by the registered rules:
  HEALTH (sec. 4.8.5)  per (cell, mass) the trajectories expected, present, excluded (clean != 1, or refused: failed_runs.csv), missing, and the
                       nu_d fit failures; the GROUP RULE (more than 2 of 25 lost, more than 1 of 10 in block S: the (cell, mass) is flagged and
                       left out of that cell's fits); the STOP RULE (more than 1 % of a block's trajectories excluded: that block's analysis
                       stops and is reported).
  E1 (sec. 4.8.6.2)    per (cell, mass): nu_d-bar over the used seeds (nud_runs.csv: the item-4a fit, nu_ref = the folder's own argmax mean,
                       decision 39 item 5 d), tau_r from the item-4a fit of the seed-averaged ACF (the ACFs cut to the shortest, the same
                       nu_ref), nu_0 = sqrt(nu_d-bar^2 + (1 / (2 pi tau_r))^2), Q = pi nu_d-bar tau_r, its SE the delete-one jackknife over the
                       trajectories (nu_d-bar and tau_r recomputed without each); c_k = nu_0 / x_k, x_k = K(alpha) / (2 pi L_eff), L_eff = L0 -
                       1.025; per cell c_s = the inverse-variance mean over the eight masses alpha >= 1, chi2 with (masses - 1) degrees of
                       freedom, sigma = (sum w)^(-1/2) max(1, sqrt(chi2 / dof)); no E1 with fewer than 5 masses; delta_0.5 = c_0.5 / c_s - 1;
                       the low-Q flag (a used mass with Q < 5).
  E0 (sec. 4.8.6.3)    the canonical estimator (paper1_populate_cs_err_20261002: argmax, TD = 200, X_EDGE = 2.5, the unweighted through-origin
                       slope over nine masses; T.slope and slope_with_errors), here with x_M = K(M/N) / (2 pi L_eff) (the exact box: no
                       truncation); printed next to E1, E1 - E0 the estimator difference. It decides nothing.
  ORIENTATION (4.8.3.1) per cell the crystal orientation the start chose, from the run record ([EDMD3-RUN]'s left/right fields) by M4c's
                       g3_lat_cryst (f72d764): close-packed lines along x or y in each compartment.
  Z, psi6 (4.8.6.4)    per trajectory s1_runs.csv's Z (the last quarter of the hold) and psi6_rec; per cell the mean and SE over every used
                       trajectory; Gamma = 2 / tau_r, Q, the mode power B C(0) per (cell, mass); the interior density from the snapshots.
  STIFFNESS (4.8.6.5)  k_S^dyn = N_s (2 pi nu_0)^2 / K(alpha)^2 over the five masses alpha >= 5, inverse-variance; the static side and rho_I
                       need Part A (they are printed as waiting).
  FITS (4.8.7.1)       per eta over the block B cells: Delta_N = c_s(N) / c_s^KR(eta) - 1 (E1; KR the rho_max = 0.90 fit of the project's
                       module) against u = N^(-1/2): F2 a + b u + d u^2 (primary), F1 a + b u, F1 without N = 100; weighted least squares
                       with the seed errors; F2 p < 0.01 flags "form not adequate"; fewer than four usable sizes: no F2, listed.
  VERDICT K (4.8.7.2)  on Phi (the 12 densities <= 0.60): z = a / sigma_a from F2, X^2 = sum z^2; KR-CONSISTENT if p(X^2) >= 0.01 and
                       max |z| < 3.5, else KR-INCONSISTENT with the densities; the weighted mean a-bar over Phi (F2 and F1), also without
                       flagged cells and densities; the approach cells' fits without a verdict.
  PRESSURE (4.8.7.3)   the same fits on Delta^Z_N = Z(N) / Z_KR(eta) - 1: a_Z, b_Z, the chi2 of F2; the expectation |a_Z| < 0.3 % on Phi.
  VERDICT I (4.8.7.4)  waits for Part A (decision 39 item 5 e): printed as waiting.
  VERDICT S (4.8.7.5)  per block S cell (N, eta): per start (i) Z and (ii) psi6 over the 20 trajectories of both masses, (iii) nu_d-bar at
                       alpha = 1 and (iv) at alpha = 10 over 10 each; z = (mean_0 - mean_turned) / sqrt(SE_0^2 + SE_turned^2); START-
                       INDEPENDENT if all four |z| < 3.5, else START-DEPENDENT (that (N, eta) is reported per start and left out of every
                       pooled statement); the main block's default start next to them (information).
  STATIONARITY (4.8.7.6) per cell: the hold, d = Z(4th quarter) - Z(3rd quarter) per trajectory, t = mean / SE, flag |t| >= 3; the record,
                       dnu = nu2 / nu1 - 1 (the halves of sec. 4.7.7), pooled over the masses, flag |z| >= 3.
  CONTINUITY (4.8.7.7) C-a E0 at (100, pi/8) against A1 v2's 3.80884 +- 0.011193 (the canonical table, c_s_err_scaled); C-b mean nu_arg at
                       (100, pi/8), M = 300 and 1500, against Test G on KOA's gen3 arm 0.0370087 +- 2.35e-05 and 0.0172607 +- 5.87e-06
                       (verdict_testG_koa_output.txt); C-c Z at (400, pi/8) against 2.81977 +- 0.00013 (notes sec. 4.7.28). Expected |z| < 3.
  PREDICTIONS (4.8.7.8) printed with their outcome where the data decide them (4 waits for Part A).
INTERIM (sec. 4.8.9): --interim prints the per-cell tables of the groups present, labelled INTERIM, without the fits and verdicts; the
fits, Verdict K and the pressure fits are printed once, when all four sizes are present.
usage (from hspist3/):
  python3 cluster/gen3_stage1_261012/s1_verdict_261012.py --design
  python3 cluster/gen3_stage1_261012/s1_verdict_261012.py --selftest [--tmp DIR]
  python3 cluster/gen3_stage1_261012/s1_verdict_261012.py --data <root> [--interim]
"""
import argparse, glob, math, os, re, sys, tempfile, warnings
import numpy as np, pandas as pd
from scipy.stats import chi2 as CHI2, norm, binom
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import gen_stage1_koa_261012 as GEN
import tests_20260913 as T
import resched_testT_followup_261007 as FU
import paper1_populate_cs_err_20261002 as E0
import plot_speed_of_sound_edmd as SOS
warnings.filterwarnings("ignore")
PHI = [e for e in GEN.ETA_B if e <= 0.60 + 1e-12]; APPROACH = [e for e in GEN.ETA_B if e > 0.60 + 1e-12]
NS = (100, 400, 900, 1600)
REF_CA = (3.80884, 0.011193)                                   # A1 v2 at eta_rec 0.392699: c_s, c_s_err_scaled (methods sec. 14.2's table)
REF_CB = {300: (0.0370087, 2.35e-05), 1500: (0.0172607, 5.87e-06)}   # Test G on KOA, gen3 arm, argmax mean +- SE (n = 400 each)
REF_CC = (2.81977, 0.00013)                                     # Z at (N 400, pi/8), the symmetry test (notes sec. 4.7.28)
Z_K, P_K, Z_S, T_STAT, Z_CONT = 3.5, 0.01, 3.5, 3.0, 3.0


# ---------------------------------------------------------------------------------------------------------------- the registered cells
def registered():
    rows, _ = GEN.plan()
    C = {}
    for r in rows:
        c = C.setdefault(r["cid"], dict(cid=r["cid"], g=r["g"], blk=r["blk"], N=r["N"], eta=r["eta"], start=r["start"], H=r["H"], L0=r["L0"],
                                        Ns=r["Ns"], Teq=r["Teq"], masses={}, nseed=r["nseed"]))
        c["masses"][r["M"]] = r["a"]
    return C


def kr_cs(eta): return T.kr_cs(eta)
def z_kr(eta): return float(SOS.Z_kolafa_rottner_2006(np.array([eta]))[0])


# ---------------------------------------------------------------------------------------------------------------- per (cell, mass)
def read_mass(d):
    f = lambda n: os.path.join(d, n)
    D = dict(red=pd.read_csv(f("red_nu.csv")) if os.path.exists(f("red_nu.csv")) else pd.DataFrame(columns=["run", "seed", "nu", "nu_pred", "dt", "n"]),
             nud=pd.read_csv(f("nud_runs.csv")) if os.path.exists(f("nud_runs.csv")) else None,
             health=pd.read_csv(f("runs_health.csv")) if os.path.exists(f("runs_health.csv")) else None,
             failed=pd.read_csv(f("failed_runs.csv")) if os.path.exists(f("failed_runs.csv")) else pd.DataFrame(columns=["folder", "run"]),
             s1=pd.read_csv(f("s1_runs.csv")) if os.path.exists(f("s1_runs.csv")) else None,
             snap=pd.read_csv(f("snap_runs.csv")) if os.path.exists(f("snap_runs.csv")) else None,
             acf=dict(np.load(f("acf_runs.npz"))) if os.path.exists(f("acf_runs.npz")) else {})
    D["cryst"] = None                                                # the default start's crystal orientation (sec. 4.8.3.1), from the run record
    if os.path.exists(f("run.log")):
        with open(f("run.log"), errors="replace") as fh:
            for l in fh:
                if l.startswith("[EDMD3-RUN] "):
                    sd = {m.group(1): cryst(m.group(2), int(m.group(4)), int(m.group(5)), float(m.group(6)), float(m.group(7))) for m in SIDE.finditer(l)}
                    D["cryst"] = (sd.get("left"), sd.get("right")); break
    return D


SIDE = re.compile(r"(left|right)=\[orient (parallel|perpendicular) rows (alternate|equal) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) ")


def cryst(orient, ns, nr, a_s, a_r):
    """f72d764 00ALLINONE.c's g3_lat_cryst (M4c): 0 = close-packed lines along x, 1 = along y (the turned crystal)"""
    along_rows = 1 if ns < 2 else (0 if nr < 2 else int(a_r < 2.0 * a_s))
    return 0 if along_rows == int(orient == "perpendicular") else 1


def seedavg_fit(acfs, dt, nu_ref):
    nl = min(len(a) for a in acfs); m = np.mean([np.asarray(a[:nl], dtype=float) for a in acfs], axis=0)
    try:
        p = FU.fit_acf(m, dt, nu_ref); return dict(A=p[0], tau_T=p[1], B=p[2], tau_r=p[3], omega=p[4], phi=p[5], C0=float(m[0]))
    except Exception:
        return None


def per_mass(cell, M, D, jackknife=True):
    a = cell["masses"][M]; exp = cell["nseed"]; red = D["red"]
    h = D["health"]; clean_runs = set(red["run"].astype(int))
    if h is not None and len(h):
        ok = {int(r.run) for r in h.itertuples() if int(r.gen3_clean) == 1 and int(r.gen3_records) == 1}
        clean_runs &= ok
    excluded = (len(red) - len(clean_runs)) + len(D["failed"])
    nud = D["nud"]
    fit_ok = {int(r.run) for r in nud.itertuples() if int(r.fit_ok) == 1} if nud is not None else set()
    used = sorted(clean_runs & fit_ok); fitfail = len(clean_runs - fit_ok)
    missing = max(0, exp - len(red) - len(D["failed"]))
    lost = exp - len(used)
    flagged = lost > (1 if cell["blk"] == "S" else 2)
    out = dict(cid=cell["cid"], M=M, alpha=a, expected=exp, present=len(red), excluded=excluded, missing=missing, fitfail=fitfail,
               used=len(used), lost=lost, flagged=flagged)
    rc = red[red["run"].isin(sorted(clean_runs))]
    out["nu_arg"] = float(rc["nu"].mean()) if len(rc) else math.nan
    out["nu_arg_se"] = float(rc["nu"].std(ddof=1) / math.sqrt(len(rc))) if len(rc) > 1 else math.nan
    out["nu_arg_sd"] = float(rc["nu"].std(ddof=1)) if len(rc) > 1 else math.nan; out["n_arg"] = len(rc)
    if len(used) < 3 or not D["acf"]:
        out.update(nud=math.nan, nud_se=math.nan, tau_r=math.nan, nu0=math.nan, nu0_se=math.nan, Q=math.nan, c=math.nan, sig_c=math.nan); return out
    nd = nud.set_index("run").loc[used, "nu_d"].to_numpy(float)
    dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    acfs = [D["acf"][f"run{r}"] for r in used]
    f = seedavg_fit(acfs, dt, nu_ref)
    out["nud"] = float(nd.mean()); out["nud_se"] = float(nd.std(ddof=1) / math.sqrt(len(nd)))
    if f is None:
        out.update(tau_r=math.nan, nu0=math.nan, nu0_se=math.nan, Q=math.nan, c=math.nan, sig_c=math.nan, meanfit_ok=0); return out
    nu0 = math.sqrt(out["nud"] ** 2 + (1 / (2 * math.pi * f["tau_r"])) ** 2)
    jk = []
    if jackknife:
        for i in range(len(used)):
            keep = [k for k in range(len(used)) if k != i]
            fi = seedavg_fit([acfs[k] for k in keep], dt, nu_ref)
            if fi is None: continue
            jk.append(math.sqrt(nd[keep].mean() ** 2 + (1 / (2 * math.pi * fi["tau_r"])) ** 2))
    n = len(jk); se = math.sqrt((n - 1) / n * sum((v - np.mean(jk)) ** 2 for v in jk)) if n >= 3 else math.nan
    x = T.k_root(a) / (2 * math.pi * (cell["L0"] - 1.025))
    out.update(tau_r=f["tau_r"], nu0=nu0, nu0_se=se, jk_fits=n, Q=math.pi * out["nud"] * f["tau_r"], Gamma=2 / f["tau_r"], B=f["B"],
               P1=f["B"] * f["C0"], x=x, c=nu0 / x, sig_c=se / x, meanfit_ok=1)
    return out


# ---------------------------------------------------------------------------------------------------------------- per cell
def e1(rows):
    u = [r for r in rows if r["alpha"] >= 1 and not r["flagged"] and math.isfinite(r.get("c", math.nan)) and math.isfinite(r.get("sig_c", math.nan)) and r["sig_c"] > 0]
    if len(u) < 5: return None
    c = np.array([r["c"] for r in u]); s = np.array([r["sig_c"] for r in u]); w = 1 / s ** 2
    cs = float((w * c).sum() / w.sum()); ch = float((w * (c - cs) ** 2).sum()); dof = len(u) - 1
    return dict(cs=cs, sig=float(w.sum() ** -0.5 * max(1.0, math.sqrt(ch / dof))), chi2=ch, dof=dof, n=len(u),
                lowQ=any(r["Q"] < 5 for r in u), masses=[r["M"] for r in u])


def e0(cell, rows):
    cs = [r for r in rows if r["n_arg"] > 0 and math.isfinite(r["nu_arg"])]
    if len(cs) < 3: return None
    cs = sorted(cs, key=lambda r: r["M"])
    x = np.array([T.k_root(r["alpha"]) / (2 * math.pi * (cell["L0"] - 1.025)) for r in cs]); y = np.array([r["nu_arg"] for r in cs])
    sy = np.array([r["nu_arg_sd"] / math.sqrt(r["n_arg"]) if r["n_arg"] > 1 else np.nan for r in cs])
    sT, scat, nm = T.slope(x, y); s2, err, errs, chi2r = E0.slope_with_errors(x, y, sy)
    return dict(cs=sT, err=err, err_scaled=errs, chi2r=chi2r, n=nm)


def mean_se(v):
    v = np.asarray([x for x in v if math.isfinite(x)], dtype=float)
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v)) if len(v) > 1 else (math.nan, math.nan, len(v))


def per_cell(cell, rows, D_by_M):
    out = dict(cid=cell["cid"], g=cell["g"], blk=cell["blk"], N=cell["N"], eta=cell["eta"], start=cell["start"])
    cr = {D["cryst"] for D in D_by_M.values() if D.get("cryst")}
    out["cryst"] = ("/".join("x" if c == 0 else "y" for c in next(iter(cr))) if len(cr) == 1 else ("differs between masses" if cr else "-"))
    out["E1"] = e1(rows); out["E0"] = e0(cell, rows)
    r05 = next((r for r in rows if r["alpha"] == 0.5), None)
    if out["E1"] and r05 and math.isfinite(r05.get("c", math.nan)) and not r05["flagged"]:
        cs, sg = out["E1"]["cs"], out["E1"]["sig"]; d = r05["c"] / cs - 1
        out["delta05"] = (d, math.hypot(r05["sig_c"] / cs, r05["c"] * sg / cs ** 2))
    Z, Ps, dh, dn = [], [], [], []
    used = {r["M"]: r for r in rows}
    for M, D in D_by_M.items():
        if D["s1"] is None or not len(D["s1"]) or used[M]["flagged"]: continue
        ok = set(D["red"]["run"].astype(int))
        if D["health"] is not None and len(D["health"]):
            ok &= {int(r.run) for r in D["health"].itertuples() if int(r.gen3_clean) == 1}
        s = D["s1"][D["s1"]["run"].isin(sorted(ok))]
        Z += list(s["Z"]); Ps += list(s["psi6_rec"]); dh += list(s["Z_q4"] - s["Z_q3"]); dn += list(s["dnu"])
    out["Z"] = mean_se(Z); out["psi6"] = mean_se(Ps)
    m, se, n = mean_se(dh); out["hold_t"] = (m / se if se and math.isfinite(se) and se > 0 else math.nan, n)
    m, se, n = mean_se(dn); out["rec_z"] = (m / se if se and math.isfinite(se) and se > 0 else math.nan, n, m)
    snaps = [D["snap"] for D in D_by_M.values() if D["snap"] is not None and len(D["snap"])]
    if snaps:
        sn = pd.concat(snaps); out["eta_int"] = (float(np.nanmean(sn[["eta_int_record_0", "eta_int_record_1"]].to_numpy())), len(sn))
    ks = [r for r in rows if r["alpha"] >= 5 and not r["flagged"] and math.isfinite(r.get("nu0", math.nan)) and math.isfinite(r.get("nu0_se", math.nan))]
    if len(ks) >= 3:
        k = np.array([cell["Ns"] * (2 * math.pi * r["nu0"]) ** 2 / T.k_root(r["alpha"]) ** 2 for r in ks])
        sk = np.array([cell["Ns"] * (2 * math.pi) ** 2 * 2 * r["nu0"] * r["nu0_se"] / T.k_root(r["alpha"]) ** 2 for r in ks]); w = 1 / sk ** 2
        out["kS_dyn"] = (float((w * k).sum() / w.sum()), float(w.sum() ** -0.5), len(ks))
    return out


# ---------------------------------------------------------------------------------------------------------------- fits and verdicts
def wls(u, y, s, deg):
    X = np.vstack([u ** k for k in range(deg + 1)]).T; W = 1 / s ** 2
    A = X.T @ (X * W[:, None]); beta = np.linalg.solve(A, X.T @ (W * y)); cov = np.linalg.inv(A)
    ch = float((W * (y - X @ beta) ** 2).sum()); dof = len(y) - (deg + 1)
    return dict(beta=beta, se=np.sqrt(np.diag(cov)), chi2=ch, dof=dof, p=float(CHI2.sf(ch, dof)) if dof > 0 else math.nan)


def fits_for(points):
    """points: [(N, Delta, sigma)] -> F2, F1, F1 without N = 100"""
    pts = sorted(p for p in points if all(math.isfinite(v) for v in p) and p[2] > 0)
    out = dict(n=len(pts), Ns=[p[0] for p in pts])
    if len(pts) >= 2:
        u = np.array([p[0] ** -0.5 for p in pts]); y = np.array([p[1] for p in pts]); s = np.array([p[2] for p in pts])
        if len(pts) >= 4: out["F2"] = wls(u, y, s, 2)
        if len(pts) >= 3: out["F1"] = wls(u, y, s, 1)
        q = [p for p in pts if p[0] != 100]
        if len(q) >= 3:
            out["F1x"] = wls(np.array([p[0] ** -0.5 for p in q]), np.array([p[1] for p in q]), np.array([p[2] for p in q]), 1)
    return out


def verdict_k(fits):
    z = {e: f["F2"]["beta"][0] / f["F2"]["se"][0] for e, f in fits.items() if "F2" in f}
    if not z: return None
    X2 = sum(v * v for v in z.values()); dof = len(z); p = float(CHI2.sf(X2, dof)); mz = max(abs(v) for v in z.values())
    ok = p >= P_K and mz < Z_K
    return dict(z=z, X2=X2, dof=dof, p=p, maxz=mz, consistent=ok, bad=[e for e, v in z.items() if abs(v) >= Z_K])


def wmean(vals):
    v = [(m, s) for m, s in vals if math.isfinite(m) and math.isfinite(s) and s > 0]
    if not v: return (math.nan, math.nan, 0)
    w = np.array([1 / s ** 2 for _, s in v]); m = np.array([x for x, _ in v])
    return (float((w * m).sum() / w.sum()), float(w.sum() ** -0.5), len(v))


def zdiff(a, sa, b, sb): return (a - b) / math.hypot(sa, sb) if math.isfinite(a) and math.isfinite(b) and math.hypot(sa, sb) > 0 else math.nan


# ---------------------------------------------------------------------------------------------------------------- the analysis
def analyse(root, interim=False, jackknife=True, quiet=False):
    C = registered(); P = print if not quiet else (lambda *a, **k: None)
    res = dict(mass={}, cell={}, groups_present=[], stopped=[])
    for g in GEN.GROUPS:
        base = os.path.join(root, g, "merged", "stage1", "gen3")
        if not os.path.isdir(base): continue
        res["groups_present"].append(g)
        cells = [c for c in C.values() if c["g"] == g]
        exp = sum(c["nseed"] * len(c["masses"]) for c in cells); excl = 0; mrows = {}; Ds = {}
        for c in cells:
            for M in c["masses"]:
                d = os.path.join(base, c["cid"], f"m_{M}")
                D = read_mass(d) if os.path.isdir(d) else dict(red=pd.DataFrame(columns=["run", "nu", "dt"]), nud=None, health=None,
                                                                failed=pd.DataFrame(columns=["folder", "run"]), s1=None, snap=None, acf={})
                Ds[(c["cid"], M)] = D
                r = per_mass(c, M, D, jackknife=jackknife); mrows[(c["cid"], M)] = r; excl += r["excluded"]
        if excl > 0.01 * exp:
            res["stopped"].append((g, excl, exp)); P(f"STOP RULE (sec. 4.8.5): group {g}: {excl} of {exp} trajectories excluded (> 1 %); its analysis stops"); continue
        res["mass"].update(mrows)
        for c in cells:
            rows = [mrows[(c["cid"], M)] for M in sorted(c["masses"])]
            res["cell"][c["cid"]] = per_cell(c, rows, {M: Ds[(c["cid"], M)] for M in c["masses"]})
    # fits per eta (block B); START-DEPENDENT (N, eta) left out
    S = {}
    for c in C.values():
        if c["blk"] != "S" or c["start"] != "0": continue
        a = res["cell"].get(c["cid"]); b = res["cell"].get(c["cid"].replace("_start0", "_turned"))
        if not a or not b: continue
        if not (math.isfinite(a["Z"][0]) and math.isfinite(b["Z"][0])):             # a start without data: no verdict, "no data"
            S[(c["N"], c["eta"])] = dict(z=dict(Z=math.nan, psi6=math.nan, nud_a1=math.nan, nud_a10=math.nan), independent=None); continue
        zs = dict(Z=zdiff(a["Z"][0], a["Z"][1], b["Z"][0], b["Z"][1]), psi6=zdiff(a["psi6"][0], a["psi6"][1], b["psi6"][0], b["psi6"][1]))
        for al, lab in ((1, "nud_a1"), (10, "nud_a10")):
            Ma = int(round(al * c["N"])); ra = res["mass"].get((a["cid"], Ma)); rb = res["mass"].get((b["cid"], Ma))
            zs[lab] = zdiff(ra["nud"], ra["nud_se"], rb["nud"], rb["nud_se"]) if ra and rb else math.nan
        S[(c["N"], c["eta"])] = dict(z=zs, independent=all(math.isfinite(v) and abs(v) < Z_S for v in zs.values()))
    res["S"] = S
    dep = {k for k, v in S.items() if v["independent"] is False}
    fitsC, fitsZ = {}, {}
    for e in GEN.ETA_B:
        pc, pz = [], []
        for N in NS:
            if (N, e) in dep: continue
            cell = res["cell"].get(GEN.cell_id("B", N, e, "default"))
            if not cell: continue
            if cell["E1"]:
                kr = kr_cs(e); pc.append((N, cell["E1"]["cs"] / kr - 1, cell["E1"]["sig"] / kr))
            zk = z_kr(e)
            if math.isfinite(cell["Z"][0]): pz.append((N, cell["Z"][0] / zk - 1, cell["Z"][1] / zk))
        fitsC[e] = fits_for(pc); fitsZ[e] = fits_for(pz)
    res["fitsC"], res["fitsZ"] = fitsC, fitsZ
    all4 = all(g in res["groups_present"] for g in ("B100", "B400", "B900", "B1600"))
    res["K"] = verdict_k({e: fitsC[e] for e in PHI}) if all4 and not interim else None
    return res


def report(res, interim=False):
    lab = "INTERIM -- " if interim else ""
    print(f"# {lab}STAGE 1, PART B: the tables and verdicts by the registered rules (s1_verdict_261012.py)\n")
    print(f"groups present: {', '.join(res['groups_present']) or 'none'}; stopped by the stop rule: {res['stopped'] or 'none'}\n")
    print("## Health per group (sec. 4.8.5)\n\n| group | (cell, mass) | expected | present | excluded | missing | fit failures | flagged (group rule) |\n|---|---|---|---|---|---|---|---|")
    for g in res["groups_present"]:
        rr = [r for (cid, M), r in res["mass"].items() if res["cell"].get(cid, {}).get("g") == g]
        print(f"| {g} | {len(rr)} | {sum(r['expected'] for r in rr)} | {sum(r['present'] for r in rr)} | {sum(r['excluded'] for r in rr)} | "
              f"{sum(r['missing'] for r in rr)} | {sum(r['fitfail'] for r in rr)} | {sum(r['flagged'] for r in rr)} |")
    print("\n## Per cell: E1, E0, Z, psi6, the stationarity tests, delta_0.5, k_S^dyn\n")
    print("| cell | lines (left/right) | E1 c_s +- sigma (masses, chi2/dof) | low-Q | E0 c_s +- err_scaled | E1 - E0 | KR | Z +- SE (n) | psi6 | hold t | record z | "
          "delta_0.5 [%] | k_S^dyn | eta_int |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for cid, c in sorted(res["cell"].items(), key=lambda kv: (GEN.GROUPS.index(kv[1]["g"]), kv[1]["N"], kv[1]["eta"], kv[1]["start"])):
        e1v, e0v = c["E1"], c["E0"]; kr = kr_cs(c["eta"])
        s1 = f"{e1v['cs']:.5f} +- {e1v['sig']:.5f} ({e1v['n']}, {e1v['chi2']:.1f}/{e1v['dof']})" if e1v else "none (< 5 masses)"
        s0 = f"{e0v['cs']:.5f} +- {e0v['err_scaled']:.5f}" if e0v else "-"
        dd = f"{e1v['cs'] - e0v['cs']:+.5f}" if e1v and e0v else "-"
        d5 = f"{100 * c['delta05'][0]:+.3f} +- {100 * c['delta05'][1]:.3f}" if c.get("delta05") else "-"
        ks = f"{c['kS_dyn'][0]:.4g} +- {c['kS_dyn'][1]:.2g}" if c.get("kS_dyn") else "-"
        ei = f"{c['eta_int'][0]:.4f} ({c['eta_int'][1]})" if c.get("eta_int") else "-"
        print(f"| {cid} | {c['cryst']} | {s1} | {'FLAG' if e1v and e1v['lowQ'] else '-'} | {s0} | {dd} | {kr:.5f} | {c['Z'][0]:.5f} +- {c['Z'][1]:.5f} ({c['Z'][2]}) | "
              f"{c['psi6'][0]:.4f} | {c['hold_t'][0]:+.2f}{' FLAG' if abs(c['hold_t'][0]) >= T_STAT else ''} | "
              f"{c['rec_z'][0]:+.2f}{' FLAG' if abs(c['rec_z'][0]) >= T_STAT else ''} | {d5} | {ks} | {ei} |")
    if interim:
        print("\nINTERIM: the fits, Verdict K, the pressure fits and Verdict I are printed once, when all four sizes are present (sec. 4.8.9)"); return
    print("\n## Verdict S (sec. 4.8.7.5), per block S cell\n\n| N | eta | z Z | z psi6 | z nu_d (alpha 1) | z nu_d (alpha 10) | verdict |\n|---|---|---|---|---|---|---|")
    for (N, e), v in sorted(res["S"].items()):
        z = v["z"]; print(f"| {N} | {e:.2f} | {z['Z']:+.2f} | {z['psi6']:+.2f} | {z['nud_a1']:+.2f} | {z['nud_a10']:+.2f} | "
                         f"{'no data' if v['independent'] is None else ('START-INDEPENDENT' if v['independent'] else '**START-DEPENDENT**')} |")
    print("\n## The finite-size fits per eta (sec. 4.8.7.1), E1: Delta_N = c_s(N) / c_s^KR - 1 [%]\n")
    print("| eta | sizes | F2 a | F2 b | F2 d | F2 p | form | F1 a | F1 p | F1 without N 100: a |\n|---|---|---|---|---|---|---|---|---|---|")
    for e, f in res["fitsC"].items():
        F2 = f.get("F2"); F1 = f.get("F1"); F1x = f.get("F1x")
        print(f"| {e:.4g} | {f['n']} | " + (f"{100 * F2['beta'][0]:+.3f} +- {100 * F2['se'][0]:.3f} | {100 * F2['beta'][1]:+.2f} | {100 * F2['beta'][2]:+.2f} | {F2['p']:.3f} | "
              f"{'not adequate' if F2['p'] < 0.01 else 'ok'}" if F2 else "no F2 (fewer than four usable sizes) | | | |") + " | "
              + (f"{100 * F1['beta'][0]:+.3f} +- {100 * F1['se'][0]:.3f} | {F1['p']:.3f}" if F1 else "- | -") + " | "
              + (f"{100 * F1x['beta'][0]:+.3f} +- {100 * F1x['se'][0]:.3f}" if F1x else "-") + " |")
    K = res["K"]
    if K:
        print(f"\nVERDICT K (sec. 4.8.7.2) on Phi ({K['dof']} densities with F2): X^2 = {K['X2']:.2f} ({K['dof']} dof), p = {K['p']:.4f}; max |z| = {K['maxz']:.2f}: "
              + ("**KR-CONSISTENT**" if K["consistent"] else f"**KR-INCONSISTENT** ({', '.join(f'{e:.4g}' for e in K['bad']) or 'by X^2'})"))
        a2 = wmean([(res["fitsC"][e]["F2"]["beta"][0], res["fitsC"][e]["F2"]["se"][0]) for e in PHI if "F2" in res["fitsC"][e]])
        a1 = wmean([(res["fitsC"][e]["F1"]["beta"][0], res["fitsC"][e]["F1"]["se"][0]) for e in PHI if "F1" in res["fitsC"][e]])
        print(f"weighted mean a over Phi: F2 {100 * a2[0]:+.3f} +- {100 * a2[1]:.3f} % ({a2[2]}); F1 {100 * a1[0]:+.3f} +- {100 * a1[1]:.3f} % ({a1[2]})")
    else:
        print("\nVERDICT K: not printed (all four sizes are needed)")
    print("\n## Pressure (sec. 4.8.7.3), no verdict: Delta^Z_N = Z(N) / Z_KR - 1 [%]\n\n| eta | a_Z | b_Z | F2 chi2 | |a_Z| < 0.3 % |\n|---|---|---|---|---|")
    for e, f in res["fitsZ"].items():
        F2 = f.get("F2")
        if F2: print(f"| {e:.4g} | {100 * F2['beta'][0]:+.4f} +- {100 * F2['se'][0]:.4f} | {100 * F2['beta'][1]:+.3f} | {F2['chi2']:.2f} | "
                     f"{'yes' if abs(F2['beta'][0]) < 0.003 else ('**no**' if e in PHI else 'approach cell')} |")
    print("\nVERDICT I (sec. 4.8.7.4): WAITS FOR PART A (decision 39 item 5 e); k_S^dyn is printed above, the static side needs Part A")
    print("\n## Continuity (sec. 4.8.7.7; expected |z| < 3, a miss is reported)\n")
    c = res["cell"].get(GEN.cell_id("B", 100, math.pi / 8, "default"))
    if c and c["E0"]: print(f"- C-a: E0 at (100, pi/8) {c['E0']['cs']:.5f} +- {c['E0']['err_scaled']:.5f} against A1 v2 {REF_CA[0]} +- {REF_CA[1]}: "
                            f"z = {zdiff(c['E0']['cs'], c['E0']['err_scaled'], *REF_CA):+.2f}")
    for M, (m, s) in REF_CB.items():
        r = res["mass"].get((GEN.cell_id("B", 100, math.pi / 8, "default"), M))
        if r: print(f"- C-b: nu_arg at (100, pi/8), M = {M}: {r['nu_arg']:.7f} +- {r['nu_arg_se']:.2e} against Test G on KOA {m} +- {s}: "
                    f"z = {zdiff(r['nu_arg'], r['nu_arg_se'], m, s):+.2f}")
    c = res["cell"].get(GEN.cell_id("B", 400, math.pi / 8, "default"))
    if c: print(f"- C-c: Z at (400, pi/8) {c['Z'][0]:.5f} +- {c['Z'][1]:.5f} against {REF_CC[0]} +- {REF_CC[1]}: z = {zdiff(c['Z'][0], c['Z'][1], *REF_CC):+.2f}")


# ---------------------------------------------------------------------------------------------------------------- --design
def design():
    print("# STAGE 1, PART B -- the verdict script's design (gate G5; sec. 4.8.7: every threshold with its false-alarm rate, before any data)\n")
    pz = 2 * norm.sf(Z_K)
    fk = 1 - (1 - P_K) * (1 - pz) ** 12
    print(f"- Verdict K: p(X^2, 12 dof) >= {P_K} and max over the 12 densities of |z| < {Z_K}: false alarm <= 1 - (1 - {P_K}) (1 - {pz:.3g})^12 = {fk:.4f} "
          f"(the registration: about 1.6 % = 1 % + 0.56 %; P(max |z| >= 3.5) = {1 - (1 - pz) ** 12:.4f})")
    fs = 1 - (1 - pz) ** 4
    print(f"- Verdict S: all four |z| < {Z_S}: false alarm {fs:.5f} per cell, {1 - (1 - fs) ** 20:.4f} for at least one of the 20 cells (the registration: 0.19 %, 3.7 %)")
    ps = 2 * norm.sf(T_STAT)
    print(f"- Stationarity: |t| >= {T_STAT} (hold) and |z| >= {T_STAT} (record): {ps:.4f} per test; about 220 tests: {220 * ps:.2f} false flags expected (registration: about 0.6)")
    print(f"- F2 'form not adequate': p < 0.01 (1 dof): 1 % per density")
    print(f"- Verdict I: p(chi2, 4 dof) >= 0.01 for rho_I = 0; else the c_I / N_s fit (3 dof) p >= 0.01: 1 % each (waits for Part A)")
    print(f"- Continuity C-a, C-b (two masses), C-c: |z| < {Z_CONT} expected: {2 * norm.sf(Z_CONT):.4f} per check, {1 - (1 - 2 * norm.sf(Z_CONT)) ** 4:.4f} for at least one of four")
    p2 = 2 * norm.sf(2.0); pf = 1 - binom.cdf(2, 12, p2)
    print(f"- Prediction 2: |z| < 2 at 10 or more of 12 densities: fails by chance with P(3 or more of 12 with |z| >= 2) = {pf:.4f} (registration: 1.5 %)")
    print("- The group rule: more than 2 of 25 (more than 1 of 10 in block S) lost in a (cell, mass): flagged, left out of the cell's fits; "
          "the stop rule: more than 1 % of a block's trajectories excluded: the block's analysis stops")
    print("- E1 needs 5 or more of the eight masses alpha >= 1; the low-Q flag at Q < 5; E0 needs 3 or more masses (the canonical function)")
    print(f"- The continuity references: C-a A1 v2 (N 100, pi/8) {REF_CA[0]} +- {REF_CA[1]} (methods sec. 14.2's table, c_s_err_scaled); C-b Test G on KOA, "
          f"gen3 arm, M 300 {REF_CB[300][0]} +- {REF_CB[300][1]}, M 1500 {REF_CB[1500][0]} +- {REF_CB[1500][1]}; C-c {REF_CC[0]} +- {REF_CC[1]} (notes sec. 4.7.28)")
    print(f"- KR: the project's module (plot_speed_of_sound_edmd: the rho_max = 0.90 fit); c_s^KR at the 16 densities: "
          + ", ".join(f"{e:.4g}: {kr_cs(e):.5f}" for e in GEN.ETA_B))
    print(f"- The fluid set Phi: {', '.join(f'{e:.4g}' for e in PHI)}; the approach cells: {', '.join(f'{e:.4g}' for e in APPROACH)}")


# ---------------------------------------------------------------------------------------------------------------- --selftest
def synth_mass(d, cell, M, rng, nu0, tau_r, Q_low=False, Zt=None, inject=None, nfail=0):
    """one synthetic (cell, mass) folder in the fetched layout: ACFs from the item-4a model with noise, the per-trajectory tables"""
    os.makedirs(d, exist_ok=True); n = cell["nseed"]; dt = 1 / (nu0 * 32); nl = int(20 / (nu0 * dt))
    t = np.arange(nl + 1) * dt; om = 2 * math.pi * math.sqrt(max(nu0 ** 2 - (1 / (2 * math.pi * tau_r)) ** 2, (0.5 * nu0) ** 2))
    acf, red, nud, health, s1 = {}, [], [], [], []
    for r in range(n - nfail):
        g = 0.3 * np.exp(-t / (2 / nu0)) + 0.7 * np.exp(-t / tau_r) * np.cos(om * t - 0.05) + rng.normal(0, 0.01, len(t))
        acf[f"run{r}"] = (g * 1e-3).astype(np.float32)
        nu_d = om / (2 * math.pi) * (1 + rng.normal(0, 0.003)); nu_a = nu_d * (1 - 1 / (8 * (math.pi * nu0 * tau_r) ** 2)) * (1 + rng.normal(0, 0.003))
        red.append(dict(run=r, seed=1000 + r, nu=nu_a, nu_pred=nu0, dt=dt, n=6400)); nud.append(dict(run=r, seed=1000 + r, nu_arg=nu_a, nu_d=nu_d, fit_ok=1, nu_ref=nu0, dt=dt, n=6400))
        health.append(dict(run=r, seed=1000 + r, gen3_clean=1, gen3_records=1))
        z = Zt * (1 + rng.normal(0, 2e-4)); zq3 = z * (1 + rng.normal(0, 2e-4)); dz = (inject or {}).get("hold", 0.0)
        s1.append(dict(run=r, seed=1000 + r, T_eq=cell["Teq"], Z=z + dz, Z_q3=zq3, Z_q4=z + dz, psi6_rec=0.1 + rng.normal(0, 0.01) + (inject or {}).get("psi6", 0.0),
                       dnu=rng.normal(0, 0.004)))
    np.savez_compressed(os.path.join(d, "acf_runs.npz"), **acf)
    pd.DataFrame(red).to_csv(os.path.join(d, "red_nu.csv"), index=False); pd.DataFrame(nud).to_csv(os.path.join(d, "nud_runs.csv"), index=False)
    pd.DataFrame(health).to_csv(os.path.join(d, "runs_health.csv"), index=False); pd.DataFrame(s1).to_csv(os.path.join(d, "s1_runs.csv"), index=False)
    pd.DataFrame([dict(folder=f".failed_run{n - 1 - k}_x", run=n - 1 - k) for k in range(nfail)], columns=["folder", "run"]).to_csv(os.path.join(d, "failed_runs.csv"), index=False)


def selftest(tmp):
    print("# STAGE 1 verdict script: --selftest\n")
    ok_all = []
    # (1) the decision rules by simulation
    rng = np.random.default_rng(20261012); R = 20000; fa = 0; pw = 0
    us = np.array([n ** -0.5 for n in NS]); s = np.array([0.0006, 0.0008, 0.0010, 0.0012])
    for k in range(R):
        f0 = {e: {"F2": wls(us, 0.01 * us + 0.02 * us ** 2 + rng.normal(0, s), s, 2)} for e in range(12)}
        f1 = {e: {"F2": wls(us, 0.004 + 0.01 * us + rng.normal(0, s), s, 2)} for e in range(12)}
        fa += not verdict_k(f0)["consistent"]; pw += not verdict_k(f1)["consistent"]
    se = math.sqrt(0.016 * 0.984 / R)
    ok1 = abs(fa / R - 0.016) < 4 * se
    print(f"1. Verdict K by simulation ({R} replicas, 12 densities x 4 sizes, a = 0): KR-INCONSISTENT in {fa / R:.4f} (design about 0.016, +- {4 * se:.4f}): "
          f"{'PASS' if ok1 else 'FAIL'}; power at a = +0.4 % everywhere: {pw / R:.3f}"); ok_all.append(ok1)
    # (2) end to end on a synthetic fetched tree
    C = registered(); root = os.path.join(tmp, "synthetic"); rng = np.random.default_rng(7)
    pick = [GEN.cell_id("B", N, e, "default") for N in NS for e in (0.10, math.pi / 8, 0.60, 0.66)] + [GEN.cell_id("S", N, 0.63, st) for N in (100, 400) for st in ("0", "turned")]
    truth = {}
    for cid in pick:
        c = C[cid]; kr = kr_cs(c["eta"]); u = c["N"] ** -0.5; cs_true = kr * (1 + 0.01 * u)
        truth[cid] = cs_true; Zt = z_kr(c["eta"]) * (1 + 0.02 * u)
        if cid == GEN.cell_id("S", 400, 0.63, "turned"): Zt *= 1.003                        # injected: a start difference in Z
        for M, al in c["masses"].items():
            x = T.k_root(al) / (2 * math.pi * (c["L0"] - 1.025)); nu0 = cs_true * x; Q = 20 + 3 * al
            if cid == GEN.cell_id("B", 400, 0.10, "default") and al == 0.5: Q = 3.0         # injected: low Q at the light mass (not used in E1)
            inj = {"hold": 0.002 * z_kr(c["eta"])} if cid == GEN.cell_id("B", 900, 0.60, "default") else None   # injected: a hold drift
            nfail = 3 if (cid == GEN.cell_id("B", 100, 0.60, "default") and al == 2) else 0                       # injected: the group rule
            synth_mass(os.path.join(root, c["g"], "merged", "stage1", "gen3", cid, f"m_{M}"), c, M, rng, nu0, Q / (math.pi * nu0), Zt=Zt, inject=inj, nfail=nfail)
    res = analyse(root, jackknife=True, quiet=True)
    chk = []
    pb = [cid for cid in pick if cid.startswith("B_")]
    zs = [(res["cell"][cid]["E1"]["cs"] - truth[cid]) / res["cell"][cid]["E1"]["sig"] for cid in pb if res["cell"].get(cid, {}).get("E1")]
    chk.append(("E1 recovers the constructed c_s (|z| < 3 in every synthetic block B cell; block S has 2 masses, no E1)", len(zs) == len(pb) and all(abs(z) < 3 for z in zs),
                f"{len(zs)} cells, max |z| {max(abs(z) for z in zs):.2f}"))
    fl = [(cid, M) for (cid, M), r in res["mass"].items() if r["flagged"] and (r["present"] or r["excluded"])]
    chk.append(("the group rule flags exactly the injected (cell, mass)", fl == [(GEN.cell_id("B", 100, 0.60, "default"), 200)], str(fl)))
    S = res["S"]; dep = [k for k, v in S.items() if v["independent"] is False]
    chk.append(("Verdict S: START-DEPENDENT exactly at the injected (400, 0.63)", dep == [(400, 0.63)], str(dep)))
    hf = [cid for cid, c in res["cell"].items() if abs(c["hold_t"][0]) >= T_STAT]
    chk.append(("the hold stationarity flags exactly the injected cell", hf == [GEN.cell_id("B", 900, 0.60, "default")], str(hf)))
    lq = [cid for cid, c in res["cell"].items() if c["E1"] and c["E1"]["lowQ"]]
    chk.append(("no low-Q flag from alpha = 0.5 (it is not used in E1)", lq == [], str(lq)))
    d5 = [c["delta05"][0] / c["delta05"][1] for c in res["cell"].values() if c.get("delta05")]
    chk.append(("delta_0.5 consistent with 0 (|z| < 3) where it exists", all(abs(z) < 3 for z in d5), f"{len(d5)} cells"))
    f2 = res["fitsC"][0.10].get("F2"); chk.append(("F2 at eta 0.10 recovers a = 0 (|z| < 3)", f2 is not None and abs(f2["beta"][0] / f2["se"][0]) < 3, f"{f2['beta'][0]:+.2e} +- {f2['se'][0]:.1e}" if f2 else "-"))
    print("\n2. End to end on a synthetic fetched tree (16 block B cells and 2 block S cells; injected: a start difference, a hold drift, 3 refused runs, a low Q)\n")
    print("| check | result | detail |\n|---|---|---|")
    for lab, ok, det in chk: print(f"| {lab} | {'PASS' if ok else 'FAIL'} | {det} |"); ok_all.append(ok)
    print(f"\nSELFTEST: {'PASS' if all(ok_all) else 'FAIL'} ({sum(ok_all)} of {len(ok_all)})")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--design", action="store_true"); ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--data"); ap.add_argument("--interim", action="store_true"); ap.add_argument("--tmp"); a = ap.parse_args()
    if a.design: design()
    if a.selftest: selftest(a.tmp or tempfile.mkdtemp(prefix="s1_verdict_selftest_"))
    if a.data: report(analyse(a.data, interim=a.interim), interim=a.interim)


if __name__ == "__main__":
    main()
