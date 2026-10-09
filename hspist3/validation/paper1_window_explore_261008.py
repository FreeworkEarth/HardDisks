#!/usr/bin/env python3
"""##CHRIS 2026-10-08 (261012 sec. 4.6; the plan author's CC task of 2026-10-08, Cowork clock): EXPLORATORY, POST HOC analysis of
the existing data in and near the melting window. No new runs, no engine work, nothing registered: it generates hypotheses for the
stage-1 pre-registration and DECIDES NOTHING. Every table below is printed by this script.
  item 1  inventory: every speed-of-sound cell under 00_eta_sweep_ROMAN with eta_true in [0.675, 0.735] (any build): campaign, N,
          H, L0, masses, seeds, record length, psi6 output, health lines, run date and build. Campaigns are never pooled; every
          figure shows ONE campaign (one binary).
  item 2  per cell and mass: c_app,M = nu_M / x_M with SE (x_M = K(alpha)/(2 pi L_eff,true), cot K = alpha K, alpha = M/(2 N_s);
          L_eff,true = L0 - 2r - t/2 - delta/2 as the registered method B), with
            nu_M  the argmax estimator of reduce_B.py / paper1_populate_cs_err_20261002.cell (TD = 200 where the record has 200
                  periods -- the canonical A1 v2 estimator, gated against the canonical table below -- and TD = the record's
                  planned periods, rounded down, for the shorter campaigns: an exploratory variant, labelled);
            nu_d  the refined frequency of 261012 sec. 4.4.13 item 4a (damped cosine with free phase plus the slow mode, fitted to
                  each trajectory's mean-removed position ACF, start = the mass's mean argmax nu);
            Gamma_M = 2/tau_r from the registered damping model (methods sec. 13; paper1_confinement_results_261004._fit_acf) fitted
                  to the seed-averaged ACF, 10-group jackknife SE; Q = pi nu tau_r.
          The deviation from the KR fluid value (ρmax = 0.90 fit, the module; validated to eta 0.7069, compared with data only to
          0.69: every KR value above 0.69 is printed but flagged, above 0.7069 it is extrapolation = INFERENCE).
          The per-mass dip D_M = c_app,M(cell)/c_app,M(entry) - 1, entry = the campaign's last cell below the window (it removes
          each mass's own offset, e.g. the plain-fluid light-mass deficit of methods sec. 16). Spearman of c_app,M and of D_M with
          alpha per cell; the slope of D_M against ln alpha (weighted) with z. Window = eta_true in [0.6995, 0.7175]
          (plot_speed_of_sound_edmd.ETA_COEX_LO..HI = 0.700..0.716, widened by the grid rounding of eta_true).
          The single-c_s chi2 (8 dof) at the unweighted (registered) and the weighted slope, next to the Test T plain-fluid
          baseline at the same 25 seeds per mass (four blocks of Test T legacy; the 279282b campaign anchor). Per campaign, the
          depth of the dip by the design's definition (261012 sec. 4.2) next to the design's M1 expectation.
  item 4b psi6 against mass and against the per-trajectory frequency (the record is a fixed number of a trajectory's OWN periods,
          so heavy dividers are measured over longer times while the structure is still changing: an aging confound).
  item 3  equilibrium prediction: c_s^2 = (kT/m)(Z + eta Z' + Z^2) for the fluid from KR (ρmax 0.90, the module; and the
          ρmax 0.88 fit, typed from the paper in paper1_kr_sanity_261002.KR) and from Henderson; the ideal coexistence plateau
          (Z + eta Z' = 0): c_0 = Z sqrt(kT/m), Z = P* pi/(4 eta), P* = beta P (2 sigma)^2 = 9.17 and 9.19 (sigma = disk radius;
          Engel et al., Eq. (1) and Table I / Figs. 3-4, PDF in ZZZ_PAPER/EOS/). The step at 0.700 for each fluid variant and
          the plan author's construction (Henderson's d ln Z/d eta at the plateau Z); the variation of c_0 over the window.
  item 4  divider periods per mass in the window (sigma-time); psi6: only per-run summaries exist (64 samples per run, mean/SD/
          min/max; 00ALLINONE.c:15307, :16108), no time series; their hold/end/run values per cell; the cost of a psi6(t) series
          and of a positions trace per trajectory.
Figures (PNG; PDF for A1 v2) in 0000_PLAN_OVERALL/paper1_speedofsound/experiments/exploratory_261008_window/.
usage (from hspist3/): python3 validation/paper1_window_explore_261008.py [--workers 10]
"""
import contextlib, csv, glob, io, math, os, re, sys, warnings
from collections import defaultdict
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import tests_20260913 as T
import plot_speed_of_sound_edmd as SOS
from paper1_populate_cs_err_20261002 import box_delta, cell as canon_cell, slope_with_errors
import paper1_confinement_results_261004 as R
import resched_testT_followup_261007 as F

ROOT = T.ROOT
OUT = os.path.join(T.PAPER1, "exploratory_261008_window")
CANON = os.path.join(T.PAPER1, "final", "260919_A1v2_final_cs_vs_eta.csv")
ETA_LO, ETA_HI = 0.675, 0.735
WIN_LO, WIN_HI = 0.6995, 0.7175
X_EDGE, KR_VALID, KR_CMP = 2.5, math.pi / 4 * 0.90, SOS.KR2006_PLOT_ETA_MAX
PSTAR = (9.17, 9.19)
SKIP = re.compile(r"/(merged|analysis[^/]*|an|replot[^/]*|_orchestration[^/]*|\.run\d*|\.failed[^/]*|\.stale[^/]*|[^/]*ABANDONED[^/]*)(/|$)")
PSI6_LOG = re.compile(r"psi6: hold=([\d.]+)\s+end=([\d.]+)\s+run_mean=([\d.]+)\+-([\d.]+)")
NEW_CORE = "2026-08-20"            # the r10 "newcore" date; older cells are listed, not analysed
INK, MUTED, GRID = "#1f1f1e", "#8a8984", "#e4e3df"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"     # dataviz reference palette, categorical slots 1, 2, 3, 7
RAMP4 = ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")                        # ordinal, validated (4 steps, light surface)


# ------------------------------------------------------------------------------------------------------------- data access
def header(p):
    with open(p, newline="") as fh:
        rd = csv.DictReader(fh); return next(rd)


def last_time(p):
    with open(p, "rb") as fh:
        fh.seek(0, 2); size = fh.tell(); fh.seek(max(0, size - 8192))
        return float(fh.read().decode(errors="replace").strip().splitlines()[-1].split(",")[0])


def meta(cell, p):
    """Geometry and record of a cell from its first trace; older traces lack some columns (fallbacks: the trace tail, 00_COMMAND.md)."""
    h = header(p); L0, eta = h.get("L0"), h.get("eta")
    if L0 is None or eta is None: return None
    nup = float(h["Predicted_Frequency"]) if h.get("Predicted_Frequency") else float("nan")
    per = float(h["Planned_Duration"]) * nup if h.get("Planned_Duration") else last_time(p) * nup
    Ns = h.get("Left_Count")
    if Ns is None:
        for c in (os.path.join(cell, "00_COMMAND.md"), os.path.join(os.path.dirname(cell), "00_COMMAND.md")):
            m = re.search(r"--particles-boxes=(\d+),", open(c, errors="ignore").read()) if os.path.exists(c) else None
            if m: Ns = m.group(1); break
    return dict(h=h, L0=float(L0), eta=float(eta), per=per, Ns=int(float(Ns)) if Ns is not None else 0,
                complete="Planned_Duration" in h and "Left_Count" in h)


def campaign_of(cell):
    parts = [q for q in os.path.relpath(cell, ROOT).split(os.sep) if not re.fullmatch(r"eta_0p\d+", q)]
    return "/".join(parts) if parts else "."


def discover():
    cells = defaultdict(lambda: defaultdict(list))
    for p in glob.glob(os.path.join(ROOT, "**", "wall_x_positions_L0_*_wallmassfactor_*_run*.csv"), recursive=True):
        if SKIP.search("/" + os.path.relpath(p, ROOT)): continue
        m = re.search(r"wallmassfactor_(\d+)_run(\d+)\.csv$", p); M, r = int(m.group(1)), int(m.group(2))
        d = os.path.dirname(p); cell = d[: -len(f"/m_{M}")] if d.endswith(f"/m_{M}") else d
        cells[cell][M].append((r, p))
    return cells


def date_of(cell):
    for lg in [os.path.join(cell, "run.log")] + sorted(glob.glob(os.path.join(cell, "m_*", "run.log"))):
        if os.path.exists(lg):
            m = re.search(r"##RUN (\d{4}-\d{2}-\d{2})", open(lg, errors="ignore").read(4000))
            if m: return m.group(1)
    for c in (os.path.join(cell, "00_COMMAND.md"), os.path.join(os.path.dirname(cell), "00_COMMAND.md")):
        if os.path.exists(c):
            m = re.search(r"Timestamp: (\d{4}-\d{2}-\d{2})", open(c, errors="ignore").read())
            if m: return m.group(1)
    m = re.search(r"(20\d{2})(\d{2})(\d{2})", cell)
    return f"{m.group(1)}-{m.group(2)}-{m.group(3)} (dir name)" if m else "unknown"


def logs_of(cell):
    return [os.path.join(cell, "run.log")] + sorted(glob.glob(os.path.join(cell, "m_*", "run.log")))


def build_of(cell):
    for lg in logs_of(cell):
        if os.path.exists(lg):
            s = open(lg, errors="ignore").read(200000)
            m = re.search(r"git ([0-9a-f]{7,}(?:-dirty)?)\s+target (\w+)", s) or re.search(r"build_git[\s=:,]+([0-9a-f]{7,}(?:-dirty)?)", s)
            if m: return " ".join(m.groups())
    return "unrecorded (Mac binary; no version line before the 2026-09-16 provenance change)"


def bad_runs(cell):
    bad, nlines, tele = set(), 0, False
    for lg in logs_of(cell):
        if os.path.exists(lg):
            s = open(lg, errors="ignore").read()
            tele |= "EDMD-HEALTH" in s or "health" in s.lower()
            for m in T.HEALTH_RE.finditer(s):
                nlines += 1
                if any(int(m.group(i)) for i in (4, 5, 6, 7)): bad.add((int(m.group(1)), int(m.group(2))))
    return bad, nlines


def psi6_of(cell):
    rows = []
    f = os.path.join(cell, "speed_of_sound_psi6.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        for _, q in d.iterrows():          # older files have hold/end only (no run_mean ... columns)
            rows.append(dict(M=int(q["wall_mass_factor"]), hold=q.get("psi6_global_hold", np.nan), end=q.get("psi6_global_end", np.nan),
                             mean=q.get("psi6_run_mean", np.nan), sd=q.get("psi6_run_sd", np.nan), lo=q.get("psi6_run_min", np.nan), hi=q.get("psi6_run_max", np.nan)))
        return rows, "speed_of_sound_psi6.csv (per-run summary" + ("" if "psi6_run_mean" in d.columns else ", hold/end only") + ")"
    for lg in sorted(glob.glob(os.path.join(cell, "m_*", "run.log"))):
        M = int(lg.split("m_")[-1].split("/")[0])
        for m in PSI6_LOG.finditer(open(lg, errors="ignore").read()):
            rows.append(dict(M=M, hold=float(m.group(1)), end=float(m.group(2)), mean=float(m.group(3)), sd=float(m.group(4)), lo=np.nan, hi=np.nan))
    return rows, ("run.log 'psi6:' lines (per-run summary)" if rows else "none")


def traj(p):
    """One trajectory: the argmax nu (reduce_B.py's estimator with TD periods: 200 where the trajectory's planned record has 200 of its
    own predicted periods -- the canonical estimator --, else its planned periods rounded down), and the mean-removed position ACF."""
    try:
        h = header(p); per_t = float(h["Planned_Duration"]) * float(h["Predicted_Frequency"])
        TD = 200 if per_t >= 199.999 else int(per_t + 1e-9)
        t, x, nup = T._load(p); dt = (t[-1] - t[0]) / (len(t) - 1)
        n = T._prefix(t, nup, TD)
        if n is None: return None
        P, df = T._spectrum(x[:n], dt); k = int(round(TD / X_EDGE)); nu = (k + int(np.argmax(P[k:]))) * df
        y = x[:n] - x[:n].mean(); nl = int(min(20.0, TD / 2.0) / (nup * dt))
        f = np.fft.rfft(y, 2 * len(y)); a = np.fft.irfft(f * np.conj(f))[:nl + 1] / np.arange(len(y), len(y) - nl - 1, -1)
        return dict(nu=nu, nup=nup, dt=dt, n=n, acf=a.astype(np.float64), T_rec=float(t[n - 1] - t[0]), TD=TD)
    except Exception:
        return None


# ------------------------------------------------------------------------------------------------------------- per cell
def analyse_cell(cell, masses, res, hdr):
    L0, eta_nom, Ns = float(hdr["L0"]), float(hdr["eta"]), int(float(hdr["Left_Count"]))
    delta = box_delta(L0); eta_t = eta_nom * L0 / (L0 - delta / 2); LeT = T.l_eff(L0) - delta / 2
    rows = []
    for M in sorted(masses):
        tr = [q for q in res[M] if q is not None]
        if len(tr) < 3: continue
        al = M / (2.0 * Ns); K = T.k_root(al); x = K / (2 * math.pi * LeT)
        nu = np.array([q["nu"] for q in tr]); dt = tr[0]["dt"]
        L = min(len(q["acf"]) for q in tr); A = np.array([q["acf"][:L] for q in tr])
        warnings.filterwarnings("ignore")
        nd = []
        for q in tr:
            try: nd.append(F.fit_acf(q["acf"], dt, nu.mean())[4] / (2 * math.pi))
            except Exception: nd.append(np.nan)
        nd = np.array(nd); nd = nd[np.isfinite(nd)]
        try:
            pm = R._fit_acf(A.mean(0), dt, nu.mean()); tau = pm[3]
            jk = []
            if len(A) >= 10:
                for g in range(10):
                    keep = [i for i in range(len(A)) if i % 10 != g]
                    try: jk.append(R._fit_acf(A[keep].mean(0), dt, nu.mean())[3])
                    except Exception: pass
            s_tau = math.sqrt((len(jk) - 1) / len(jk) * float(np.sum((np.array(jk) - np.mean(jk)) ** 2))) if len(jk) >= 3 else float("nan")
        except Exception:
            tau, s_tau = float("nan"), float("nan")
        se = nu.std(ddof=1) / math.sqrt(len(nu)); sed = nd.std(ddof=1) / math.sqrt(len(nd)) if len(nd) > 2 else float("nan")
        rows.append(dict(M=M, alpha=al, x=x, n=len(nu), TD=tr[0]["TD"], nu=nu.mean(), se=se, c=nu.mean() / x, sc=se / x,
                         nud=nd.mean() if len(nd) else np.nan, cd=(nd.mean() / x) if len(nd) else np.nan, scd=sed / x,
                         tau=tau, Gam=2 / tau if tau == tau else np.nan, sGam=2 * s_tau / tau ** 2 if tau == tau else np.nan,
                         Q=math.pi * nu.mean() * tau, period=1 / nu.mean(), T_rec=float(np.mean([q["T_rec"] for q in tr]))))
    if len(rows) < 3: return None
    x = np.array([r["x"] for r in rows]); y = np.array([r["nu"] for r in rows]); sy = np.array([r["se"] for r in rows])
    cu = float((x * y).sum() / (x * x).sum()); w = 1 / sy ** 2; cw = float((w * x * y).sum() / (w * x * x).sum())
    chu = float((((y - cu * x) / sy) ** 2).sum()); chw = float((((y - cw * x) / sy) ** 2).sum())
    cu_errs = slope_with_errors(x, y, sy)[2]
    kr = T.kr_cs(eta_t)
    al = np.array([r["alpha"] for r in rows]); c = np.array([r["c"] for r in rows]); cd = np.array([r["cd"] for r in rows])
    sp = spearmanr(al, c); spd = spearmanr(al[np.isfinite(cd)], cd[np.isfinite(cd)])
    return dict(cell=cell, camp=campaign_of(cell), L0=L0, eta_nom=eta_nom, eta=eta_t, Ns=Ns, N=2 * Ns, cu_errs=cu_errs,
                H=Ns * math.pi * 0.25 / (eta_nom * L0), rows=rows, cu=cu, cw=cw, sw=1 / math.sqrt(float((w * x * x).sum())),
                chu=chu, chw=chw, dof=len(rows) - 1, kr=kr, rho=sp.correlation, prho=sp.pvalue, rhod=spd.correlation, prhod=spd.pvalue,
                inside=WIN_LO <= eta_t <= WIN_HI)


def dips(cells):
    """Per campaign: D_M = c_app,M / c_app,M(entry) - 1, entry = the last cell below the window; Spearman and ln(alpha) slope."""
    for camp in sorted({c["camp"] for c in cells}):
        cc = sorted([c for c in cells if c["camp"] == camp], key=lambda c: c["eta"])
        below = [c for c in cc if c["eta"] < WIN_LO]
        if not below: continue
        ent = below[-1]; ref = {r["M"]: r for r in ent["rows"]}
        for c in cc:
            c["entry"] = ent["eta"]; D, sD, al = [], [], []
            for r in c["rows"]:
                if r["M"] not in ref: continue
                q = ref[r["M"]]; v = r["c"] / q["c"]; s_v = v * math.hypot(r["sc"] / r["c"], q["sc"] / q["c"])
                D.append(v - 1); sD.append(s_v); al.append(r["alpha"]); r["D"], r["sD"] = v - 1, s_v
            D, sD, al = np.array(D), np.array(sD), np.array(al)
            if c is ent or len(D) < 3:
                c["rhoD"] = c["pD"] = c["slope"] = c["sslope"] = np.nan; continue
            s = spearmanr(al, D); c["rhoD"], c["pD"] = s.correlation, s.pvalue
            X = np.vstack([np.ones_like(al), np.log(al)]).T; W = np.diag(1 / sD ** 2)
            cov = np.linalg.inv(X.T @ W @ X); b = cov @ X.T @ W @ D
            c["slope"], c["sslope"] = float(b[1]), float(math.sqrt(cov[1, 1]))
            c["chi_lin"] = float((((D - X @ b) / sD) ** 2).sum())


# ------------------------------------------------------------------------------------------------------------- item 3
def kr88(eta):
    import paper1_kr_sanity_261002 as KRS
    return KRS.Z(eta, 0.88), KRS.dZ(eta, 0.88)


def kr90(eta):
    h = 1e-5; a = np.array([eta])
    z = float(SOS.Z_kolafa_rottner_2006(a)[0]); dz = float((SOS.Z_kolafa_rottner_2006(a + h) - SOS.Z_kolafa_rottner_2006(a - h))[0] / (2 * h))
    return z, dz


def hend(eta):
    h = 1e-6; z = float(SOS.Z_henderson_eos(np.array([eta]))[0])
    dz = float((SOS.Z_henderson_eos(np.array([eta + h])) - SOS.Z_henderson_eos(np.array([eta - h])))[0] / (2 * h))
    return z, dz


def cs_of(eta, z, dz):
    v = z + eta * dz + z * z
    return math.sqrt(v) if v > 0 else float("nan")


def zplat(eta, ps): return ps * math.pi / (4 * eta)


def psi6_runs(cell):
    """Per run (M, run index) -> (hold, end, run mean) of psi6: A1 v2 run.log sections ('##RUN ... run r seed s' then 'psi6:'),
    or the per-run CSV (wall_mass_factor, repeat) of the flat campaigns that wrote run means."""
    out = {}
    f = os.path.join(cell, "speed_of_sound_psi6.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        if "psi6_run_mean" in d.columns:
            for _, q in d.iterrows(): out[(int(q["wall_mass_factor"]), int(q["repeat"]))] = (q["psi6_global_hold"], q["psi6_global_end"], q["psi6_run_mean"])
        return out
    for lg in sorted(glob.glob(os.path.join(cell, "m_*", "run.log"))):
        M = int(lg.split("m_")[-1].split("/")[0])
        for sec in open(lg, errors="ignore").read().split("##RUN")[1:]:
            m = re.search(r"run (\d+) seed", sec.split("\n", 1)[0]); q = PSI6_LOG.search(sec)
            if m and q: out[(M, int(m.group(1)))] = (float(q.group(1)), float(q.group(2)), float(q.group(3)))
    return out


def structure(cells, res):
    """Item 4 (information): is the per-mass ordering in the window a frequency effect or a record-length (aging) effect? The record of
    each trajectory is a fixed number of ITS OWN periods, so heavy dividers run 4-5x longer in time while psi6 is still changing."""
    print("\n## Item 4b -- structure against mass and against the measured frequency (EXPLORATORY; per-run psi6 summaries)\n")
    print("| campaign | eta_true | psi6 run mean, M ascending (seed means) | psi6 at the end, M ascending | record, lightest - heaviest [sigma-time] | "
          "Spearman(psi6 run mean, alpha) | Spearman(c_app,M, psi6 run mean) over masses | per trajectory: Spearman(nu/nu_M - 1, psi6 run mean - mean_M), n |\n"
          "|---|---|---|---|---|---|---|---|")
    for c in cells:
        if not (0.69 <= c["eta"] <= 0.725): continue
        pr = psi6_runs(c["cell"])
        if not pr: continue
        pm, pe, x_nu, x_ps = [], [], [], []
        for r in c["rows"]:
            v = [pr[(r["M"], o["r"])] for o in res[c["cell"]][r["M"]] if o is not None and (r["M"], o["r"]) in pr]
            if not v: pm.append(np.nan); pe.append(np.nan); continue
            v = np.array(v, float); pm.append(v[:, 2].mean()); pe.append(v[:, 1].mean())
            nus = [(o["nu"], pr[(r["M"], o["r"])][2]) for o in res[c["cell"]][r["M"]] if o is not None and (r["M"], o["r"]) in pr]
            a = np.array(nus, float); x_nu += list(a[:, 0] / a[:, 0].mean() - 1); x_ps += list(a[:, 1] - a[:, 1].mean())
        al = np.array([r["alpha"] for r in c["rows"]]); cc = np.array([r["c"] for r in c["rows"]]); pm = np.array(pm); pe = np.array(pe)
        ok = np.isfinite(pm)
        if ok.sum() < 3: continue
        s1 = spearmanr(al[ok], pm[ok]); s2 = spearmanr(cc[ok], pm[ok]); s3 = spearmanr(x_nu, x_ps)
        print(f"| {c['camp']} | {c['eta']:.4f} | {', '.join('%.2f' % v for v in pm)} | {', '.join('%.2f' % v for v in pe)} | "
              f"{c['rows'][0]['T_rec']:.0f} - {c['rows'][-1]['T_rec']:.0f} | {s1.correlation:+.2f} (p {s1.pvalue:.2f}) | {s2.correlation:+.2f} (p {s2.pvalue:.2f}) | "
              f"{s3.correlation:+.2f} (p {s3.pvalue:.1e}), {len(x_nu)} |")


def item3(canon):
    print("\n## Item 3 -- equilibrium prediction [DERIVATION from the cited EOS; SOURCE: Engel et al. Eq. (1), Table I; INFERENCE where flagged]\n")
    print("c_s^2 = (kT/m)(Z + eta Z' + Z^2), kT = m = 1. Fluid curves; KR = Kolafa-Rottner 2006: 'ρmax 0.90' is the module (fitted to "
          f"eta {KR_VALID:.4f}, compared with data only to {KR_CMP}); 'ρmax 0.88' is fitted to eta 0.6912 (beyond = extrapolation).\n")
    print("| eta | KR ρmax 0.90: Z | eta Z' | c_s | KR ρmax 0.88: Z | eta Z' | c_s | Henderson: Z | eta Z' | c_s | flag |\n|---|---|---|---|---|---|---|---|---|---|---|")
    for e in (0.66, 0.67, 0.68, 0.69, 0.695, 0.700, 0.7034, 0.7069):
        a, b, h = kr90(e), kr88(e), hend(e)
        flag = "fit range, compared" if e <= KR_CMP else ("ρmax 0.90 fit range, not compared; ρmax 0.88 extrapolated" if e <= KR_VALID else "")
        print(f"| {e:.4f} | {a[0]:.4f} | {e * a[1]:+.3f} | {cs_of(e, *a):.4f} | {b[0]:.4f} | {e * b[1]:+.3f} | {cs_of(e, *b):.4f} | "
              f"{h[0]:.4f} | {e * h[1]:+.3f} | {cs_of(e, *h):.4f} | {flag} |")
    print("\nIdeal coexistence plateau (Z + eta Z' = 0 => c_0^2 = Z^2): c_0 = Z = P* pi / (4 eta)\n")
    print("| eta | Z = c_0, P* = 9.17 | Z = c_0, P* = 9.19 |\n|---|---|---|")
    for e in (0.700, 0.702, 0.704, 0.706, 0.708, 0.710, 0.712, 0.714, 0.716):
        print(f"| {e:.3f} | {zplat(e, 9.17):.4f} | {zplat(e, 9.19):.4f} |")
    e = 0.700; zp = zplat(e, 9.17)
    zh, dzh = hend(e); hyb = math.sqrt(zp + e * (dzh / zh) * zp + zp * zp)
    print("\n### The step at eta = 0.700 (fluid -> plateau) and the variation of c_0 across the window\n")
    print("| fluid reference at 0.700 | Z | eta Z' | Z^2 | c_fluid | c_0 (P* 9.17) | step c_0/c_fluid - 1 [%] | c_0 (P* 9.19) | step [%] |\n|---|---|---|---|---|---|---|---|---|")
    for lab, (z, dz) in (("KR ρmax 0.90 (the module; inside its fit range, Z' already in the loop)", kr90(e)),
                         ("KR ρmax 0.88 (extrapolated from 0.6912: INFERENCE)", kr88(e)), ("Henderson (fluid approximant)", hend(e))):
        cf = cs_of(e, z, dz)
        print(f"| {lab} | {z:.4f} | {e * dz:+.3f} | {z * z:.3f} | {cf:.4f} | {zplat(e, 9.17):.4f} | {100 * (zplat(e, 9.17) / cf - 1):+.2f} | "
              f"{zplat(e, 9.19):.4f} | {100 * (zplat(e, 9.19) / cf - 1):+.2f} |")
    print(f"| plan author's construction: Z = the plateau Z, eta Z' = eta (Z'/Z)_Henderson Z | {zp:.4f} | {e * dzh / zh * zp:+.3f} | {zp * zp:.3f} | "
          f"{hyb:.4f} | {zp:.4f} | {100 * (zp / hyb - 1):+.2f} | {zplat(e, 9.19):.4f} | {100 * (zplat(e, 9.19) / hyb - 1):+.2f} |")
    print(f"\nvariation of c_0 across the window: c_0(0.716)/c_0(0.700) - 1 = {100 * (0.700 / 0.716 - 1):+.2f} % (either P*)")
    print("\n### The N = 100 data against these levels (canonical A1 v2 table, 260919_A1v2_final_cs_vs_eta.csv)\n")
    print("| eta_true | c_s (canonical) | c_s_err_scaled | KR ρmax 0.90 c_s | data/KR - 1 [%] | c_0 (P* 9.17) | data/c_0 - 1 [%] |\n|---|---|---|---|---|---|---|")
    for _, q in canon[(canon["eta"] >= 0.65) & (canon["eta"] <= 0.735)].iterrows():
        e = float(q["eta"]); k = cs_of(e, *kr90(e)) if e <= KR_VALID else float("nan"); c0 = zplat(e, 9.17) if 0.6995 <= e <= 0.7175 else float("nan")
        print(f"| {e:.4f} | {q['c_s']:.4f} | {q['c_s_err_scaled']:.4f} | {k:.4f} | {100 * (q['c_s'] / k - 1):+.1f} | {c0:.4f} | {100 * (q['c_s'] / c0 - 1):+.1f} |")
    win = canon[(canon["eta"] >= 0.69) & (canon["eta"] <= 0.725)]
    cmax = win.loc[win["c_s"].idxmax()]; cmin = win[win["eta"] > cmax["eta"]].loc[lambda d: d["c_s"].idxmin()]
    print(f"\ncanonical N = 100: maximum {cmax['c_s']:.4f} at {cmax['eta']:.4f}, minimum {cmin['c_s']:.4f} at {cmin['eta']:.4f}: "
          f"depth {100 * (1 - cmin['c_s'] / cmax['c_s']):.1f} %; ratio data/c_0 at the minimum {cmin['c_s'] / zplat(cmin['eta'], 9.17):.3f}, "
          f"data/KR(ρmax 0.90) at {cmax['eta']:.4f} {cmax['c_s'] / cs_of(cmax['eta'], *kr90(cmax['eta'])):.3f}")


# ------------------------------------------------------------------------------------------------------------- baseline
def baseline():
    """Test T legacy (pi/8, the plain fluid), four blocks of 25 seeds per mass, and the 279282b campaign anchor (25 seeds)."""
    import resched_testT_261007 as TT
    with contextlib.redirect_stdout(io.StringIO()):
        c = {x["cid"]: x for x in R.cells()}[TT.CID]; co = dict(c); R.method_B(co)
    Ms, x = c["Ms"], np.asarray(co["x"], float); out = []
    nu = {M: pd.read_csv(os.path.join(TT.LOC, TT.REL_T, "legacy", TT.CID, f"m_{M}", "red_nu.csv")).sort_values("run")["nu"].to_numpy(float) for M in Ms}
    sets = [(f"Test T legacy, seeds {25 * b}-{25 * b + 24}", {M: nu[M][25 * b:25 * b + 25] for M in Ms}) for b in range(4)]
    sets.append(("campaign anchor (279282b), 25 seeds", {M: np.array([r for r in co["B"] if r["M"] == M][0]["nus"]) for M in Ms}))
    for lab, v in sets:
        y = np.array([v[M].mean() for M in Ms]); sy = np.array([v[M].std(ddof=1) / math.sqrt(len(v[M])) for M in Ms])
        cu = float((x * y).sum() / (x * x).sum()); w = 1 / sy ** 2; cw = float((w * x * y).sum() / (w * x * x).sum())
        out.append((lab, float((((y - cu * x) / sy) ** 2).sum()), float((((y - cw * x) / sy) ** 2).sum()), float(spearmanr(c["Ms"], y / x).correlation)))
    return out


# ------------------------------------------------------------------------------------------------------------- figures
def setup_mpl():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                         "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 1.6, "figure.dpi": 150})
    return plt


def fig_cell(plt, c, path, pdf):
    al = [r["alpha"] for r in c["rows"]]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.6))
    a1.errorbar(al, [r["c"] for r in c["rows"]], yerr=[r["sc"] for r in c["rows"]], fmt="o", ms=5, color=INK, capsize=2, lw=1, label="argmax nu (registered estimator)")
    a1.errorbar(al, [r["cd"] for r in c["rows"]], yerr=[r["scd"] for r in c["rows"]], fmt="s", ms=4.5, mfc="none", color=BLUE, capsize=2, lw=1, label="nu_d (ACF fit, sec. 4.4.13 4a)")
    a1.axhline(c["cw"], color=MUTED, ls=":", lw=1.2, label=f"weighted single c_s {c['cw']:.3f}")
    if c["eta"] <= KR_VALID:
        a1.axhline(c["kr"], color=ORANGE, ls="-" if c["eta"] <= KR_CMP else "--", label=f"KR fluid {c['kr']:.3f}" + ("" if c["eta"] <= KR_CMP else " (fit range, not compared)"))
    else:
        a1.text(0.02, 0.04, f"KR at eta {c['eta']:.4f}: beyond its fit range (0.7069), not drawn", transform=a1.transAxes, color=MUTED, fontsize=7.5)
    a1.set_xscale("log"); a1.set_xlabel("alpha = M / (2 N_s m)"); a1.set_ylabel("c_app,M = nu_M / x_M  [sqrt(kT/m)]"); a1.legend(fontsize=7.5)
    a2.errorbar(al, [r["Gam"] for r in c["rows"]], yerr=[r["sGam"] for r in c["rows"]], fmt="o", ms=5, color=INK, capsize=2, lw=1)
    a2.set_xscale("log"); a2.set_yscale("log"); a2.set_xlabel("alpha = M / (2 N_s m)"); a2.set_ylabel("Gamma_M = 2 / tau_r  [1/sigma-time]")
    fig.suptitle(f"{c['camp']}  |  N = {c['N']}, H = {c['H']:.2f}, L0 = {c['L0']:.4f}, eta_true = {c['eta']:.4f}"
                 f"{'  (inside the window)' if c['inside'] else ''}  |  EXPLORATORY, post hoc", fontsize=9)
    fig.tight_layout(); fig.savefig(path + ".png")
    if pdf: fig.savefig(path + ".pdf")
    plt.close(fig)


def fig_dips(plt, camp, cc, path, pdf):
    ins = [c for c in cc if c["inside"]]; outs = [c for c in cc if not c["inside"] and c["eta"] != c.get("entry")]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for c in outs:
        ax.plot([r["alpha"] for r in c["rows"] if "D" in r], [100 * r["D"] for r in c["rows"] if "D" in r], "--", color=MUTED, lw=1, marker=".", ms=4)
        r = [q for q in c["rows"] if "D" in q][-1]; ax.annotate(f"{c['eta']:.4f}", (r["alpha"], 100 * r["D"]), xytext=(4, 0), textcoords="offset points", fontsize=7, color=MUTED, va="center")
    cols = RAMP4 if len(ins) <= 4 else None
    for i, c in enumerate(ins):
        col = cols[i] if cols else BLUE
        rr = [q for q in c["rows"] if "D" in q]
        ax.errorbar([q["alpha"] for q in rr], [100 * q["D"] for q in rr], yerr=[100 * q["sD"] for q in rr], color=col, marker="o", ms=4.5, capsize=2, lw=1.6,
                    label=f"eta {c['eta']:.4f}: Spearman {c['rhoD']:+.2f}")
    ax.axhline(0, color=INK, lw=0.8); ax.set_xscale("log"); ax.set_xlabel("alpha = M / (2 N_s m)")
    ax.set_ylabel(f"dip D_M = c_app,M / c_app,M(eta {cc[0].get('entry', float('nan')):.4f}) - 1  [%]")
    ax.set_title(f"{camp}: per-mass dip, window cells (blue, light -> dark = rising eta); grey dashed = outside\nEXPLORATORY, post hoc", fontsize=8.5)
    ax.legend(fontsize=7.5); fig.tight_layout(); fig.savefig(path + ".png")
    if pdf: fig.savefig(path + ".pdf")
    plt.close(fig)


def fig_eta(plt, camp, cc, path, pdf):
    Ms = sorted({r["M"] for c in cc for r in c["rows"]}); n = len(Ms); ncol = 3; nrow = math.ceil(n / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(9.6, 2.3 * nrow), sharex=True, sharey=True); axs = np.atleast_1d(axs).ravel()
    for ax, M in zip(axs, Ms):
        pts = [(c["eta"], r["D"], r["sD"], r["alpha"]) for c in cc for r in c["rows"] if r["M"] == M and "D" in r]
        if not pts: continue
        e, D, sD, al = zip(*sorted(pts))
        ax.axvspan(0.700, 0.716, color=GRID, lw=0); ax.axhline(0, color=MUTED, lw=0.8)
        ax.errorbar(e, 100 * np.array(D), yerr=100 * np.array(sD), color=INK, marker="o", ms=3.5, capsize=1.5, lw=1.2)
        ax.set_title(f"M = {M} (alpha {al[0]:g})", fontsize=8)
    for ax in axs[n:]: ax.axis("off")
    fig.supxlabel("eta_true"); fig.supylabel("dip D_M [%] relative to the entry cell")
    fig.suptitle(f"{camp}: per-mass dip against eta (grey band = Engel's coexistence 0.700-0.716) | EXPLORATORY, post hoc", fontsize=9)
    fig.tight_layout(); fig.savefig(path + ".png")
    if pdf: fig.savefig(path + ".pdf")
    plt.close(fig)


def fig_equilibrium(plt, canon, a1v2, path):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    e = np.linspace(0.66, KR_VALID, 200)
    k90 = np.array([cs_of(q, *kr90(q)) for q in e]); ax.plot(e[e <= KR_CMP], k90[e <= KR_CMP], color=ORANGE, label="KR ρmax 0.90 (fluid fit, compared to 0.69)")
    ax.plot(e[e >= KR_CMP], k90[e >= KR_CMP], color=ORANGE, ls="--", label="KR ρmax 0.90, 0.69-0.7069 (fit range, not compared)")
    eh = np.linspace(0.66, 0.716, 200); ax.plot(eh, [cs_of(q, *hend(q)) for q in eh], color=AQUA, label="Henderson (fluid approximant)")
    ep = np.linspace(0.700, 0.716, 50)
    ax.plot(ep, [zplat(q, 9.17) for q in ep], color=VIOLET, lw=2.2, label="ideal plateau c_0 = Z, P* = 9.17 (Engel)")
    ax.plot(ep, [zplat(q, 9.19) for q in ep], color=VIOLET, lw=1.2, ls="--", label="ideal plateau, P* = 9.19")
    for c in a1v2:
        ax.plot([c["eta"]] * len(c["rows"]), [r["c"] for r in c["rows"]], ".", color=MUTED, ms=3.5, zorder=2)
    ax.plot([], [], ".", color=MUTED, label="N = 100 per-mass c_app,M (A1 v2, item 2)")
    q = canon[(canon["eta"] >= 0.66) & (canon["eta"] <= 0.735)]
    ax.errorbar(q["eta"], q["c_s"], yerr=q["c_s_err_scaled"], fmt="o", ms=4.5, color=INK, capsize=2, lw=1, zorder=3, label="N = 100 canonical c_s (A1 v2)")
    ax.axvspan(0.700, 0.716, color=GRID, lw=0, zorder=0)
    ax.set_xlabel("eta_true"); ax.set_ylabel("sound speed [sqrt(kT/m)]"); ax.set_xlim(0.655, 0.737)
    ax.set_title("Item 3: equilibrium levels against the N = 100 data (grey band = coexistence 0.700-0.716) | EXPLORATORY", fontsize=8.5)
    ax.legend(fontsize=7.2, loc="upper left"); fig.tight_layout(); fig.savefig(path + ".png"); fig.savefig(path + ".pdf"); plt.close(fig)


# ------------------------------------------------------------------------------------------------------------- main
def main():
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 10
    print("# Exploratory: per-mass response in and near the melting window (261012 sec. 4.6) -- POST HOC, decides nothing\n")
    cells = discover(); canon = pd.read_csv(CANON)
    # ---------------------------------------------------------------- item 1
    inv = []; unparsed = 0
    for cell, ms in cells.items():
        M0 = sorted(ms)[0]; mt = meta(cell, sorted(ms[M0])[0][1])
        if mt is None: unparsed += 1; continue
        L0, en, h, per, Ns = mt["L0"], mt["eta"], mt["h"], mt["per"], mt["Ns"]
        et = en * L0 / (L0 - box_delta(L0) / 2)
        if not ETA_LO <= et <= ETA_HI: continue
        bad, nh = bad_runs(cell); ps, psrc = psi6_of(cell); d = date_of(cell)
        inv.append(dict(cell=cell, camp=campaign_of(cell), eta=et, eta_nom=en, L0=L0, Ns=Ns, H=Ns * math.pi * 0.25 / (en * L0) if Ns else float("nan"), ms=ms,
                        seeds=(min(len(v) for v in ms.values()), max(len(v) for v in ms.values())), per=per, bad=bad, nh=nh,
                        psi6=psrc, date=d, build=build_of(cell), hdr=h,
                        analysed=d[:10] >= NEW_CORE and len(ms) >= 3 and mt["complete"] and Ns > 0))
    inv.sort(key=lambda q: (q["camp"], q["eta"]))
    print("## Item 1 -- inventory: cells with eta_true in [0.675, 0.735] under 00_eta_sweep_ROMAN (any build)\n")
    print("| campaign | eta_true | eta (header) | N | H | L0 | masses (count: range) | seeds per mass | record [planned periods] | psi6 saved | "
          "health lines (bad runs) | run date | build | analysed |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for q in inv:
        Ms = sorted(q["ms"])
        print(f"| {q['camp']} | {q['eta']:.4f} | {q['eta_nom']:.4f} | {2 * q['Ns']} | {q['H']:.2f} | {q['L0']:.4f} | {len(Ms)}: {Ms[0]}-{Ms[-1]} | "
              f"{q['seeds'][0]}-{q['seeds'][1]} | {q['per']:.1f} | {q['psi6']}; no time series | {q['nh']} ({len(q['bad'])}) | {q['date']} | {q['build']} | "
              f"{'yes' if q['analysed'] else 'no (before the 2026-08-20 new core)' if q['date'][:10] < NEW_CORE else 'no (< 3 masses or an older trace format)'} |")
    print(f"\n(cells whose traces carry no L0/eta columns, not placed in eta: {unparsed})")
    # ---------------------------------------------------------------- per-trajectory work
    jobs = []
    for q in inv:
        if not q["analysed"]: continue
        for M, runs in q["ms"].items():
            for r, p in sorted(runs):
                if (M, r) in q["bad"] or not T.trace_check(p)[0]: continue
                jobs.append((q["cell"], M, r, p))
    with Pool(workers) as pool:
        out = pool.map(traj, [p for _, _, _, p in jobs], chunksize=8)
    res = defaultdict(lambda: defaultdict(list))
    for (cell, M, r, p), o in zip(jobs, out):
        if o is not None: o["r"] = r
        res[cell][M].append(o)
    cells_a = []
    for q in inv:
        if not q["analysed"]: continue
        a = analyse_cell(q["cell"], q["ms"], res[q["cell"]], q["hdr"])
        if a: a["TD"] = sorted({r["TD"] for r in a["rows"]}); cells_a.append(a)
    dips(cells_a)
    # ---------------------------------------------------------------- estimator gate (A1 v2 = canonical)
    print("\n### Estimator gate: the A1 v2 cells recomputed here against the canonical table (TD = 200, the registered unweighted slope)\n")
    print("| eta_true | c_s here | c_s canonical | rel. difference | per-mass nu identical to paper1_populate_cs_err_20261002.cell |\n|---|---|---|---|---|")
    worst = 0.0
    for c in [c for c in cells_a if c["camp"].startswith("A1v2")]:
        k = canon.iloc[(canon["eta"] - c["eta"]).abs().argmin()]; rd = c["cu"] / k["c_s"] - 1; worst = max(worst, abs(rd))
        same = all(abs(canon_cell((c["eta"], c["L0"], r["M"], T.cell_runs(os.path.join(c["cell"], f"m_{r['M']}"), r["M"])))["nu"] - r["nu"]) < 1e-12 for r in c["rows"])
        print(f"| {c['eta']:.4f} | {c['cu']:.6f} | {k['c_s']:.6f} | {rd:+.1e} | {'yes' if same else '**NO**'} |")
    print(f"\nestimator gate: largest relative difference {worst:.1e} -> {'PASS (the A1 v2 numbers below are the canonical ones)' if worst < 1e-6 else '**FAIL: stop**'}")
    if worst >= 1e-6: sys.exit(1)
    # ---------------------------------------------------------------- item 2
    print("\n## Item 2 -- per-mass apparent sound speed, damping and dip, per cell (EXPLORATORY)\n")
    for c in cells_a:
        tag = "INSIDE the window" if c["inside"] else "outside"
        krf = "" if c["eta"] <= KR_CMP else (" [KR inside its fit range, not compared with data]" if c["eta"] <= KR_VALID else " [KR extrapolated: INFERENCE]")
        tds = f"{c['TD'][0]}" if len(c["TD"]) == 1 else f"{c['TD'][0]}-{c['TD'][-1]} by mass"
        print(f"\n### {c['camp']} | eta_true {c['eta']:.4f} ({tag}) | N {c['N']}, H {c['H']:.2f}, L0 {c['L0']:.4f} | estimator TD = {tds}"
              f"{' (canonical)' if c['TD'] == [200] else ' (exploratory: records shorter than 200 periods)'} | KR c_s {c['kr']:.4f}{krf}\n")
        print("| M | alpha | n | c_app,M (argmax) +- SE | c_app,M/KR - 1 [%] | c_app,M (nu_d) +- SE | dip D_M vs entry [%] | Gamma_M [1/sigma-time] +- SE | Q | period [sigma-time] |\n"
              "|---|---|---|---|---|---|---|---|---|---|")
        for r in c["rows"]:
            dD = f"{100 * r['D']:+.2f} +- {100 * r['sD']:.2f}" if "D" in r else ""
            print(f"| {r['M']} | {r['alpha']:g} | {r['n']} | {r['c']:.4f} +- {r['sc']:.4f} | {100 * (r['c'] / c['kr'] - 1):+.1f} | {r['cd']:.4f} +- {r['scd']:.4f} | {dD} | "
                  f"{r['Gam']:.4g} +- {r['sGam']:.2g} | {r['Q']:.1f} | {r['period']:.2f} |")
        print(f"\nsingle c_s: unweighted {c['cu']:.4f} (chi2 {c['chu']:.1f}, {c['dof']} dof), weighted {c['cw']:.4f} +- {c['sw']:.4f} (chi2 {c['chw']:.1f}); "
              f"Spearman(c_app, alpha) {c['rho']:+.2f} (p {c['prho']:.2f}), with nu_d {c['rhod']:+.2f} (p {c['prhod']:.2f}); "
              + (f"dip vs entry eta {c['entry']:.4f}: Spearman(D, alpha) {c['rhoD']:+.2f} (p {c['pD']:.2f}), slope dD/d ln(alpha) "
                 f"{100 * c['slope']:+.3f} +- {100 * c['sslope']:.3f} %/e-fold (z {c['slope'] / c['sslope']:+.2f}, chi2 of the line {c['chi_lin']:.1f})"
                 if c.get("slope") == c.get("slope") and "slope" in c else "entry cell (D = 0 by construction)" if "entry" in c else "no entry cell below the window"))
    print("\n### Summary per campaign (EXPLORATORY): the mass dependence inside and outside the window\n")
    print("| campaign | eta_true | inside | c_s unweighted | chi2 unweighted (dof) | chi2 weighted | Spearman(c_app, alpha) | Spearman(D_M, alpha) | "
          "dD/d ln alpha [%] (z) |\n|---|---|---|---|---|---|---|---|---|")
    for c in cells_a:
        sl = f"{100 * c['slope']:+.3f} ({c['slope'] / c['sslope']:+.1f})" if c.get("slope") == c.get("slope") and "slope" in c else "-"
        rd = f"{c['rhoD']:+.2f}" if c.get("rhoD") == c.get("rhoD") and "rhoD" in c else "-"
        print(f"| {c['camp']} | {c['eta']:.4f} | {'yes' if c['inside'] else 'no'} | {c['cu']:.4f} | {c['chu']:.1f} ({c['dof']}) | {c['chw']:.1f} | {c['rho']:+.2f} | {rd} | {sl} |")
    print("\n### Pooled over the window cells of each campaign (EXPLORATORY): inside vs outside\n")
    print("| campaign | window cells | mean Spearman(D_M, alpha) inside | outside cells | mean Spearman(D_M, alpha) outside | window cells with dD/d ln alpha "
          "< 0 at z < -2 | with z > +2 | with abs(z) < 2 |\n|---|---|---|---|---|---|---|---|")
    for camp in sorted({c["camp"] for c in cells_a}):
        cc = [c for c in cells_a if c["camp"] == camp and c.get("slope") == c.get("slope") and "slope" in c]
        ins = [c for c in cc if c["inside"]]; outs = [c for c in cc if not c["inside"]]
        z = [c["slope"] / c["sslope"] for c in ins]
        print(f"| {camp} | {len(ins)} | {np.mean([c['rhoD'] for c in ins]) if ins else float('nan'):+.2f} | {len(outs)} | "
              f"{np.mean([c['rhoD'] for c in outs]) if outs else float('nan'):+.2f} | {sum(v < -2 for v in z)} | {sum(v > 2 for v in z)} | {sum(abs(v) < 2 for v in z)} |")
    print("\n### The depth of the dip per campaign, the design's definition (261012 sec. 4.2): D = (c_max - c_min)/c_max, c_min after c_max, "
          "eta_true in [0.685, 0.725], c = the unweighted slope +- c_s_err_scaled (EXPLORATORY; campaigns with >= 4 cells there)\n")
    print("| campaign | N | H | cells | c_s by eta_true | c_max at | c_min (after it) at | D [%] +- | design M1: 14.71 % (N/100)^(-1/2) |\n|---|---|---|---|---|---|---|---|---|")
    for camp in sorted({c["camp"] for c in cells_a}):
        cc = sorted([c for c in cells_a if c["camp"] == camp and 0.685 <= c["eta"] <= 0.725], key=lambda c: c["eta"])
        if len(cc) < 4: continue
        im = int(np.argmax([c["cu"] for c in cc])); after = cc[im + 1:]
        lst = ", ".join(f"{c['eta']:.4f}: {c['cu']:.2f}" for c in cc)
        if not after:
            print(f"| {camp} | {cc[0]['N']} | {cc[0]['H']:.0f} | {len(cc)} | {lst} | {cc[im]['eta']:.4f} | (none after the maximum) | 0 (no max-then-min) | "
                  f"{14.71 * (cc[0]['N'] / 100) ** -0.5:.2f} |"); continue
        jn = int(np.argmin([c["cu"] for c in after])); a, b = cc[im], after[jn]
        D = 1 - b["cu"] / a["cu"]; sD = (b["cu"] / a["cu"]) * math.hypot(a["cu_errs"] / a["cu"], b["cu_errs"] / b["cu"])
        print(f"| {camp} | {a['N']} | {a['H']:.0f} | {len(cc)} | {lst} | {a['eta']:.4f} | {b['eta']:.4f} | {100 * D:.1f} +- {100 * sD:.1f} | "
              f"{14.71 * (a['N'] / 100) ** -0.5:.2f} |")
    print("\n### The Test T plain-fluid baseline at the same seed number (25 per mass), pi/8\n")
    print("| data | chi2 unweighted (8 dof) | chi2 weighted (8 dof) | Spearman(c_app, alpha) |\n|---|---|---|---|")
    for lab, cu, cw, rho in baseline():
        print(f"| {lab} | {cu:.1f} | {cw:.1f} | {rho:+.2f} |")
    structure(cells_a, res)
    # ---------------------------------------------------------------- item 3
    item3(canon)
    # ---------------------------------------------------------------- item 4
    print("\n## Item 4 -- divider periods per mass, and psi6\n")
    print("| campaign | eta_true | N | period per mass [sigma-time], M ascending | record per trajectory [sigma-time], lightest - heaviest |\n|---|---|---|---|---|")
    for c in cells_a:
        if not (0.69 <= c["eta"] <= 0.725): continue
        per = ", ".join("%.1f" % r["period"] for r in c["rows"])
        print(f"| {c['camp']} | {c['eta']:.4f} | {c['N']} | {per} | "
              f"{c['rows'][0]['T_rec']:.0f} - {c['rows'][-1]['T_rec']:.0f} |")
    print("\npsi6 per cell (per-run summaries, averaged over all runs of the cell; no psi6 time series exists in any campaign):\n")
    print("| campaign | eta_true | runs | psi6 at release (hold), mean | at the end, mean | run mean (64 samples), mean | within-run SD, mean | "
          "SD of the run means across runs |\n|---|---|---|---|---|---|---|---|")
    for q in inv:
        if not q["analysed"] or not (0.675 <= q["eta"] <= 0.735): continue
        ps, src = psi6_of(q["cell"])
        if not ps: continue
        d = pd.DataFrame(ps)
        print(f"| {q['camp']} | {q['eta']:.4f} | {len(d)} | {d['hold'].mean():.3f} | {d['end'].mean():.3f} | {d['mean'].mean():.3f} | {d['sd'].mean():.3f} | {d['mean'].std(ddof=1):.3f} |")
    print("\nCost of a time-resolved structural clock per trajectory [DERIVATION: text output, ~18 bytes per number incl. separator]:\n")
    print("| N | record of the heaviest divider [sigma-time] (A1 v2, eta 0.7007) | frames at 4 per sigma-time | global psi6(t) series | "
          "positions (x, y per disk) per frame | positions trace per trajectory |\n|---|---|---|---|---|---|")
    ref = [c for c in cells_a if c["camp"].startswith("A1v2") and abs(c["eta"] - 0.7007) < 0.001]
    Trec = ref[0]["rows"][-1]["T_rec"] if ref else float("nan")
    for N in (100, 400, 900):
        fr = 4 * Trec * math.sqrt(N / 100)          # the record scales with L0 ~ sqrt(N/100) at the design's aspect ratio
        print(f"| {N} | {Trec * math.sqrt(N / 100):.0f} | {fr:.0f} | {fr * 2 * 18 / 1e3:.0f} kB | {N * 2 * 18 / 1e3:.1f} kB | {fr * N * 2 * 18 / 1e6:.0f} MB |")
    print("(record scales with L0, i.e. with sqrt(N/100) at the design's H = 10 sqrt(N/100); 2 frames per sigma-time is the Nyquist minimum "
          "for a 1 sigma-time clock, 4 is used here)")
    # ---------------------------------------------------------------- figures
    os.makedirs(OUT, exist_ok=True); plt = setup_mpl(); paths = []
    for c in cells_a:
        sub = os.path.join(OUT, re.sub(r"[^A-Za-z0-9_]+", "_", c["camp"])); os.makedirs(sub, exist_ok=True)
        p = os.path.join(sub, f"cell_eta{c['eta']:.4f}_N{c['N']}"); fig_cell(plt, c, p, c["camp"].startswith("A1v2")); paths.append(p)
    for camp in sorted({c["camp"] for c in cells_a}):
        cc = sorted([c for c in cells_a if c["camp"] == camp], key=lambda c: c["eta"])
        if not any("D" in r for c in cc for r in c["rows"]): continue
        sub = os.path.join(OUT, re.sub(r"[^A-Za-z0-9_]+", "_", camp)); pdf = camp.startswith("A1v2")
        fig_dips(plt, camp, cc, os.path.join(sub, "summary_dip_vs_alpha"), pdf); fig_eta(plt, camp, cc, os.path.join(sub, "summary_dip_vs_eta"), pdf)
        paths += [os.path.join(sub, "summary_dip_vs_alpha"), os.path.join(sub, "summary_dip_vs_eta")]
    fig_equilibrium(plt, canon, [c for c in cells_a if c["camp"].startswith("A1v2")], os.path.join(OUT, "item3_equilibrium_levels"))
    paths.append(os.path.join(OUT, "item3_equilibrium_levels"))
    print(f"\nfigures: {len(paths)} (PNG; PDF too for A1 v2 and item 3) under {os.path.relpath(OUT, os.path.dirname(HS))}/")


if __name__ == "__main__":
    main()
