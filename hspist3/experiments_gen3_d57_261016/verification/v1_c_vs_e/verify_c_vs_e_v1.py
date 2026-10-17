#!/usr/bin/env python3
"""Independent verification of Report D57-3 task 1 (c against the relaxed strain e). Read-only on the repo; writes only next to this file."""
import glob, math, os, sys, json
import numpy as np, pandas as pd

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage2_261016"))
import p2_lattice as LT

OUT = os.path.dirname(os.path.abspath(__file__))
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016/P2/merged/stage2/gen3")
E1F = os.path.join(HS, "experiments_gen3_p2_261016/report_261016/p2_cells_e1.csv")
CC = os.path.join(HS, "experiments_gen3_d57_261016/c_vs_e")
PIL = os.path.join(HS, "experiments_gen3_d56_261015/data/pilot/P2_N0900_eta0.760000/m_4500")
PILC = os.path.join(HS, "experiments_gen3_d56_261015/c_analysis_b/c_cells.csv")
MASSES = [900, 1800, 4500, 9000, 18000]
CELLS = [f"P2_N0900_eta{eta}_{g}_{p}" for eta in ("0.720000", "0.780000") for g in ("G1", "G2") for p in ("pin0", "pinturned")]

log = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)

# ---------------- (a) e_rec per run and per cell, straight from p2_frames.csv
runs = []; anomalies = []
for cell in CELLS:
    for M in MASSES:
        d = os.path.join(P2, cell, f"m_{M}")
        F = pd.read_csv(os.path.join(d, "p2_frames.csv")); R = pd.read_csv(os.path.join(d, "p2_runs.csv"))
        if F.duplicated(["run", "frame"]).any(): anomalies.append(f"{cell} m_{M}: duplicated (run, frame) rows")
        if sorted(F["run"].unique()) != list(range(10)): anomalies.append(f"{cell} m_{M}: frame runs {sorted(F['run'].unique())}")
        if sorted(R["run"]) != list(range(10)): anomalies.append(f"{cell} m_{M}: p2_runs runs {sorted(R['run'])}")
        for r in range(10):
            f = F[(F["run"] == r) & (F["phase"] == 1)]
            e0 = f["e_0"].to_numpy(float); e1 = f["e_1"].to_numpy(float)
            nnan = int(np.isnan(e0).sum() + np.isnan(e1).sum())
            if nnan: anomalies.append(f"{cell} m_{M} run {r}: {nnan} NaN strains in the record")
            if len(f) == 0: anomalies.append(f"{cell} m_{M} run {r}: no record frames")
            # phase must be contiguous: hold then record
            ph = F[F["run"] == r].sort_values("frame")["phase"].to_numpy()
            if np.any(np.diff(ph) < 0): anomalies.append(f"{cell} m_{M} run {r}: phase goes 1 -> 0")
            m = 0.5 * (e0 + e1)
            rr = R[R["run"] == r].iloc[0]
            runs.append(dict(cell=cell, M=M, run=r, nrec=len(f), nhold=int((F[(F["run"] == r)]["phase"] == 0).sum()),
                             snapshots_col=int(rr["snapshots"]),
                             e_rec=float(np.nanmean(m)), e_rec_0=float(np.nanmean(e0)), e_rec_1=float(np.nanmean(e1)),
                             e_end=0.5 * (rr["end_e_0"] + rr["end_e_1"]), e_hold=0.5 * (rr["hold_e_0"] + rr["hold_e_1"]),
                             lines_y=int((rr["end_line_0"] == 90) + (rr["end_line_1"] == 90)),
                             psi6=0.5 * (rr["psi6_rec_0"] + rr["psi6_rec_1"]), Z=float(rr["Z"]), clean=int(rr["clean"])))
RU = pd.DataFrame(runs); RU.to_csv(os.path.join(OUT, "my_runs.csv"), index=False, float_format="%.10g")
say("anomalies:", anomalies if anomalies else "none")
say("record frames per run: min", RU["nrec"].min(), "max", RU["nrec"].max(), "; hold frames min", RU["nhold"].min(), "max", RU["nhold"].max(),
    "; frames total == snapshots column:", bool(((RU["nrec"] + RU["nhold"]) == RU["snapshots_col"]).all()), "; clean all 1:", bool((RU["clean"] == 1).all()))

E1 = pd.read_csv(E1F).set_index("cell")
cells = {}
for cell in CELLS:
    g = RU[RU["cell"] == cell]; n = len(g); assert n == 50
    cells[cell] = dict(n=n, e=g["e_rec"].mean(), se=g["e_rec"].std(ddof=1) / math.sqrt(n), sd=g["e_rec"].std(ddof=1),
                       e_left=g["e_rec_0"].mean(), e_right=g["e_rec_1"].mean(), e_end=g["e_end"].mean(), e_hold=g["e_hold"].mean(),
                       lines_y=int(g["lines_y"].sum()), psi6=g["psi6"].mean(), Z=g["Z"].mean(),
                       c=float(E1.loc[cell, "E1_c"]), sc=float(E1.loc[cell, "E1_sig"]))

# compare with CC's cells table
CCc = pd.read_csv(os.path.join(CC, "c_vs_e_cells.csv")).set_index("cell")
cmp_rows = []
def rel(a, b):
    if isinstance(a, (int, np.integer)) and isinstance(b, (int, np.integer)): return 0.0 if a == b else math.inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
colmap = dict(e_rec="e", e_rec_se="se", e_rec_sd="sd", e_rec_left="e_left", e_rec_right="e_right", e_end="e_end", e_hold="e_hold",
              lines_y_end="lines_y", psi6="psi6", Z="Z", c="c", sig_c="sc", runs="n")
for cell in CELLS:
    for ccol, mcol in colmap.items():
        a = CCc.loc[cell, ccol]; b = cells[cell][mcol]
        cmp_rows.append(dict(cell=cell, qty=ccol, cc=float(a), mine=float(b), rel=rel(float(a), float(b))))
CM = pd.DataFrame(cmp_rows); CM.to_csv(os.path.join(OUT, "cmp_cells.csv"), index=False, float_format="%.10g")
say("cells: max relative difference over all compared columns =", CM["rel"].max(), "(note: CC's csv has 8 significant digits)")
for cell in CELLS:
    c = cells[cell]
    say(f"  {cell}: e_rec {100*c['e']:+.6f}% (CC {100*CCc.loc[cell,'e_rec']:+.6f}%), SE {100*c['se']:.6f}% (CC {100*CCc.loc[cell,'e_rec_se']:.6f}%), "
        f"sd {100*c['sd']:.6f}%, c {c['c']:.6f} +- {c['sc']:.6f}")

# also compare per-run with CC's c_vs_e_runs.csv
CR = pd.read_csv(os.path.join(CC, "c_vs_e_runs.csv"))
mg = RU.merge(CR, on=["cell", "M", "run"], suffixes=("", "_cc"))
say("runs merged:", len(mg), "of", len(RU), "; max rel diff e_rec", float((abs(mg["e_rec"] - mg["e_rec_cc"]) / abs(mg["e_rec_cc"])).max()),
    "; frames equal:", bool((mg["nrec"] == mg["frames"]).all()))

# ---------------- fits
def wls(e, y, s2):
    """weighted least squares y = a + b e with variances s2 -> a, b, cov"""
    A = np.column_stack([np.ones_like(e), e]); W = np.diag(1.0 / s2)
    cov = np.linalg.inv(A.T @ W @ A); a, b = cov @ (A.T @ W @ y)
    return a, b, cov

def effvar_fit(e, se, c, sc, tol=1e-14, itmax=1000):
    e, se, c, sc = (np.asarray(v, float) for v in (e, se, c, sc))
    y = np.log(c); sy = sc / c; b = 0.0; it = 0
    for it in range(1, itmax + 1):
        s2 = sy ** 2 + b ** 2 * se ** 2
        a, bn, cov = wls(e, y, s2)
        if abs(bn - b) <= tol * max(1.0, abs(bn)): b = bn; break
        b = bn
    s2 = sy ** 2 + b ** 2 * se ** 2; a, b, cov = wls(e, y, s2)
    chi2 = float(((y - a - b * e) ** 2 / s2).sum()); dof = len(e) - 2
    sb = math.sqrt(cov[1, 1])
    return dict(a=a, b=b, sb=sb, chi2=chi2, dof=dof, sb_scaled=sb * math.sqrt(chi2 / dof) if dof > 0 else math.nan,
                sb_scaled_max1=sb * max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else sb, iters=it, c0=math.exp(a), sa=math.sqrt(cov[0, 0]))

def two_point(p1, p2):
    y1, y2 = math.log(p1["c"]), math.log(p2["c"]); sy1, sy2 = p1["sc"] / p1["c"], p2["sc"] / p2["c"]
    de = p2["e"] - p1["e"]; b = (y2 - y1) / de
    sb = math.sqrt((sy1 ** 2 + sy2 ** 2) / de ** 2 + b ** 2 * (p1["se"] ** 2 + p2["se"] ** 2) / de ** 2)
    a = y1 - b * p1["e"]; a2 = y2 - b * p2["e"]
    return dict(a=a, a_check=a2, b=b, sb=sb, c0=math.exp(a))

pt = lambda cell: dict(e=cells[cell]["e"], se=cells[cell]["se"], c=cells[cell]["c"], sc=cells[cell]["sc"], name=cell)
g1a, g1b = pt("P2_N0900_eta0.720000_G1_pin0"), pt("P2_N0900_eta0.720000_G1_pinturned")
w = np.array([1 / g1a["sc"] ** 2, 1 / g1b["sc"] ** 2])
P1 = dict(e=0.5 * (g1a["e"] + g1b["e"]), se=math.sqrt(g1a["se"] ** 2 + g1b["se"] ** 2) / 2,
          c=float((w * np.array([g1a["c"], g1b["c"]])).sum() / w.sum()), sc=float(w.sum() ** -0.5), name="G1 combined")
P2t, P3 = pt("P2_N0900_eta0.720000_G2_pinturned"), pt("P2_N0900_eta0.720000_G2_pin0")
say(f"point 1 (G1 combined): e {P1['e']:.10g} se {P1['se']:.10g} c {P1['c']:.10g} sc {P1['sc']:.10g}")
say(f"point 2 (G2 turned):   e {P2t['e']:.10g} se {P2t['se']:.10g} c {P2t['c']:.10g} sc {P2t['sc']:.10g}")
say(f"point 3 (G2 pin0):     e {P3['e']:.10g} se {P3['se']:.10g} c {P3['c']:.10g} sc {P3['sc']:.10g}")
three = sorted([P1, P2t, P3], key=lambda q: q["e"]); say("order by e:", [q["name"] for q in three])
four = [pt(c) for c in CELLS[:4]]
fits = {}
f = effvar_fit(*[[q[k] for q in three] for k in ("e", "se", "c", "sc")]); fits["eta 0.72, the three boxes (H 30 combined, 29.833, 29.708)"] = f
f = effvar_fit(*[[q[k] for q in four] for k in ("e", "se", "c", "sc")]); fits["eta 0.72, the four cells (check)"] = f
fits["eta 0.72, the segment H 30 combined -> 29.833"] = two_point(three[0], three[1])
fits["eta 0.72, the segment 29.833 -> 29.708"] = two_point(three[1], three[2])
t1, t2 = pt("P2_N0900_eta0.780000_G1_pinturned"), pt("P2_N0900_eta0.780000_G2_pinturned")
fits["eta 0.78, the two turned cells (H 30, 28.667)"] = two_point(t1, t2)
# order-independence check of the two-point formula
alt = two_point(t2, t1); say("two-point order check (0.78): b", alt["b"], "sb", alt["sb"])

CS = pd.read_csv(os.path.join(CC, "c_vs_e_slopes.csv")).set_index("fit")
srows = []
for k, f in fits.items():
    r = CS.loc[k]
    for q_cc, q_me in (("slope", "b"), ("slope_se", "sb"), ("chi2", "chi2"), ("slope_se_scaled", "sb_scaled_max1"), ("ln_c_at_e0", "a"), ("c_at_e0", "c0")):
        if q_me not in f: continue
        srows.append(dict(fit=k, qty=q_cc, cc=float(r[q_cc]), mine=float(f[q_me]), rel=rel(float(r[q_cc]), float(f[q_me])) if r[q_cc] != 0 or f[q_me] != 0 else 0.0))
    say(f"FIT {k}: " + ", ".join(f"{kk}={vv:.10g}" for kk, vv in f.items() if isinstance(vv, (float, int))))
SM = pd.DataFrame(srows); SM.to_csv(os.path.join(OUT, "cmp_slopes.csv"), index=False, float_format="%.10g")
say("slopes: max relative difference =", SM["rel"].max())

# robustness of the 3-point fit: true minimum of the effective-variance chi2 (b in the weights too) and scipy ODR
from scipy.optimize import minimize_scalar
from scipy import odr
e3 = np.array([q["e"] for q in three]); se3 = np.array([q["se"] for q in three]); y3 = np.log([q["c"] for q in three]); sy3 = np.array([q["sc"] / q["c"] for q in three])
def chi2_of_b(b):
    w = 1 / (sy3 ** 2 + b ** 2 * se3 ** 2); a = (w * (y3 - b * e3)).sum() / w.sum(); return float((w * (y3 - a - b * e3) ** 2).sum())
mn = minimize_scalar(chi2_of_b, bracket=(-5, -3, -1)); say(f"3-point: true min of effective-variance chi2: b {mn.x:.6f}, chi2 {mn.fun:.4f}")
o = odr.ODR(odr.RealData(e3, y3, sx=se3, sy=sy3), odr.unilinear, beta0=[2.6, -3.4]).run()
say(f"3-point: scipy ODR: b {o.beta[1]:.6f} +- {o.sd_beta[1]:.6f} (sd_beta is chi2-scaled), unscaled {math.sqrt(o.cov_beta[1,1]):.6f}, chi2 {o.sum_square:.4f}")
# chi2 of the three points against the line through the outer two, and the deviation of the middle point
f2 = two_point(three[0], three[2]); ymid = f2["a"] + f2["b"] * three[1]["e"]
say(f"middle point vs the line through the outer two: ln c dev {y3[1]-ymid:+.6f}, in units of sqrt(sy^2 + b^2 se^2) of the middle point alone: "
    f"{(y3[1]-ymid)/math.sqrt(sy3[1]**2 + f2['b']**2*se3[1]**2):+.3f}")
# chi2 without the e errors (c errors only)
a_, b_, cov_ = wls(e3, y3, sy3 ** 2); say(f"3-point, c errors only: b {b_:.6f} +- {math.sqrt(cov_[1,1]):.6f}, chi2 {float(((y3-a_-b_*e3)**2/sy3**2).sum()):.4f}")
from scipy.stats import chi2 as CHI2
say(f"p-value of chi2 {fits['eta 0.72, the three boxes (H 30 combined, 29.833, 29.708)']['chi2']:.4f} for 1 dof: "
    f"{CHI2.sf(fits['eta 0.72, the three boxes (H 30 combined, 29.833, 29.708)']['chi2'], 1):.3g}")

# ---------------- (g) the pilot
prow = []
for p in sorted(glob.glob(os.path.join(PIL, "snapshots_*_run*.bin.gz"))):
    r = int(p.split("_run")[1].split(".")[0]); S = LT.read_snapshots(p)
    fr1 = [fr for fr in S["frames"] if fr["phase"] == 1]
    ee = [LT.config_strain(fr["xy"], fr["cen"], S["W"], S["H"]) for fr in fr1]
    e0 = np.array([x[0]["e"] for x in ee]); e1 = np.array([x[1]["e"] for x in ee])
    prow.append(dict(run=r, nframes=len(S["frames"]), nrec=len(fr1), nnan=int(np.isnan(e0).sum() + np.isnan(e1).sum()), e_rec=float(np.mean((e0 + e1) / 2)),
                     e0=float(np.mean(e0)), e1=float(np.mean(e1)), line_last=(ee[-1][0]["line"], ee[-1][1]["line"]), W=S["W"], H=S["H"], N=S["N"]))
PR = pd.DataFrame(prow); say("pilot runs:\n" + PR.to_string())
pe = PR["e_rec"].mean(); pse = PR["e_rec"].std(ddof=1) / math.sqrt(len(PR)); psd = PR["e_rec"].std(ddof=1)
PC = pd.read_csv(PILC); pc = float(PC[PC["cell"] == 76]["c"].iloc[0]); pcrows = len(PC[PC["cell"] == 76])
pr = CCc.loc["D56-1 pilot N0900 eta0.76 (H 30, turned 17x27, alpha 5 only)"]
for q, mine in (("e_rec", pe), ("e_rec_se", pse), ("e_rec_sd", psd), ("e_rec_left", PR["e0"].mean()), ("e_rec_right", PR["e1"].mean()), ("c", pc), ("runs", len(PR))):
    say(f"pilot {q}: CC {float(pr[q]):.10g} mine {mine:.10g} rel {rel(float(pr[q]), float(mine)):.3g}")
say(f"pilot c rows for cell 76 in c_cells.csv: {pcrows} (values {PC[PC['cell']==76]['c'].tolist()})")
json.dump(dict(log=log), open(os.path.join(OUT, "verify_log.json"), "w"), indent=1)
