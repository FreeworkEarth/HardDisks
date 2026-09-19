#!/usr/bin/env python3
"""##CHRIS 2026-09-08 -- TASK 5 analysis for the pressure validation campaign.
Parts A (numerical validity), B (Z_pair vs KR, 1/sqrt(N) extrapolation),
C (wall route: isotropy, wall-pair gap, Z_wall_inf). Part D (structural
boundary / ladder) is a separate script because it needs the ladder blocks.
Usage: analyze_pressure_campaign.py CAMPAIGN_DIR  -> CAMPAIGN_DIR/analysis/
All numbers come from the CSVs on disk. KR (Kolafa-Rottner 2006) is quoted only
for eta <= 0.69: above ~0.695 the fit turns over and is not a reference."""
import sys, os, glob, csv, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import plot_speed_of_sound_edmd as sos
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

R = sys.argv[1]; OUT = os.path.join(R, "analysis"); os.makedirs(OUT, exist_ok=True)
COLS = "eta N seed boxW boxH eq mt pair_ev wL wR wB wT Zp Zp_sem Zx Zx_sem Zy Zy_sem T psi6g psi6l nb fa orep crep wovd valid".split()
KR_MAX = 0.69

rows = []
# ##CHRIS: the N=2500 cells live in their own subdirectory (separate driver,
# separate calibration table); include them so the fit has four sizes.
for f in (glob.glob(f"{R}/preserved/traj_*.csv") + glob.glob(f"{R}/runs/*_traj_*.csv")
          + glob.glob(f"{R}/N2500_*/runs/*_traj_*.csv")):
    for l in open(f):
        if l.strip():
            r = dict(zip(COLS, l.strip().split(","))); r["_file"] = os.path.basename(f); rows.append(r)
for r in rows:
    for k in ("eta","boxW","boxH","Zp","Zp_sem","Zx","Zx_sem","Zy","Zy_sem","T","psi6g","psi6l"): r[k] = float(r[k])
    for k in ("N","nb","fa","orep","crep","wovd","valid","pair_ev","wL","wR","wB","wT"): r[k] = int(float(r[k]))
    r["Zw"] = 0.5*(r["Zx"]+r["Zy"])
    r["diag"] = r["_file"].startswith("diag_")
cal = {}
for r in csv.DictReader(open(f"{R}/chunk_calibration.csv")):
    cal[(round(float(r["eta"]),3), int(r["N"]))] = r
for cf in glob.glob(f"{R}/N2500_*/chunk_calibration_N2500.csv"):
    for r in csv.DictReader(open(cf)):
        cal[(round(float(r["eta"]),3), int(r["N"]))] = r
def KR(eta): return float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0]) if eta <= KR_MAX+1e-9 else float("nan")

cells = {}
for r in rows: cells.setdefault((round(r["eta"],3), r["N"]), []).append(r)
def stat(g, key, semkey):
    v = np.array([x[key] for x in g]); s = np.array([x[semkey] for x in g])
    m = v.mean(); sd = v.std(ddof=1) if len(v) > 1 else float("nan")
    bsem = math.sqrt((s**2).mean()/len(v))        # block-sem combined over seeds
    ssem = sd/math.sqrt(len(v)) if len(v) > 1 else float("nan")
    return m, sd, bsem, ssem, max(bsem, ssem if ssem==ssem else 0.0)

md = []
# ---------------- A ----------------
md.append("## A. Numerical validity\n")
md.append(f"Accepted trajectories on disk: **{len(rows)}** (all valid=1: {all(r['valid']==1 for r in rows)}; "
          f"any nonzero health counter: {sum(1 for r in rows if r['fa'] or r['orep'] or r['crep'] or r['wovd'])}).\n")
md.append("| eta | N | seeds | chunk (prod) | calib A / B | T_mean range | sum health |\n|---|---|---|---|---|---|---|")
for (e,N),g in sorted(cells.items()):
    c = cal.get((e,N)); T = [x["T"] for x in g]
    h = sum(x["fa"]+x["orep"]+x["crep"]+x["wovd"] for x in g)
    md.append(f"| {e:.3f} | {N} | {len(g)} | {c['production_chunk'] if c else '—'} | "
              f"{(c['safe_chunk_seedA']+' / '+c['safe_chunk_seedB']) if c else '—'} | {min(T):.3f}–{max(T):.3f} | {h} |")
md.append("\nT_mean varies between seeds (velocities are drawn from a unit Gaussian and not rescaled; relative "
          "scatter ~1/sqrt(N)). Z = P/(rho kB T) is formed with the measured KE per block, so it is T-independent.\n")

# ---------------- B ----------------
md.append("\n## B. Z_pair vs Kolafa–Rottner 2006, and N -> inf extrapolation\n")
md.append("| eta | N | seeds | Z_pair | seed sd | block sem | KR | dev |\n|---|---|---|---|---|---|---|---|")
Binf = {}
for (e,N),g in sorted(cells.items()):
    m,sd,bsem,ssem,_ = stat(g,"Zp","Zp_sem"); kr = KR(e)
    dev = f"{100*(m-kr)/kr:+.2f}%" if kr==kr else "— (KR invalid)"
    md.append(f"| {e:.3f} | {N} | {len(g)} | {m:.4f} | {sd:.4f} | {bsem:.4f} | {kr:.4f} | {dev} |" if kr==kr else
              f"| {e:.3f} | {N} | {len(g)} | {m:.4f} | {sd:.4f} | {bsem:.4f} | — | {dev} |")
def fit_inf_pow(e, key, semkey, power):
    """##CHRIS: same weighted straight line but in 1/N**power, so the a+b/sqrt(N)
    and a+b/N forms can be compared. Above eta ~ 0.1 the two intercepts differ by
    more than their statistical errors, and that spread is a genuine model
    (form) uncertainty that the claim has to carry."""
    pts = [(N,)+stat(g,key,semkey) for (ee,N),g in cells.items() if ee==e]
    pts = [(N,m,err) for (N,m,sd,bsem,ssem,err) in pts if err>0]
    if len(pts) < 3: return None
    x = np.array([N**(-power) for N,_,_ in pts]); y = np.array([m for _,m,_ in pts])
    w = 1/np.array([err for *_,err in pts])**2
    A = np.vstack([np.ones_like(x), x]).T; W = np.diag(w)
    cov = np.linalg.inv(A.T@W@A); coef = cov@A.T@W@y
    chi2 = float(((y - A@coef)**2 * w).sum())
    return coef[0], math.sqrt(cov[0,0]), coef[1], len(pts), chi2, len(pts)-2

def fit_inf(e, key, semkey, only=None):
    """Weighted straight line Z = Z_inf + a/sqrt(N). Returns
    (Z_inf, sigma, a, n_points, chi2, dof). `only` restricts the N used."""
    pts = [(N,)+stat(g,key,semkey) for (ee,N),g in cells.items() if ee==e and (only is None or N in only)]
    pts = [(N,m,err) for (N,m,sd,bsem,ssem,err) in pts if err>0]
    if len(pts) < 2: return None
    x = np.array([1/math.sqrt(N) for N,_,_ in pts]); y = np.array([m for _,m,_ in pts]); w = 1/np.array([err for *_,err in pts])**2
    A = np.vstack([np.ones_like(x), x]).T; W = np.diag(w)
    cov = np.linalg.inv(A.T@W@A); coef = cov@A.T@W@y
    chi2 = float(((y - A@coef)**2 * w).sum()); dof = len(pts)-2
    return coef[0], math.sqrt(cov[0,0]), coef[1], len(pts), chi2, dof
def chi2_sf(chi2, dof):
    """Survival function of chi-square (dof = 1 or 2 only needed here)."""
    if dof <= 0: return float("nan")
    if dof == 1: return math.erfc(math.sqrt(chi2/2.0))
    if dof == 2: return math.exp(-chi2/2.0)
    # generic series fallback (regularized upper incomplete gamma) for small dof
    k = dof/2.0; x = chi2/2.0; term = 1.0; ssum = 1.0
    for n in range(1, 200):
        term *= x/(k+n); ssum += term
    return min(1.0, max(0.0, 1.0 - math.exp(-x)*x**k/math.gamma(k+1)*ssum))
md.append("\nFit Z(N) = Z_inf + a/sqrt(N) (weighted, error = max(block sem, seed sem)), per eta with >= 2 N. "
          "chi2 has (#N - 2) degrees of freedom; p = P(chi2_dof >= observed). "
          "An eta with chi2 > 4 at 1 dof (p < 0.046) is flagged: the 1/sqrt(N) form is rejected there and Z_inf is not a bulk value.\n")
md.append("| eta | #N | Z_pair_inf | sigma | slope a | chi2 | dof | p | Z_inf from N=900/1600 only | KR | Z_inf vs KR | form |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
FORMROWS=[]
for e in sorted({e for e,_ in cells}):
    f = fit_inf(e,"Zp","Zp_sem")
    if not f: continue
    zi,sg,a,n,chi2,dof = f; kr = KR(e); Binf[e]=(zi,sg)
    f2 = fit_inf(e,"Zp","Zp_sem",only={900,1600})
    z2 = f"{f2[0]:.4f} (a={f2[2]:+.2f})" if f2 else "—"
    p = chi2_sf(chi2,dof) if dof>0 else float("nan")
    form = "n/a (2 N)" if dof < 1 else ("rejected (chi2>4)" if chi2 > 4 else "ok")
    dev = f"{100*(zi-kr)/kr:+.2f}%" if kr==kr else "—"
    g1=fit_inf_pow(e,"Zp","Zp_sem",0.5); g2=fit_inf_pow(e,"Zp","Zp_sem",1.0)
    if g1 and g2:
        spread=abs(g1[0]-g2[0]); comb=math.hypot(g1[1],g2[1])
        FORMROWS.append((e,g1[0],g1[1],g1[4],g2[0],g2[1],g2[4],spread,comb,kr))
    md.append(f"| {e:.3f} | {n} | {zi:.4f} | {sg:.4f} | {a:+.3f} | {chi2:.2f} | {dof} | {p:.3f} | {z2} | {kr:.4f} | {dev} | {form} |" if kr==kr else
              f"| {e:.3f} | {n} | {zi:.4f} | {sg:.4f} | {a:+.3f} | {chi2:.2f} | {dof} | {p:.3f} | {z2} | — | — | {form} |")

md.append("\n### Extrapolation-form comparison: Z_inf from a + b/sqrt(N) vs a + b/N\n")
md.append("Both are weighted straight lines through the same per-cell means. Where the two "
          "intercepts differ by more than their combined statistical error, the choice of form "
          "is the dominant uncertainty and the claim must carry that spread.\n")
md.append("| eta | Z_inf (1/sqrtN) | sigma | chi2 | Z_inf (1/N) | sigma | chi2 | spread | comb. sigma | dev sqrtN | dev 1/N | form dominates |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")
for (e,z1,s1,c1,z2,s2,c2,sp,cb,kr) in FORMROWS:
    d1=f"{100*(z1-kr)/kr:+.2f}%" if kr==kr else "—"
    d2=f"{100*(z2-kr)/kr:+.2f}%" if kr==kr else "—"
    md.append(f"| {e:.3f} | {z1:.4f} | {s1:.4f} | {c1:.2f} | {z2:.4f} | {s2:.4f} | {c2:.2f} | "
              f"{sp:.4f} | {cb:.4f} | {d1} | {d2} | {'**yes**' if sp>cb else 'no'} |")

# ---------------- C ----------------
md.append("\n## C. Wall route: isotropy, wall–pair gap, Z_wall_inf\n")
md.append("Note: the pair virial and the wall momentum flux are tied by momentum balance within the same trajectories; "
          "their agreement is an internal consistency check of the pressure bookkeeping, not an independent measurement of Z. "
          "Z_wall_x vs Z_wall_y is the isotropy check of the wall route.\n")
md.append("| eta | N | Z_wall_x | Z_wall_y | (x−y)/sem | Z_wall | wall−pair | gap ratio vs N=400 | sqrt(400/N) |\n|---|---|---|---|---|---|---|---|---|")
gap400 = {}
for (e,N),g in sorted(cells.items()):
    mx,_,bx,sx,ex = stat(g,"Zx","Zx_sem"); my,_,by,sy,ey = stat(g,"Zy","Zy_sem")
    mp = stat(g,"Zp","Zp_sem")[0]; mw = 0.5*(mx+my); gap = 100*(mw-mp)/mp
    if N==400: gap400[e]=gap
    ratio = f"{gap/gap400[e]:.2f}" if e in gap400 and gap400[e] else "—"
    md.append(f"| {e:.3f} | {N} | {mx:.4f} | {my:.4f} | {(mx-my)/math.hypot(ex,ey):+.2f} | {mw:.4f} | {gap:+.2f}% | {ratio} | {math.sqrt(400/N):.2f} |")
md.append("\n| eta | Z_wall_inf | sigma | Z_pair_inf | wall−pair at inf | KR | Z_wall_inf vs KR |\n|---|---|---|---|---|---|---|")
for e in sorted({e for e,_ in cells}):
    fw = fit_inf(e,"Zw","Zx_sem")
    if not fw or e not in Binf: continue
    zw,sw = fw[0],fw[1]; zp,sp = Binf[e]; kr = KR(e)
    md.append(f"| {e:.3f} | {zw:.4f} | {sw:.4f} | {zp:.4f} | {100*(zw-zp)/zp:+.2f}% | {kr:.4f} | {100*(zw-kr)/kr:+.2f}% |" if kr==kr else
              f"| {e:.3f} | {zw:.4f} | {sw:.4f} | {zp:.4f} | {100*(zw-zp)/zp:+.2f}% | — | — |")
open(f"{OUT}/tables_ABC.md","w").write("\n".join(md)+"\n")

# ---------------- plots ----------------
etas = sorted({e for e,_ in cells if KR(e)==KR(e)})
fig, ax = plt.subplots(figsize=(8.6,5.6)); cmap = plt.cm.viridis(np.linspace(0.05,0.9,len(etas)))
for c,e in zip(cmap,etas):
    kr = KR(e); pts = sorted((N,)+stat(g,"Zp","Zp_sem")[:1]+(stat(g,"Zp","Zp_sem")[4],) for (ee,N),g in cells.items() if ee==e)
    x=[1/math.sqrt(N) for N,_,_ in pts]; y=[100*(m-kr)/kr for _,m,_ in pts]; ye=[100*er/kr for *_,er in pts]
    ax.errorbar(x,y,yerr=ye,fmt="o",color=c,ms=5,capsize=3,label=f"η={e:.3f}")
    f = fit_inf(e,"Zp","Zp_sem")
    if f:
        xs=np.linspace(0,max(x)*1.05,20); ax.plot(xs,100*((f[0]+f[2]*xs)-kr)/kr,"-",color=c,lw=1.2,alpha=0.8)
        ax.plot([0],[100*(f[0]-kr)/kr],"*",color=c,ms=12)
ax.axhline(0,color="k",lw=0.8); ax.axvline(0,color="k",lw=0.6)
ax.set_xlabel(r"$1/\sqrt{N}$   (N = 1600, 900, 400 →)"); ax.set_ylabel(r"$Z_{\rm pair}$ deviation from Kolafa–Rottner  [%]")
ax.set_title("Pair-virial route: finite-size deviation and $1/\\sqrt{N}$ extrapolation (★)"); ax.grid(True,ls=":",alpha=0.6)
ax.legend(fontsize=7,ncol=3); fig.tight_layout(); fig.savefig(f"{OUT}/dev_vs_invsqrtN.pdf"); fig.savefig(f"{OUT}/dev_vs_invsqrtN.png",dpi=200)

fig, ax = plt.subplots(figsize=(8.6,5.6))
eg = np.linspace(0.005,0.69,400); ax.plot(eg, sos.Z_kolafa_rottner_2006(eg), "-", color="C3", lw=2, label="Kolafa–Rottner 2006 (clipped at η=0.69)")
for N,mk,cl in ((400,"s","0.6"),(900,"^","0.35"),(1600,"o","k")):
    pe=[e for (e,n) in sorted(cells) if n==N]; 
    ax.plot(pe,[stat(cells[(e,N)],"Zp","Zp_sem")[0] for e in pe],mk,color=cl,ms=4,alpha=0.8,label=f"Z_pair N={N}")
    ax.plot(pe,[0.5*(stat(cells[(e,N)],"Zx","Zx_sem")[0]+stat(cells[(e,N)],"Zy","Zy_sem")[0]) for e in pe],mk,mfc="none",color="C0",ms=4,alpha=0.8,label=f"Z_wall N={N}")
be=sorted(Binf); ax.plot(be,[Binf[e][0] for e in be],"*",color="C1",ms=13,label=r"$Z_{\rm pair,\infty}$ (1/√N fit)")
ax.set_xlabel(r"packing fraction $\eta$"); ax.set_ylabel(r"$Z = P/(\rho k_B T)$"); ax.set_yscale("log")
ax.set_title("Equilibrium EOS from EDMD: two independent routes"); ax.grid(True,ls=":",alpha=0.6); ax.legend(fontsize=7,ncol=2)
fig.tight_layout(); fig.savefig(f"{OUT}/Z_vs_eta.pdf"); fig.savefig(f"{OUT}/Z_vs_eta.png",dpi=200)
print(f"wrote {OUT}/tables_ABC.md, dev_vs_invsqrtN.pdf, Z_vs_eta.pdf  ({len(rows)} trajectories)")
