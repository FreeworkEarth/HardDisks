#!/usr/bin/env python3
"""##CHRIS: GPT handoff item 4 -- fixed-aspect family A (finitesize_aspect_20260826):
same manifest treatment as the route-A refit (eligibility, health classes, exact
cot K, through-origin and free intercept) per (eta, N) cell, then per eta the
finite-size fit c_s(N) = c_inf + a/sqrt(N) with chi2, the 900/1600-only
intercept, and the comparison to KR (<= 0.69) and Liu. Analysis only."""
import csv, glob, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos
F = sys.argv[1]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
R = 0.5; KR_MAX = 0.69

def chi2_sf(c, dof):
    if dof <= 0: return float("nan")
    if dof == 1: return math.erfc(math.sqrt(c/2))
    if dof == 2: return math.exp(-c/2)
    return float("nan")

man, cells = [], {}
for cell in sorted(glob.glob(f"{F}/eta_*/N*")):
    eta_tag = os.path.basename(os.path.dirname(cell)); N = int(os.path.basename(cell)[1:])
    ana = glob.glob(f"{cell}/an/analysis/*/speed_of_sound_runs.csv")
    if len(ana) != 1: continue
    runs = list(csv.DictReader(open(ana[0])))
    seed_of, health = {}, {}
    for ml in glob.glob(f"{cell}/m_*/run.log"):
        for ln in open(ml):
            m = re.search(r"Running: L0 = ([\d.]+), M = (\d+)\*m, run = (\d+), seed = (\d+)", ln)
            if m: seed_of[(int(m.group(2)), int(m.group(3)))] = int(m.group(4))
            h = re.search(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): forced_advance=(\d+) wall_clamp_repairs=(\d+) overlap_repairs=(\d+) wall_overdue=(\d+)", ln)
            if h: health[(int(h.group(1)), int(h.group(2)))] = tuple(int(h.group(i)) for i in (4,5,6,7))
    L0 = float(runs[0]["L0"]); H = 10.0*math.sqrt(N/100); eta = N*math.pi*R*R/(2*L0*H)
    N_side = N//2; per_mass = {}
    for r in runs:
        M, run = int(float(r["wall_mass_factor"])), int(float(r["run"]))
        fa, cr, orp, wo = health.get((M, run), (0,0,0,0))
        qual = r["frequency_quality_pass"] in ("1","True","true"); bnd = r["peak_on_search_boundary"] in ("1","True","true")
        reasons = []
        if not qual: reasons.append("quality_fail")
        if bnd: reasons.append("peak_on_search_boundary")
        if fa or cr or orp: reasons.append(f"health fa={fa} crep={cr} orep={orp}")
        if wo and (fa or cr or orp): reasons.append(f"wall_overdue={wo} with other health")
        elig = not reasons; nu = float(r["nu"]) if r["nu"] not in ("","nan") else float("nan")
        man.append(dict(family="A", eta=f"{eta:.6f}", N=N, L0=L0, H=H, M=M, run=run, seed=seed_of.get((M,run)), nu=r["nu"],
                        peak_on_search_boundary=int(bnd), frequency_quality_pass=int(qual),
                        forced_advance=fa, clamp_repair=cr, overlap_repair=orp, wall_overdue=wo,
                        eligible=int(elig), exclusion_reasons=";".join(reasons)))
        if elig and nu == nu: per_mass.setdefault(M, []).append(nu)
    xs, ys, ws, nused = [], [], [], 0; L_eff = L0 - 2*R
    for M, nus in sorted(per_mass.items()):
        if len(nus) < 3: continue
        K = sos.k_root_bisect(M/(2.0*N_side))
        xs.append(K/(2*math.pi*L_eff)); ys.append(float(np.mean(nus))); ws.append(max(float(np.std(nus,ddof=1))/math.sqrt(len(nus)),1e-12)); nused += len(nus)
    if len(xs) < 3: continue
    x,y,s = map(np.array,(xs,ys,ws))
    cs, cse = sos.weighted_linreg(x,y,s,force_zero_intercept=True)[:2]
    cf, cfe, b0, b0e = sos.weighted_linreg(x,y,s,force_zero_intercept=False)[:4]
    cells[(round(eta,3), N)] = dict(eta=eta, N=N, L0=L0, H=H, cs=cs, cse=cse, cf=cf, cfe=cfe, b0=b0,
                                    n_masses=len(xs), n_all=len(runs), n_elig=sum(1 for m_ in man if m_["N"]==N and abs(float(m_["eta"])-eta)<1e-6 and m_["eligible"]), n_used=nused)

with open(f"{OUT}/famA_fit_input_manifest.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(man[0].keys())); w.writeheader(); w.writerows(man)

def eos(eta):
    a = np.array([eta]); Zl = sos.Z_liu_global(a); dZl = sos.dZ_liu_global(a)
    cl = float(sos.cs_adiabatic_2d_monatomic(Zl,dZl,a,kbt=1,m=1)[0])
    if eta <= KR_MAX:
        h=1e-5; Zk = sos.Z_kolafa_rottner_2006(a); dZk = (sos.Z_kolafa_rottner_2006(a+h)-sos.Z_kolafa_rottner_2006(a-h))/(2*h)
        ck = float(sos.cs_adiabatic_2d_monatomic(Zk,dZk,a,kbt=1,m=1)[0])
    else: ck = float("nan")
    return ck, cl

def fit(pts):
    x = np.array([1/math.sqrt(N) for N,_,_ in pts]); y = np.array([c for _,c,_ in pts]); w = 1/np.array([e for *_,e in pts])**2
    A = np.vstack([np.ones_like(x),x]).T; W = np.diag(w); cov = np.linalg.inv(A.T@W@A); coef = cov@A.T@W@y
    chi2 = float(((y-A@coef)**2*w).sum()); return coef[0], math.sqrt(cov[0,0]), coef[1], chi2, len(pts)-2

L = ["## famA: per-cell fits (through origin; free-intercept in brackets)\n",
     "| eta | N | L0 | H | eligible/all | masses | c_s (origin) | ± | c_s (free) | intercept |", "|---|---|---|---|---|---|---|---|---|---|"]
for (e,N),c in sorted(cells.items()):
    L.append(f"| {c['eta']:.4f} | {N} | {c['L0']:.3f} | {c['H']:.2f} | {c['n_elig']}/{c['n_all']} | {c['n_masses']} | {c['cs']:.4f} | {c['cse']:.4f} | {c['cf']:.4f} | {c['b0']:+.5f} |")
L += ["\n## famA: finite-size fit per eta, c_s(N) = c_inf + a/sqrt(N) (weights 1/err^2, err = fit error of the through-origin c_s)\n",
      "| eta | #N | c_inf | sigma | a | chi2 | dof | p | c_inf (900/1600 only) | KR (<=0.69) | dev KR | Liu | dev Liu | form |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
summ = []
for e in sorted({e for e,_ in cells}):
    pts = [(N, cells[(e,N)]['cs'], cells[(e,N)]['cse']) for (ee,N) in cells if ee==e]
    if len(pts) < 3: 
        L.append(f"| {e:.3f} | {len(pts)} | — | — | — | — | — | — | — | — | — | — | — | n/a (<3 N) |"); continue
    ci,sg,a,chi2,dof = fit(sorted(pts)); p = chi2_sf(chi2,dof)
    p2 = [q for q in pts if q[0] in (900,1600)]; c2 = fit(p2)[0] if len(p2)==2 else float("nan")
    ck, cl = eos(e); form = "rejected (chi2>4)" if (dof==1 and chi2>4) else ("rejected (chi2>6)" if (dof==2 and chi2>6) else "ok")
    L.append(f"| {e:.3f} | {len(pts)} | {ci:.4f} | {sg:.4f} | {a:+.3f} | {chi2:.2f} | {dof} | {p:.3f} | {c2:.4f} | {ck:.4f} | {100*(ci-ck)/ck:+.2f}% | {cl:.4f} | {100*(ci-cl)/cl:+.2f}% | {form} |" if ck==ck else
             f"| {e:.3f} | {len(pts)} | {ci:.4f} | {sg:.4f} | {a:+.3f} | {chi2:.2f} | {dof} | {p:.3f} | {c2:.4f} | — | — | {cl:.4f} | {100*(ci-cl)/cl:+.2f}% | {form} |")
    summ.append(dict(eta=f"{e:.4f}", n_N=len(pts), c_s_inf=f"{ci:.5f}", sigma=f"{sg:.5f}", a=f"{a:.4f}", chi2=f"{chi2:.3f}", dof=dof, p=f"{p:.4f}", c_s_inf_900_1600=f"{c2:.5f}", cs_KR=f"{ck:.5f}" if ck==ck else "", cs_Liu=f"{cl:.5f}", form=form))
open(f"{OUT}/famA_tables.md","w").write("\n".join(L)+"\n")
with open(f"{OUT}/famA_cs_inf_vs_eta.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
print(f"manifest rows {len(man)}  eligible {sum(m['eligible'] for m in man)}  cells {len(cells)}  eta with >=3 N: {len(summ)}")
print("exclusions:", {k: sum(1 for m in man if k in m['exclusion_reasons']) for k in ("quality_fail","peak_on_search_boundary","health fa","wall_overdue")})
