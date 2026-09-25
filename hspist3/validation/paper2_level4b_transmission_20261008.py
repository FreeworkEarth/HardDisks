#!/usr/bin/env python3
"""##CHRIS 2026-10-08: Level 4b. Applies 261008 sections 1A-1D exactly as committed (9488973,
c6191b5, f510014) BEFORE any 4b trace was analysed. Analysis only.

The observable is section 1D's EXCESS ledger, not Delta E_2/W_in: the reversible P h X term is a
factor 8 larger than the acoustic transfer at M = 200 and would have killed T(x) for the wrong
reason. X_i = Delta KE_i - E_qs,i(L_i(t)), with L_i from the recorded SegEtas and E_qs from the
same KR isentrope as W_qs. The divider mode lives entirely inside E_qs and cancels.
"""
import glob,json,math,os,sys
import numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.optimize import curve_fit
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO=os.path.dirname(os.path.dirname(HERE))
P=os.path.join(REPO,"hspist3","experiments_energy_transfer","level4b_transmission_20261008b")

NS,H,RD,LC,D0 = 50,10.0,0.5,38.75,3.875
Zf =lambda e: float(sos.Z_kolafa_rottner_2006(np.array([e]))[0])
dZf=lambda e: float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
ETA0 = NS*math.pi*RD*RD/(LC*H)
CS   = math.sqrt(Zf(ETA0)+ETA0*dZf(ETA0)+Zf(ETA0)**2)
ZAC  = NS*1.0*CS/LC
TRANSIT, ROUND = LC/CS, 2*LC/CS
W_KR, W_ID = 7.0715, 5.5556
SLOPE_KR, SLOPE_ID = ZAC*D0, (NS*CS/LC if False else NS*math.sqrt(2.0)/LC)*D0

# isentrope T_ad(L) on a grid, interpolated
_Lg=np.linspace(20.0,55.0,1400)
_Tg=np.array([math.exp(quad(lambda e: Zf(e)/e, ETA0, NS*math.pi*RD*RD/(L*H))[0]) for L in _Lg])
Tad=lambda L: np.interp(L,_Lg,_Tg)
Eqs=lambda L: NS*(Tad(L)-1.0)                       # signed energy change of one compartment

def Tspec(x): return 1.0-x*(1.0-np.exp(-1.0/x))     # section 1C
def Tspec_Z(x_over_Z, Z):                           # Z-free form for fitting
    x=x_over_Z/Z; return 1.0-x*(1.0-np.exp(-1.0/x))

CELLS=[(10,0.05),(10,0.2),(10,0.5),(10,1.0),(50,0.05),(50,0.2),(50,0.5),(50,1.0),
       (200,0.05),(200,0.2),(200,0.5),(200,1.0)]
tag=lambda M,u: f"M{M}_u{str(u).replace('0.','').replace('.','').ljust(3,'0')}" if False else \
    {0.05:"u005",0.2:"u020",0.5:"u050",1.0:"u100"}[u].join([f"M{M}_",""])

def load(M,u):
    d=os.path.join(P,f"M{M}_"+{0.05:"u005",0.2:"u020",0.5:"u050",1.0:"u100"}[u])
    out=[]
    for f in sorted(glob.glob(os.path.join(d,"red_*.csv"))):
        e=pd.read_csv(f,low_memory=False)
        se=e["SegEtas"].str.split(";",expand=True).astype(float).to_numpy()
        sc=e["SegCounts"].str.split(";",expand=True).astype(float).to_numpy()
        L=sc*math.pi*RD*RD/(se*H)                      # (n,2) lengths; col 0 = LEFT, col 1 = RIGHT
        out.append(dict(t=e["Time"].to_numpy(float),
                        KEl=e["KE_gas_left"].to_numpy(float), KEr=e["KE_gas_right"].to_numpy(float),
                        v=e["W0_v"].to_numpy(float), W=e["PistonWork"].to_numpy(float),
                        Ll=L[:,0], Lr=L[:,1]))
    return out

res={}
print(f"eta_0={ETA0:.8f}  c_s={CS:.6f}  Z={ZAC:.4f}  transit={TRANSIT:.2f}  round trip={ROUND:.2f}")
print(f"W_qs targets: KR {W_KR} slope {SLOPE_KR:.3f} | ideal {W_ID} slope {SLOPE_ID:.3f}\n")
ledger_max=0.0; ledger_where=None
for M,u in CELLS:
    S=load(M,u)
    if not S: print(f"  {M},{u}: no data"); continue
    tp=D0/u; gap=0.25/u; w0,w1 = gap+tp+TRANSIT, gap+TRANSIT+ROUND   # return reaches the divider at gap+66.45
    n=min(len(s["t"]) for s in S); t=S[0]["t"][:n]
    F,Fk,NAIVE,KICK,DD,WIN,VS=[],[],[],[],[],[],[]
    for s in S:
        dKl=s["KEl"][:n]-s["KEl"][0]; dKr=s["KEr"][:n]-s["KEr"][0]
        KEd=0.5*M*s["v"][:n]**2; KEd-=KEd[0]
        X1=dKr-Eqs(s["Lr"][:n]); X2=dKl-Eqs(s["Ll"][:n])   # gas 1 = pushed = RIGHT
        W=s["W"][:n]
        r=np.abs(W-(dKl+dKr+KEd)); 
        if r.max()>ledger_max: ledger_max, ledger_where = r.max(), (M,u)
        den=X1+X2+KEd
        if w1>w0:
            k=(t>=w0)&(t<=w1)
            if k.sum()>=3:
                F.append(np.mean(np.where(np.abs(den[k])>1e-9, X2[k]/den[k], np.nan)))
                NAIVE.append(np.mean(dKl[k]/max(W[k].mean(),1e-9)))
            i=int(np.searchsorted(t,w1)); i=min(i,n-1)
            if abs(den[i])>1e-9: Fk.append((X2[i]+0.5*KEd[i])/den[i])
        pass
        pl=t>t[-1]*0.6; DD.append(((s["KEr"][:n]-s["KEl"][:n])/NS)[pl].mean())
        WIN.append(W[-1]); VS.append(s["v"][:n])
    # FIRST KICK: the COHERENT velocity, i.e. the seed MEAN, not max|v| per seed (which at M = 10 is
    # dominated by thermal motion: rms v_th = sqrt(kT/M) = 0.316 there against a predicted kick of 0.87).
    vm=np.mean(np.array(VS),axis=0); e=(t>=gap)&(t<=gap+TRANSIT+ROUND)
    kick=float(vm[e][np.argmax(np.abs(vm[e]))]) if e.sum() else float("nan")
    kick_pred=2*ZAC*D0/M
    sem=lambda a: (np.std(a,ddof=1)/math.sqrt(len(a))) if len(a)>1 else float("nan")
    x=M/(2*ZAC*tp)
    res[(M,u)]=dict(M=M,u=u,tau_push=tp,gap=gap,x=x,T=float(Tspec(x)),n=len(S),
        f=float(np.nanmean(F)) if F else float("nan"), fe=float(sem(np.array(F))) if len(F)>1 else float("nan"),
        fk=float(np.nanmean(Fk)) if Fk else float("nan"), fke=float(sem(np.array(Fk))) if len(Fk)>1 else float("nan"),
        naive=float(np.nanmean(NAIVE)) if NAIVE else float("nan"),
        kick=kick, kick_pred=float(kick_pred),
        W=float(np.mean(WIN)), We=float(sem(np.array(WIN))),
        D=float(np.mean(DD))*NS/float(np.mean(WIN)), De=float(sem(np.array(DD)))*NS/float(np.mean(WIN)),
        win=(w0,w1), has_window=bool(w1>w0))
    r=res[(M,u)]
    print(f"  M={M:3d} u={u:4.2f}  tau_push={tp:6.3f} gap={gap:5.2f}  x={x:7.4f}  T(x)={r['T']:.4f}  "
          f"window {'[%.1f, %.1f]'%(w0,w1) if w1>w0 else 'EMPTY':>16}  W_in={r['W']:7.3f}  kick dV={r['kick']:+.4f} (pred {r['kick_pred']:+.4f})")
json.dump({f"{k[0]}_{k[1]}":v for k,v in res.items()},open(os.path.join(HERE,"261008_level4b_results.json"),"w"),indent=1)
print(f"\n### LEDGER: max |W_in - (dKE_1 + dKE_2 + KE_div)| over ALL samples, ALL seeds, ALL 12 cells")
print(f"    = **{ledger_max:.3e} kT**, worst cell {ledger_where}   (float32 rounding on quantities of order 10)")

print("\n### 1. W_in(u) LEDGER: intercept and slope per M_d, fitted on u = 0.05, 0.2, 0.5")
print("| M_d | W_in(0.05) | W_in(0.2) | W_in(0.5) | intercept | slope | vs KR 7.0715/8.747 | vs ideal 5.5556/7.071 |")
print("|---|---|---|---|---|---|---|---|")
for M in (10,50,200):
    us=[0.05,0.2,0.5]; ws=[res[(M,u)]["W"] for u in us]; we=[res[(M,u)]["We"] for u in us]
    A=np.vstack([np.ones(3),us]).T
    C=np.linalg.inv(A.T@np.diag(1/np.array(we)**2)@A)
    b=C@(A.T@np.diag(1/np.array(we)**2)@np.array(ws))
    ei=math.sqrt(C[0,0]); es=math.sqrt(C[1,1])
    print(f"| {M} | {ws[0]:.3f} ± {we[0]:.3f} | {ws[1]:.3f} ± {we[1]:.3f} | {ws[2]:.3f} ± {we[2]:.3f} | "
          f"**{b[0]:.3f} ± {ei:.3f}** | **{b[1]:.2f} ± {es:.2f}** | "
          f"{abs(b[0]-W_KR)/ei:.1f}σ / {abs(b[1]-SLOPE_KR)/es:.1f}σ | {abs(b[0]-W_ID)/ei:.1f}σ / {abs(b[1]-SLOPE_ID)/es:.1f}σ |")

print("\n### 2. FIRST KICK and the abort gate (coherent, seed-mean velocity)")
print("| M_d | u | measured first-kick dV | predicted 2Zd/M | gate dV > 0.04 |")
print("|---|---|---|---|---|")
for M,u in CELLS:
    r=res[(M,u)]
    print(f"| {M} | {u} | **{r['kick']:+.4f}** | {r['kick_pred']:+.4f} | {'PASS' if abs(r['kick'])>0.04 else '**FAIL**'} |")

print("\n### 3. f_first against T(x)")
print("| M_d | u | x | T(x) | **f_first (window)** | f at 66.45 (+½KE_div) | naive ΔKE₂/W_in | in fit? |")
print("|---|---|---|---|---|---|---|---|")
fit=[]
for M,u in CELLS:
    r=res[(M,u)]
    inf = (u in (0.2,0.5))
    if not r["has_window"]:
        print(f"| {M} | {u} | {r['x']:.4f} | {r['T']:.4f} | — no first pass — | — | — | no |"); continue
    if inf and np.isfinite(r["f"]) and np.isfinite(r["fe"]) and r["fe"]>0: fit.append((r["x"],r["f"],r["fe"]))
    print(f"| {M} | {u} | {r['x']:.4f} | {r['T']:.4f} | **{r['f']:.4f} ± {r['fe']:.4f}** | {r['fk']:.4f} ± {r['fke']:.4f} | "
          f"{r['naive']:.4f} | {'**yes**' if inf else 'plotted only'} |")
if len(fit)>=3:
    xs=np.array([f[0] for f in fit]); ys=np.array([f[1] for f in fit]); ss=np.array([f[2] for f in fit])
    g=lambda xx,Zfit: Tspec(xx*ZAC/Zfit)
    try:
        pf,cf=curve_fit(g,xs,ys,p0=[ZAC],sigma=ss,absolute_sigma=True,maxfev=40000)
        chi2=float((((ys-g(xs,pf[0]))/ss)**2).sum()/max(1,len(xs)-1))
        print(f"\n  six-cell fit, Z free: **Z = {pf[0]:.4f} ± {math.sqrt(cf[0][0]):.4f}** against {ZAC:.4f} "
              f"({100*(pf[0]/ZAC-1):+.1f} %), **chi2_red = {chi2:.2f}**")
        okZ = abs(pf[0]/ZAC-1)<0.25; okc = chi2<3
        print(f"  gate Z within 25 %: {'PASS' if okZ else 'FAIL'}   gate chi2_red < 3: {'PASS' if okc else 'FAIL'}")
        print(f"  ==> **{'COLLAPSE' if (okZ and okc) else 'NOT A COLLAPSE'}**")
    except Exception as ex: print("  fit failed:",ex)
r1,r2=res.get((10,0.2)),res.get((200,1.0))
if r1 and r2:
    print(f"\n  LIMITS: (10, 0.2) f_first = {r1['f']:.4f} vs > 0.6  -> {'PASS' if r1['f']>0.6 else '**FAIL**'}")
    print(f"          (200, 1.0) f_first = {r2['f']:.4f} vs < 0.1  -> {'PASS' if r2['f']<0.1 else '**FAIL**'}")

print("\n### 4. D(u, M_d) at the plateau")
print("| M_d | u = 0.05 | 0.2 | 0.5 | 1.0 |")
print("|---|---|---|---|---|")
for M in (10,50,200):
    print(f"| {M} | " + " | ".join(f"{res[(M,u)]['D']:+.4f} ± {res[(M,u)]['De']:.4f}" for u in (0.05,0.2,0.5,1.0)) + " |")
