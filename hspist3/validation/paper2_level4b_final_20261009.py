#!/usr/bin/env python3
"""##CHRIS 2026-10-09: Level 4b, FINAL pass. Applies 261008 sections 1A-1D as committed (9488973,
c6191b5, f510014). One run, one results section, no iteration. Every number in the printed tables
is produced by this script."""
import glob,json,math,os,sys
import numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO=os.path.dirname(os.path.dirname(HERE))
P=os.path.join(REPO,"hspist3","experiments_energy_transfer","level4b_transmission_20261008b")
OUT=os.path.join(REPO,"0000_PLAN_OVERALL","paper2_energytransfer","experiments","final")
BLUE,RED,BLACK,GREY="#1f4e9c","#c0392b","#000000","#7f7f7f"

NS,H,RD,LC,D0=50,10.0,0.5,38.75,3.875
Zf=lambda e: float(sos.Z_kolafa_rottner_2006(np.array([e]))[0])
dZf=lambda e: float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
ETA0=NS*math.pi*RD*RD/(LC*H); CS=math.sqrt(Zf(ETA0)+ETA0*dZf(ETA0)+Zf(ETA0)**2)
CSID=math.sqrt(2.0); ZAC=NS*CS/LC; ZID=NS*CSID/LC
TR,RT=LC/CS,2*LC/CS
W_KR,W_ID=7.0715,5.5556; S_KR,S_ID=ZAC*D0,ZID*D0
_Lg=np.linspace(20.,55.,1400)
_Tg=np.array([math.exp(quad(lambda e: Zf(e)/e, ETA0, NS*math.pi*RD*RD/(L*H))[0]) for L in _Lg])
Eqs=lambda L: NS*(np.interp(L,_Lg,_Tg)-1.0)
Tsp=lambda x: 1.0-x*(1.0-np.exp(-1.0/x))
ANS=lambda x: 1.0/(1.0+math.pi**2*x*x)
UD={0.05:"u005",0.2:"u020",0.5:"u050",1.0:"u100"}
CELLS=[(M,u) for M in (10,50,200) for u in (0.05,0.2,0.5,1.0)]
def sem(a): a=np.asarray(a,float); return a.std(ddof=1)/math.sqrt(len(a)) if len(a)>1 else float("nan")
def jack(num,den):
    n=len(num); j=np.array([np.delete(num,i).mean()/np.delete(den,i).mean() for i in range(n)])
    return num.mean()/den.mean(), math.sqrt((n-1)/n*((j-j.mean())**2).sum())

R={}; LEDG=0.0; LWH=None
for M,u in CELLS:
    fs=sorted(glob.glob(os.path.join(P,f"M{M}_{UD[u]}","red_*.csv")))
    if not fs: continue
    tp=D0/u; gap=0.25/u; w0,w1=gap+tp+TR, gap+TR+RT
    NUM,DEN,RAT,NAI,WIN,DD,FK,VS=[],[],[],[],[],[],[],[]
    tref=None
    for f in fs:
        d=pd.read_csv(f,low_memory=False); t=d["Time"].to_numpy(float); tref=t
        se=d["SegEtas"].str.split(";",expand=True).astype(float).to_numpy()
        sc=d["SegCounts"].str.split(";",expand=True).astype(float).to_numpy()
        L=sc*math.pi*RD*RD/(se*H)
        dKl=d["KE_gas_left"].to_numpy(float)-d["KE_gas_left"].to_numpy(float)[0]
        dKr=d["KE_gas_right"].to_numpy(float)-d["KE_gas_right"].to_numpy(float)[0]
        v=d["W0_v"].to_numpy(float); KEd=0.5*M*v*v; KEd-=KEd[0]
        W=d["PistonWork"].to_numpy(float)
        r=np.abs(W-(dKl+dKr+KEd))
        if r.max()>LEDG: LEDG,LWH=r.max(),(M,u)
        X1=dKr-Eqs(L[:,1]); X2=dKl-Eqs(L[:,0]); den=X1+X2+KEd
        if w1>w0:
            k=(t>=w0)&(t<=w1)
            NUM.append(X2[k].mean()); DEN.append(den[k].mean()); RAT.append((X2[k]/den[k]).mean())
            NAI.append((dKl[k]/max(W[k].mean(),1e-9)).mean())
            i=min(int(np.searchsorted(t,w1)),len(t)-1); FK.append((X2[i]+0.5*KEd[i])/den[i])
        WIN.append(W[-1]); pl=t>t[-1]*0.6
        DD.append(((d["KE_gas_right"].to_numpy(float)-d["KE_gas_left"].to_numpy(float))/NS)[pl].mean())
        VS.append(v)
    n=min(len(a) for a in VS); vm=np.mean([a[:n] for a in VS],axis=0); t=tref[:n]
    e=(t>=gap)&(t<=gap+TR+RT); kick=float(vm[e][np.argmax(np.abs(vm[e]))]) if e.sum() else float("nan")
    kj=[]
    for i in range(len(VS)):
        q=np.mean([VS[j][:n] for j in range(len(VS)) if j!=i],axis=0)
        kj.append(q[e][np.argmax(np.abs(q[e]))])
    kj=np.array(kj); ke=math.sqrt((len(kj)-1)/len(kj)*((kj-kj.mean())**2).sum())
    x=M/(2*ZAC*tp)
    rec=dict(M=M,u=u,tp=tp,x=x,T=float(Tsp(x)),ans=float(ANS(x)),n=len(fs),
             W=float(np.mean(WIN)),We=float(sem(WIN)),
             D=float(np.mean(DD))*NS/float(np.mean(WIN)),De=float(sem(DD))*NS/float(np.mean(WIN)),
             kick=kick,kick_e=ke,kick_KR=2*ZAC*D0/M,kick_ID=2*ZID*D0/M,has=bool(w1>w0),win=(w0,w1))
    if w1>w0:
        f_ra,e_ra=jack(np.array(NUM),np.array(DEN))
        rec.update(f=float(f_ra),fe=float(e_ra),fmr=float(np.mean(RAT)),fmre=float(sem(RAT)),
                   fk=float(np.mean(FK)),fke=float(sem(FK)),naive=float(np.mean(NAI)),
                   X2=float(np.mean(NUM)),X2e=float(sem(NUM)),den=float(np.mean(DEN)),dene=float(sem(DEN)),
                   Epulse=float(ZAC*u*D0))
    R[(M,u)]=rec

print("### (a) LEDGER")
print(f"max |W_in - (dKE_1 + dKE_2 + KE_div)| over ALL samples, ALL seeds, ALL 12 cells = **{LEDG:.3e} kT**  (worst cell M={LWH[0]}, u={LWH[1]})")
print("\n### (b) W_in(u)")
print("| M_d | W_in(0.05) | W_in(0.2) | W_in(0.5) | intercept (3pt) | slope (3pt) | intercept (0.05,0.2 only) |")
print("|---|---|---|---|---|---|---|")
WFIT={}
for M in (10,50,200):
    us=np.array([0.05,0.2,0.5]); w=np.array([R[(M,u)]["W"] for u in us]); we=np.array([R[(M,u)]["We"] for u in us])
    A=np.vstack([np.ones(3),us]).T; Wt=np.diag(1/we**2); C=np.linalg.inv(A.T@Wt@A); b=C@(A.T@Wt@w)
    sl2=(w[1]-w[0])/(us[1]-us[0]); ic2=w[0]-sl2*us[0]
    WFIT[M]=(b[0],math.sqrt(C[0,0]),b[1],math.sqrt(C[1,1]),ic2,sl2)
    print(f"| {M} | {w[0]:.3f} ± {we[0]:.3f} | {w[1]:.3f} ± {we[1]:.3f} | {w[2]:.3f} ± {we[2]:.3f} | "
          f"**{b[0]:.3f} ± {math.sqrt(C[0,0]):.3f}** | {b[1]:.2f} ± {math.sqrt(C[1,1]):.2f} | **{ic2:.3f}** (slope {sl2:.2f}) |")
print(f"targets: KR intercept {W_KR} slope {S_KR:.3f} | ideal intercept {W_ID} slope {S_ID:.3f}")
print("\n### (c) FIRST KICK at M = 200 (coherent, seed-mean velocity; jackknife error)")
print("| u | measured dV | 2Zd/M (KR) | 2Z_id d/M (ideal) | gate |dV| > 0.04 |")
print("|---|---|---|---|---|")
for u in (0.05,0.2,0.5,1.0):
    r=R[(200,u)]
    print(f"| {u} | **{r['kick']:+.4f} ± {r['kick_e']:.4f}** | {-r['kick_KR']:+.4f} | {-r['kick_ID']:+.4f} | "
          f"{'PASS' if abs(r['kick'])>0.04 else '**FAIL**'} |")
print("\n### (d) f_first")
print("| M | u | x | T(x) | ansatz | **f (ratio of means)** | f (mean of ratios) | f at 66.45 +½KE_div | naive ΔKE₂/W_in | <X₂> | <denom> | E_pulse | in fit |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
fit=[]
for M,u in CELLS:
    r=R[(M,u)]
    if not r["has"]:
        print(f"| {M} | {u} | {r['x']:.4f} | {r['T']:.4f} | {r['ans']:.4f} | — no first pass — | — | — | — | — | — | — | no |"); continue
    inf=u in (0.2,0.5)
    if inf: fit.append((r["x"],r["f"],r["fe"]))
    print(f"| {M} | {u} | {r['x']:.4f} | {r['T']:.4f} | {r['ans']:.4f} | **{r['f']:+.4f} ± {r['fe']:.4f}** | {r['fmr']:+.4f} | "
          f"{r['fk']:+.4f} ± {r['fke']:.4f} | {r['naive']:.4f} | {r['X2']:+.4f} ± {r['X2e']:.4f} | {r['den']:+.4f} ± {r['dene']:.4f} | {r['Epulse']:.3f} | {'**yes**' if inf else 'plot'} |")
xs=np.array([f[0] for f in fit]); ys=np.array([f[1] for f in fit]); ss=np.array([f[2] for f in fit])
g=lambda xx,Z: Tsp(xx*ZAC/Z)
try:
    pf,cf=curve_fit(g,xs,ys,p0=[ZAC],sigma=ss,absolute_sigma=True,maxfev=60000)
    Zf_,Ze=float(pf[0]),float(math.sqrt(cf[0][0])); chi2=float((((ys-g(xs,Zf_))/ss)**2).sum()/max(1,len(xs)-1))
except Exception as ex:
    Zf_,Ze,chi2=float("nan"),float("nan"),float("nan"); print("fit failed:",ex)
okZ=abs(Zf_/ZAC-1)<0.25 if Zf_==Zf_ else False; okc=chi2<3 if chi2==chi2 else False
print(f"\nsix-cell pre-registered fit (u = 0.2, 0.5), Z free: **Z = {Zf_:.4f} ± {Ze:.4f}** vs {ZAC:.4f}; **chi2_red = {chi2:.2f}**")
print(f"  Z within 25 %: {'PASS' if okZ else 'FAIL'};  chi2_red < 3: {'PASS' if okc else 'FAIL'}  ==> **{'COLLAPSE' if okZ and okc else 'MODEL TEST FAILED'}**")
print(f"LIMITS: (10,0.2) f = {R[(10,0.2)]['f']:+.4f} vs > 0.6 -> {'PASS' if R[(10,0.2)]['f']>0.6 else '**FAIL**'};  "
      f"(200,1.0) f = {R[(200,1.0)]['f']:+.4f} vs < 0.1 -> {'PASS' if R[(200,1.0)]['f']<0.1 else '**FAIL**'}")
print("\n### (e) D(u, M_d)")
print("| M_d | u = 0.05 | 0.2 | 0.5 | 1.0 |")
print("|---|---|---|---|---|")
for M in (10,50,200):
    print(f"| {M} | " + " | ".join(f"{R[(M,u)]['D']:+.4f} ± {R[(M,u)]['De']:.4f}" for u in (0.05,0.2,0.5,1.0)) + " |")
r=R[(200,1.0)]; print(f"the 32-seed cell: **D(200, 1.0) = {r['D']:+.4f} ± {r['De']:.4f}** on {r['n']} seeds")

fig,ax=plt.subplots(1,2,figsize=(11.5,4.6))
xg=np.logspace(-1.6,1.3,400)
p=ax[0]
p.plot(xg,Tsp(xg),"-",color=BLACK,lw=1.8,label=r"$T(x)=1-x(1-e^{-1/x})$  (prediction)")
p.plot(xg,[ANS(v) for v in xg],"--",color=GREY,lw=1.3,label=r"first-written ansatz $1/(1+\pi^2x^2)$")
for us,mk,fc,lb in (((0.2,0.5),"o",BLUE,"fitted (u = 0.2, 0.5)"),((1.0,),"s","white","u = 1.0, plotted only")):
    X=[R[(M,u)]["x"] for M,u in CELLS if u in us and R[(M,u)]["has"]]
    Y=[R[(M,u)]["f"] for M,u in CELLS if u in us and R[(M,u)]["has"]]
    E=[R[(M,u)]["fe"] for M,u in CELLS if u in us and R[(M,u)]["has"]]
    p.errorbar(X,Y,yerr=E,fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=lb)
Xn=[R[(M,u)]["x"] for M,u in CELLS if R[(M,u)]["has"]]; Yn=[R[(M,u)]["naive"] for M,u in CELLS if R[(M,u)]["has"]]
p.plot(Xn,Yn,"^",color=RED,ms=6,ls="none",label=r"naive $\Delta E_2/W_{\rm in}$ (the trap)")
p.set_xscale("log"); p.set_xlabel(r"$x=M_d/(2Z\tau_{\rm push})$"); p.set_ylabel(r"$f_{\rm first}$")
p.set_ylim(-0.2,1.05); p.axhline(0,color=GREY,lw=.7)
p.set_title("(a) first-pass transmission vs the spectral average",fontsize=10); p.legend(fontsize=7.5); p.grid(alpha=.25,which="both")
p=ax[1]
for M,mk,fc in ((10,"o",BLUE),(50,"s","white"),(200,"^",BLUE)):
    us=[0.05,0.2,0.5]; p.errorbar(us,[R[(M,u)]["W"] for u in us],yerr=[R[(M,u)]["We"] for u in us],
        fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=f"$M_d={M}$")
uu=np.linspace(0,0.55,50)
p.plot(uu,W_KR+S_KR*uu,"-",color=RED,lw=1.6,label=r"KR: $W_{qs}+Zdu$")
p.plot(uu,W_ID+S_ID*uu,"--",color=GREY,lw=1.3,label="ideal gas")
p.set_xlabel(r"$u$  [$\sigma$/$\sigma$-time]"); p.set_ylabel(r"$W_{\rm in}$  [$kT$]")
p.set_title("(b) the piston-work ledger",fontsize=10); p.legend(fontsize=7.5); p.grid(alpha=.25)
fig.suptitle(r"Level 4b: first-pass transmission and the work ledger ($\eta=0.10134170$, $Z=2.2572$)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.94])
for ext in ("png","pdf"):
    fp=os.path.join(OUT,f"261008_p2_level4b_transmission.{ext}"); fig.savefig(fp,dpi=180)
print(f"\n### (f) figure: {os.path.join(OUT,'261008_p2_level4b_transmission.png')}")
json.dump({f"{k[0]}_{k[1]}":v for k,v in R.items()} | {"_meta":dict(ledger_max=float(LEDG),Z=ZAC,Zfit=Zf_,Zfit_err=Ze,chi2=chi2,
           Wfit={str(k):list(map(float,v)) for k,v in WFIT.items()})},
          open(os.path.join(HERE,"261008_level4b_results.json"),"w"),indent=1)
