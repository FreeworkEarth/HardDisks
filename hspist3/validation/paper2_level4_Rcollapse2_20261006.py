#!/usr/bin/env python3
"""##CHRIS 2026-10-06: R-collapse, PASS 2. Analysis only; no core physics touched.

THE PREDICTIONS BELOW WERE COMMITTED BEFORE ANY RECORD OF THIS CAMPAIGN WAS FITTED
(git 3af501c, 0000_PLAN_OVERALL/ALL_MARKDOWNS/261006_paper2_level4_Rcollapse2.md, sections 4 and 5).
Nothing in this file may be edited to change them; if the verdict is uncomfortable, the verdict is
the result.

WHAT PASS 1 GOT WRONG. It phrased "the variable is R" as g = tau_T/M carrying over unchanged.
tau_GP is proportional to M x L_c, so doubling N_s at fixed eta doubles L_c and doubles tau_GP/M
(50.349 M -> 100.698 M). A hypothesis phrased in g builds the L_c scaling out of the thing being
tested. The admissible invariant is

    f = tau_T / tau_GP ,  tau_GP = (4/sqrt(2 pi)) M L_c/sqrt(m kT)/(1 + eta Z'/Z)

with each box using its own L_c. "The variable is R" means f = f(R), R = N_s m/M.

PRE-REGISTERED PREDICTIONS (interpolations quoted in the report so they cannot be re-chosen):
    A -- f = f(R):        f(4) = 0.98 +- 0.05   f(2) = 1.18 +- 0.06   f(1) = 1.52 +- 0.04
    B -- f = f_lad(M)/2:  f(25) = 0.59 +- 0.02  f(50) = 0.76 +- 0.02  f(100) = 1.22 +- 0.09
VERDICT RULE:
    gates  -- record/tau >= 60 AND calibration slope >= 0.6, at ALL THREE cells;
    A wins -- within 2 sigma of A at all three AND B excluded > 3 sigma at >= 2 of them;
    B wins -- symmetrically; otherwise NOT RESOLVED, failing diagnostic reported.
    sigma  = hypot(measurement error, prediction error). Calibration grid centred on the MEASURED
             tau (the 2026-10-02 clamping trap), never on a prediction.

TWO DELIBERATE DIFFERENCES FROM THE PASS-1 SCRIPT, both stated in the report:
  1. dt is read from the recorded Time column, not from the hardcoded SPS = 612000/10195.000532.
     That constant is 0.049 % high (the runs record dt_sigma = 1/60 exactly), so pass 1's tau was
     0.049 % low -- far inside its 3-9 % error, but there is no reason to carry an avoidable bias.
  2. The nested-record diagnostic: M = 50 and 100 reuse pass 1's seeds AND pass 1's --trace-every,
     so each pass-2 trace strictly contains the pass-1 trace. Truncating pass 2 back to pass 1's
     length and re-fitting measures the record-length bias with ZERO seed noise. Record length is
     the systematic that failed pass 1, so it is worth measuring rather than assuming.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
from scipy.signal import lfilter
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO=os.path.dirname(os.path.dirname(HERE))
ET=os.path.join(REPO,"hspist3","experiments_energy_transfer")

NS, LC, BURN, GAMMA = 100, 77.5, 2000.0, 4.118
LEFF = LC-1.0                                   # L_c - 2r, as every Level 4 mode calculation
ETA  = NS*math.pi*0.25/(LC*10.0)                # = 200 pi r^2/(2 L_c h) = 0.10134170
Zv   = float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0])
dZv  = float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
CS   = math.sqrt(Zv+ETA*dZv+Zv*Zv)
GP   = (4/math.sqrt(2*math.pi))*LC/(1+ETA*dZv/Zv)           # tau_GP per unit mass, THIS box

RUN2 = os.path.join(ET,"level4_Rcollapse2_20261006")
RUN1 = os.path.join(ET,"level4_Rcollapse_20261005")
CELLS = {25:os.path.join(RUN2,"Md25"), 50:os.path.join(RUN2,"Md50"), 100:os.path.join(RUN2,"Md100")}
PASS1 = {50:(os.path.join(RUN1,"Md50"),7253.,329426.5), 100:(os.path.join(RUN1,"Md100"),15773.,803838.5)}
PRED_A = {25:(0.98,0.05), 50:(1.18,0.06), 100:(1.52,0.04)}
PRED_B = {25:(0.59,0.02), 50:(0.76,0.02), 100:(1.22,0.09)}
MASSES = (25,50,100)

def kroot(a): return brentq(lambda k: math.cos(k)/math.sin(k)-a*k, 1e-9, math.pi-1e-9)

def acf(x,nl):
    x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m)
    c=np.fft.irfft(f*np.conj(f))[:nl].real; return c/c[0]

def load(d):
    """dT/N per seed, plus dt measured from the recorded Time column (not a hardcoded constant)."""
    D=[]; dts=[]
    for f in sorted(glob.glob(os.path.join(d,"red_*.csv"))):
        e=pd.read_csv(f,usecols=["Time","KE_gas_left","KE_gas_right"])
        t=e["Time"].to_numpy(float)
        dts.append((t[-1]-t[0])/(len(t)-1))
        D.append((e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))/NS)
    if not D: raise SystemExit(f"no reduced traces in {d}")
    dt=float(np.median(dts))
    if max(dts)-min(dts) > 1e-6*dt:
        print(f"  NOTE {os.path.basename(d)}: dt varies by {1e6*(max(dts)-min(dts))/dt:.1f} ppm across seeds")
    n=min(len(a) for a in D); lo=int(BURN/dt)
    return [a[lo:n] for a in D], dt

def mk_fit(nl,dt,om):
    lag=np.arange(nl)*dt
    def f(series,p0):
        c=np.mean([acf(s,nl) for s in series],axis=0)
        g=lambda t,A,tT,B,tr: A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(om*t)
        p,_=curve_fit(g,lag,c,p0=p0,bounds=([0,20,0,20],[1.5,1e6,1.5,5e4]),maxfev=120000)
        return p
    return f

def fit_block(series,per,dt):
    B=max(2,int(round(per/dt)))
    bl=[s[:len(s)//B*B].reshape(-1,B).mean(1) for s in series]
    bdt=B*dt; k=min(len(bl[0])//2,16)
    c=np.mean([acf(b,k) for b in bl],axis=0); lg=np.arange(k)*bdt; m=c>0.05
    if m.sum()<3: return float("nan")
    s=np.polyfit(lg[m],np.log(c[m]),1)[0]
    return -1/s if s<0 else float("nan")

def fit_S0(series,dt,nu,tau0):
    """Amplitude FREE: tau comes from the corner frequency, not from the assumed variance."""
    sp=[]
    for s in series:
        x=s-s.mean(); w=np.hanning(len(x)); x=x*w
        sp.append(np.abs(np.fft.rfft(x))**2*2.0*dt/((w**2).sum()))
    S=np.mean(sp,axis=0); fr=np.fft.rfftfreq(len(series[0]),dt); m=(fr>0)&(fr<0.25*nu)
    g=lambda f,S0,tau: S0/(1.0+(2*math.pi*f*tau)**2)
    p,cov=curve_fit(g,fr[m],S[m],p0=[S[m][0],tau0],bounds=([1e-12,1.0],[1e8,1e8]),maxfev=80000)
    return float(p[1]),float(math.sqrt(cov[1][1]))

def synth(tau,nraw,nseed,om,dt,TR0,A0,B0,rng):
    a=math.exp(-dt/tau); rho=math.exp(-dt/TR0)
    c1=2*rho*math.cos(om*dt); c2=-rho*rho
    E=rng.normal(size=(nseed,nraw))*math.sqrt(1-a*a); W=rng.normal(size=(nseed,nraw))
    S1=lfilter([1.0],[1.0,-a],E,axis=1); V=lfilter([1.0],[1.0,-c1,-c2],W,axis=1)
    sd=V[:,200:].std(axis=1,keepdims=True); V/=np.where(sd>0,sd,1.0)
    return list(math.sqrt(A0)*S1+math.sqrt(B0)*V)

def measure(D,dt,M,seed_off,label=""):
    """Three estimators + bias calibration centred on the measured tau. Returns a dict."""
    a=M/(2.0*NS); K=kroot(a); nu=CS*K/(2*math.pi*LEFF); om=2*math.pi*nu; per=1/nu
    nraw=len(D[0]); rec=nraw*dt
    tau=GP*M; TR0=0.4*per*5; A0,B0=0.50,0.47
    for _ in range(5):                                   # self-consistent 8-tau fit window
        nl=max(64,int(min(8.0*tau,rec/3)/dt)); f=mk_fit(nl,dt,om)
        try: p=f(D,[A0,tau,B0,TR0])
        except Exception: break
        if abs(p[1]-tau)<0.02*tau: tau=p[1]; break
        tau=p[1]
    nl=max(64,int(min(8.0*tau,rec/3)/dt)); f=mk_fit(nl,dt,om)
    p=f(D,[A0,tau,B0,TR0]); raw=p[1]; TR0=p[3]; A0,B0=p[0],p[2]
    rb=fit_block(D,per,dt)
    rng=np.random.default_rng(seed_off+M); tr,rc,sd=[],[],[]
    for gm in (0.5,0.7,0.9,1.0,1.2,1.5,2.0):             # grid centred on the MEASURED tau
        true=gm*raw/0.75; got=[]
        for _ in range(6):
            try:
                v=f(synth(true,nraw,len(D),om,dt,TR0,A0,B0,rng),[A0,true,B0,TR0])[1]
                if np.isfinite(v): got.append(v)
            except Exception: pass
        if len(got)>2: tr.append(true); rc.append(np.mean(got)); sd.append(np.std(got,ddof=1))
    tr,rc,sd=np.array(tr),np.array(rc),np.array(sd)
    slope=float(np.polyfit(tr,rc,1)[0])
    clamped = not (rc.min() <= raw <= rc.max())
    vm=float(np.interp(raw,rc,tr)); em=float(np.interp(raw,rc,sd))/max(slope,1e-9)
    vb=float(np.interp(rb,rc,tr)) if rc.min()<=rb<=rc.max() else float("nan")
    s0,_=fit_S0(D,dt,nu,raw)
    vs0=float(np.interp(s0,rc,tr)) if rc.min()<=s0<=rc.max() else s0
    vals=[v for v in (vm,vb,vs0) if np.isfinite(v)]
    return dict(M=M,R=NS/M,alpha=a,K=K,period=per,rec=rec,dt=dt,nseed=len(D),
                tau=vm,err=em,block=vb,S0=vs0,spread=max(vals)/min(vals),
                slope=slope,clamped=bool(clamped),Ltau=rec/vm,
                g=vm/M,ge=em/M,f=vm/(GP*M),fe=em/(GP*M),
                P_recoll=(M/GAMMA)*CS/(2*LC),tau_r=float(TR0),A=float(A0),B=float(B0),label=label)

def main():
    print(f"N_s={NS}  L_c={LC}  eta={ETA:.8f}  Z={Zv:.6f}  etaZ'={ETA*dZv:.6f}  c_s={CS:.6f}")
    print(f"L_eff={LEFF}  rho_0={NS/LC:.5f}  tau_GP = {GP:.3f} M   (ladder box: {GP/2:.3f} M)\n")

    out={}
    for M in MASSES:
        D,dt=load(CELLS[M]); out[M]=measure(D,dt,M,12000,"pass2")
        r=out[M]
        print(f"  M={M:3d} R={r['R']:.0f}  {r['nseed']} seeds  record {r['rec']:.0f}  dt {r['dt']:.4f}  "
              f"L/tau {r['Ltau']:.1f}  slope {r['slope']:.3f}  tau {r['tau']:.0f} +- {r['err']:.0f}"
              + ("   *** CALIBRATION CLAMPED ***" if r['clamped'] else ""))

    print("\n### 1. Measurements\n")
    print("| M | R | L/tau | cal slope | modelled | block | S(0) | spread | tau_T | g = tau/M | **f = tau/tau_GP** |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for M in MASSES:
        r=out[M]
        print(f"| {M} | {r['R']:.0f} | **{r['Ltau']:.0f}** | {r['slope']:.3f} | {r['tau']:.0f} ± {r['err']:.0f} | "
              f"{r['block']:.0f} | {r['S0']:.0f} | {r['spread']:.2f}× | {r['tau']:.0f} ± {r['err']:.0f} | "
              f"{r['g']:.1f} ± {r['ge']:.1f} | **{r['f']:.3f} ± {r['fe']:.3f}** |")

    print("\n### 2. The pre-registered test\n")
    print("| M | R | f measured | A: f(R) | σ from A | B: f_lad(M)/2 | σ from B | P = 0.274/R |")
    print("|---|---|---|---|---|---|---|---|")
    sA={}; sB={}
    for M in MASSES:
        r=out[M]; (fa,sa)=PRED_A[M]; (fb,sb)=PRED_B[M]
        sA[M]=abs(r['f']-fa)/math.hypot(r['fe'],sa); sB[M]=abs(r['f']-fb)/math.hypot(r['fe'],sb)
        print(f"| {M} | {r['R']:.0f} | **{r['f']:.3f} ± {r['fe']:.3f}** | {fa:.2f} ± {sa:.2f} | "
              f"**{sA[M]:.1f}σ** | {fb:.2f} ± {sb:.2f} | **{sB[M]:.1f}σ** | {r['P_recoll']:.3f} |")

    c_slope=all(out[M]['slope']>=0.6 for M in MASSES)
    c_L=all(out[M]['Ltau']>=60 for M in MASSES)
    c_clamp=not any(out[M]['clamped'] for M in MASSES)
    okA=all(sA[M]<2 for M in MASSES) and sum(1 for M in MASSES if sB[M]>3)>=2
    okB=all(sB[M]<2 for M in MASSES) and sum(1 for M in MASSES if sA[M]>3)>=2
    print(f"\n  gate: L/tau >= 60 at all three     : {'PASS' if c_L else 'FAIL'}"
          f"  ({', '.join(f'{out[M]['Ltau']:.0f}' for M in MASSES)})")
    print(f"  gate: cal slope >= 0.6 at all three : {'PASS' if c_slope else 'FAIL'}"
          f"  ({', '.join(f'{out[M]['slope']:.3f}' for M in MASSES)})")
    print(f"  gate: calibration not clamped       : {'PASS' if c_clamp else 'FAIL'}")
    print(f"  A within 2σ at all three            : {okA and 'yes' or 'no'}")
    print(f"  B within 2σ at all three            : {okB and 'yes' or 'no'}")
    if not (c_slope and c_L and c_clamp):
        verdict="NOT RESOLVED -- sampling gate failed"
    elif okA and not okB: verdict="PREDICTION A SELECTED: the variable is R = N_s m/M, not M"
    elif okB and not okA: verdict="PREDICTION B SELECTED: the variable is M, not R"
    else:                 verdict="NOT RESOLVED -- neither prediction passes the 2σ/3σ rule"
    print(f"\n  ==> **{verdict}**")

    print("\n### 3. Nested-record diagnostic (zero seed noise)\n")
    print("Pass 2 reuses pass 1's seeds and cadence at M = 50 and 100, so truncating pass 2 to")
    print("pass 1's record length isolates the record-length bias from seed noise.\n")
    print("| M | pass 1 record | pass 1 tau | pass 2 truncated to pass 1 | pass 2 full | full/short |")
    print("|---|---|---|---|---|---|")
    nest={}
    for M,(d1,tau1,rec1) in PASS1.items():
        D,dt=load(CELLS[M]); n1=int(rec1/dt)
        if n1>len(D[0]): print(f"| {M} | {rec1:.0f} | {tau1:.0f} | (pass 2 shorter?!) | | |"); continue
        t=measure([a[:n1] for a in D],dt,M,13000,"trunc")
        nest[M]=dict(tau1=tau1,trunc=t['tau'],trunc_err=t['err'],full=out[M]['tau'],
                     ratio=out[M]['tau']/t['tau'])
        print(f"| {M} | {rec1:.0f} | {tau1:.0f} | **{t['tau']:.0f} ± {t['err']:.0f}** | "
              f"{out[M]['tau']:.0f} ± {out[M]['err']:.0f} | **{out[M]['tau']/t['tau']:.3f}×** |")

    json.dump(dict(meta=dict(NS=NS,LC=LC,ETA=ETA,Z=Zv,etaZp=ETA*dZv,CS=CS,GP_per_M=GP,
                             LEFF=LEFF,rho0=NS/LC,GAMMA=GAMMA,verdict=verdict,
                             gates=dict(Ltau=c_L,slope=c_slope,unclamped=c_clamp),
                             predA={str(k):v for k,v in PRED_A.items()},
                             predB={str(k):v for k,v in PRED_B.items()}),
                   cells={str(k):v for k,v in out.items()},
                   sigma_A={str(k):sA[k] for k in sA}, sigma_B={str(k):sB[k] for k in sB},
                   nested={str(k):v for k,v in nest.items()}),
              open(os.path.join(HERE,"261006_Rcollapse2_results.json"),"w"),indent=1)
    print(f"\nwrote 261006_Rcollapse2_results.json")

if __name__=="__main__":
    main()
