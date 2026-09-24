#!/usr/bin/env python3
"""##CHRIS 2026-10-05: the R-collapse test. Is the variable R = N_s m/M, or M itself?

PRE-REGISTERED (same rules as the ladder; predictions fixed before the data):
  prediction A -- the variable is R:  g(R=2) = 59.2, g(R=1) = 76.7   (interpolated from the
                  N_s = 50 ladder's g(R), which spans R = 5..0.25)
  prediction B -- the variable is M:  g(M=50) = 76.7, g(M=100) = 123.3
VERDICT RULE: whichever prediction is within 2 sigma of the calibrated tau_T/M at BOTH masses,
with the other excluded at > 3 sigma, is selected. Anything else is "not resolved".
Required throughout: adoption slope >= 0.6, L/tau >= 60, three estimators.

Geometry: N_s = 100 per side, L_c = 2 x 38.75 = 77.5 EXACTLY (eta identical to 8 digits,
0.10134170; grid-exact), box 156.0, divider centre 78.0, t = 1.0 recorded. rho_0 = 1.2903, the
same as the ladder -- which is the path Cencini's limit (p.4, Sect. II.B) is defined along.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
from scipy.signal import lfilter
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO=os.path.dirname(os.path.dirname(HERE))
ET=os.path.join(REPO,"hspist3","experiments_energy_transfer")
NS, LC, BURN = 100, 77.5, 2000.0
LEFF = LC-1.0
SPS = 612000/10195.000532
ETA = NS*math.pi*0.25/(LC*10.0)
Zv=float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0]); dZv=float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
CS=math.sqrt(Zv+ETA*dZv+Zv*Zv)
GP=(4/math.sqrt(2*math.pi))*LC/(1+ETA*dZv/Zv)
CELLS={50:(os.path.join(ET,"level4_Rcollapse_20261005","Md50"),590),
       100:(os.path.join(ET,"level4_Rcollapse_20261005","Md100"),690)}
PRED_A={50:59.2, 100:76.7}
PRED_B={50:76.7, 100:123.3}

def kroot(a): return brentq(lambda k: math.cos(k)/math.sin(k)-a*k,1e-9,math.pi-1e-9)
def acf(x,nl):
    x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m)
    c=np.fft.irfft(f*np.conj(f))[:nl].real; return c/c[0]
def load(M):
    d,every=CELLS[M]; dt=every/SPS
    D=[]
    for f in sorted(glob.glob(os.path.join(d,"red_*.csv"))):
        e=pd.read_csv(f,usecols=["KE_gas_left","KE_gas_right"])
        D.append((e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))/NS)
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

def main():
    print(f"N_s={NS}, L_c={LC}, eta={ETA:.8f}, c_s={CS:.5f}, L_eff={LEFF}, rho_0={NS/LC:.4f}")
    print(f"tau_GP = {GP:.2f} M\n")
    print("| M | R | L/tau | slope | modelled | block | S(0) | spread | g = tau/M | pred A | pred B |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    out={}
    for M in (50,100):
        D,dt=load(M); a=M/(2.0*NS); K=kroot(a); nu=CS*K/(2*math.pi*LEFF); om=2*math.pi*nu; per=1/nu
        nraw=len(D[0]); rec=nraw*dt
        tau=GP*M; TR0=0.4*per*5; A0,B0=0.50,0.47
        for _ in range(5):
            nl=max(64,int(min(8.0*tau,rec/3)/dt)); f=mk_fit(nl,dt,om)
            try: p=f(D,[A0,tau,B0,TR0])
            except Exception: break
            if abs(p[1]-tau)<0.02*tau: tau=p[1]; break
            tau=p[1]
        nl=max(64,int(min(8.0*tau,rec/3)/dt)); f=mk_fit(nl,dt,om)
        p=f(D,[A0,tau,B0,TR0]); raw=p[1]; TR0=p[3]; A0,B0=p[0],p[2]
        rb=fit_block(D,per,dt)
        rng=np.random.default_rng(11000+M); tr,rc,sd=[],[],[]
        for gm in (0.5,0.7,0.9,1.0,1.2,1.5,2.0):
            true=gm*raw/0.75; got=[]
            for _ in range(6):
                try:
                    v=f(synth(true,nraw,len(D),om,dt,TR0,A0,B0,rng),[A0,true,B0,TR0])[1]
                    if np.isfinite(v): got.append(v)
                except Exception: pass
            if len(got)>2: tr.append(true); rc.append(np.mean(got)); sd.append(np.std(got,ddof=1))
        tr,rc,sd=np.array(tr),np.array(rc),np.array(sd)
        slope=float(np.polyfit(tr,rc,1)[0])
        vm=float(np.interp(raw,rc,tr)); em=float(np.interp(raw,rc,sd))/max(slope,1e-9)
        vb=float(np.interp(rb,rc,tr)) if rc.min()<=rb<=rc.max() else float("nan")
        s0,_=fit_S0(D,dt,nu,raw)
        vs0=float(np.interp(s0,rc,tr)) if rc.min()<=s0<=rc.max() else s0
        vals=[v for v in (vm,vb,vs0) if np.isfinite(v)]
        g=vm/M; ge=em/M
        out[M]=dict(tau=vm,err=em,g=g,ge=ge,slope=slope,Ltau=rec/vm,rec=rec,
                    block=vb,S0=vs0,spread=max(vals)/min(vals))
        print(f"| {M} | {NS/M:.0f} | **{rec/vm:.0f}** | **{slope:.3f}** | {vm:.0f} ± {em:.0f} | {vb:.0f} | {vs0:.0f} | "
              f"{max(vals)/min(vals):.2f}x | **{g:.1f} ± {ge:.1f}** | {PRED_A[M]} | {PRED_B[M]} |")
    json.dump({str(k):v for k,v in out.items()},open(os.path.join(HERE,"261005_Rcollapse_results.json"),"w"),indent=1)
    print("\n### VERDICT RULE, applied as written\n")
    print("| M | g measured | pred A (R) | σ from A | pred B (M) | σ from B |")
    print("|---|---|---|---|---|---|")
    sA={}; sB={}
    for M in (50,100):
        g=out[M]["g"]; ge=out[M]["ge"]
        sA[M]=abs(g-PRED_A[M])/ge; sB[M]=abs(g-PRED_B[M])/ge
        print(f"| {M} | {g:.1f} ± {ge:.1f} | {PRED_A[M]} | **{sA[M]:.1f}σ** | {PRED_B[M]} | **{sB[M]:.1f}σ** |")
    okA=all(sA[M]<2 for M in sA) and all(sB[M]>3 for M in sB)
    okB=all(sB[M]<2 for M in sB) and all(sA[M]>3 for M in sA)
    c_slope=all(out[M]["slope"]>=0.6 for M in out); c_L=all(out[M]["Ltau"]>=60 for M in out)
    print(f"\n  slope >= 0.6 at both masses : {'PASS' if c_slope else 'FAIL'}")
    print(f"  L/tau >= 60 at both masses  : {'PASS' if c_L else 'FAIL'}")
    if not (c_slope and c_L):
        print("\n  ==> sampling criteria failed; NOT RESOLVED")
    elif okA:
        print("\n  ==> PREDICTION A SELECTED: the variable is R = N_s m/M, not M.")
    elif okB:
        print("\n  ==> PREDICTION B SELECTED: the variable is M, not R.")
    else:
        print("\n  ==> NOT RESOLVED: neither prediction satisfies within-2σ at both masses with the other excluded at >3σ")

if __name__=="__main__":
    main()
