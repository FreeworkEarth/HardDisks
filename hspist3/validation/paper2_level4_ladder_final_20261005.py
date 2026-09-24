#!/usr/bin/env python3
"""##CHRIS 2026-10-04: analysis of the ladder RERUN at 65 tau_true.

Rules fixed in 261004_paper2_level4_ladder_rerun.md section 1, written while the runs executed:
three estimators (modelled / block / S(0) amplitude-free); calibration grid CENTRED ON THE MEASURED
tau spanning 0.5-2x it; adoption slope >= 0.6; self-consistent 8-tau fit window; exponent from all
five, {10,20,50}, {50,100,200}.

VERDICT RULE, fixed before the data: b is a result only if EVERY mass has L/tau >= 60 AND
calibration slope >= 0.6 AND the light and heavy subsets agree within 2 sigma. Otherwise report
the failing diagnostic and stop.

eta_phys = 0.101342, from the RECORDED wall_thickness_sigma = 1.0.
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
from scipy.signal import lfilter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos

REPO = os.path.dirname(os.path.dirname(HERE))
ET = os.path.join(REPO, "hspist3", "experiments_energy_transfer")
NS, NTOT, LEFF, BURN = 50, 100, 37.75, 2000.0
SPS = 612000 / 10195.000532
ETA = 0.101342
Zv = float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0])
dZv = float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
CS = math.sqrt(Zv + ETA * dZv + Zv * Zv)
A_HARD = (4/math.sqrt(2*math.pi)) * 38.75 / (1 + ETA*dZv/Zv)
A_IDEAL = (4/math.sqrt(2*math.pi)) * 38.75

CELLS = {
 10:  (os.path.join(ET, "level4_topup_20261005", "Md10"), 200),
 20:  (os.path.join(ET, "level4_topup_20261005", "Md20"), 50),
 50:  (os.path.join(ET, "level4_ladder_rerun_20261004", "Md50"), 350),
 100: (os.path.join(ET, "level4_ladder_rerun_20261004", "Md100"), 450),
 200: (os.path.join(ET, "level4_ladder_rerun_20261004", "Md200"), 600),
}

def kroot(a): return brentq(lambda k: math.cos(k)/math.sin(k)-a*k, 1e-9, math.pi-1e-9)
def acf(x, nl):
    x = x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m)
    c=np.fft.irfft(f*np.conj(f))[:nl].real; return c/c[0]

def load(M):
    d, every = CELLS[M]; dt = every/SPS
    fs = sorted(glob.glob(os.path.join(d,"red_*.csv")))
    D=[]
    for f in fs:
        e=pd.read_csv(f, usecols=["KE_gas_left","KE_gas_right"])
        D.append((e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))/NS)
    n=min(len(a) for a in D); lo=int(BURN/dt)
    return [a[lo:n] for a in D], dt

def mk_fit(nl, dt, om):
    lag=np.arange(nl)*dt
    def f(series, p0):
        c=np.mean([acf(s,nl) for s in series],axis=0)
        g=lambda t,A,tT,B,tr: A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(om*t)
        p,_=curve_fit(g,lag,c,p0=p0,bounds=([0,20,0,20],[1.5,1e6,1.5,5e4]),maxfev=120000)
        return p
    return f

def fit_block(series, per, dt):
    B=max(2,int(round(per/dt)))
    bl=[s[:len(s)//B*B].reshape(-1,B).mean(1) for s in series]
    bdt=B*dt; k=min(len(bl[0])//2,16)
    c=np.mean([acf(b,k) for b in bl],axis=0); lg=np.arange(k)*bdt
    m=c>0.05
    if m.sum()<3: return float("nan")
    s=np.polyfit(lg[m],np.log(c[m]),1)[0]
    return -1/s if s<0 else float("nan")

def fit_S0(series, dt, nu, tau0):
    sp=[]
    for s in series:
        x=s-s.mean(); w=np.hanning(len(x)); x=x*w
        sp.append(np.abs(np.fft.rfft(x))**2*2.0*dt/((w**2).sum()))
    S=np.mean(sp,axis=0); fr=np.fft.rfftfreq(len(series[0]),dt)
    m=(fr>0)&(fr<0.25*nu)
    g=lambda f,S0,tau: S0/(1.0+(2*math.pi*f*tau)**2)
    p,cov=curve_fit(g,fr[m],S[m],p0=[S[m][0],tau0],bounds=([1e-12,1.0],[1e8,1e8]),maxfev=80000)
    return float(p[1]), float(math.sqrt(cov[1][1]))

def synth(tau,nraw,nseed,om,dt,TR0,A0,B0,rng):
    a=math.exp(-dt/tau); rho=math.exp(-dt/TR0)
    c1=2*rho*math.cos(om*dt); c2=-rho*rho
    E=rng.normal(size=(nseed,nraw))*math.sqrt(1-a*a); W=rng.normal(size=(nseed,nraw))
    S1=lfilter([1.0],[1.0,-a],E,axis=1); V=lfilter([1.0],[1.0,-c1,-c2],W,axis=1)
    sd=V[:,200:].std(axis=1,keepdims=True); V/=np.where(sd>0,sd,1.0)
    return list(math.sqrt(A0)*S1+math.sqrt(B0)*V)

def powerlaw(Ms,tau,err):
    Ms=np.asarray(Ms,float); tau=np.asarray(tau,float); err=np.asarray(err,float)
    w=1.0/(err/tau)**2; A=np.vstack([np.log(Ms),np.ones_like(Ms)]).T
    cov=np.linalg.inv(A.T@(A*w[:,None])); b=cov@(A.T@(w*np.log(tau)))
    return float(b[0]), float(math.sqrt(cov[0,0])), float(math.exp(b[1])), float(math.exp(b[1])*math.sqrt(cov[1,1]))

def main():
    print(f"eta = {ETA}, c_s = {CS:.5f}, tau_T = {A_HARD:.3f} M (hard disk), {A_IDEAL:.3f} M (ideal)\n")
    out={}
    print("| M | seeds | record | L/tau actual | cal slope | modelled | block | S(0) | spread | in range |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for M in (10,20,50,100,200):
        D,dt=load(M); K=kroot(M/100.0); nu=CS*K/(2*math.pi*LEFF); om=2*math.pi*nu; per=1/nu
        nraw=len(D[0]); rec=nraw*dt
        tau=A_HARD*M; TR0=0.4*per*5; A0,B0=0.50,0.47
        for _ in range(5):
            nl=max(64,int(min(8.0*tau,rec/3)/dt)); fitf=mk_fit(nl,dt,om)
            try: p=fitf(D,[A0,tau,B0,TR0])
            except Exception: break
            if abs(p[1]-tau)<0.02*tau: tau=p[1]; break
            tau=p[1]
        nl=max(64,int(min(8.0*tau,rec/3)/dt)); fitf=mk_fit(nl,dt,om)
        p=fitf(D,[A0,tau,B0,TR0]); raw_m=p[1]; TR0=p[3]; A0,B0=p[0],p[2]
        raw_b=fit_block(D,per,dt)
        rng=np.random.default_rng(9000+M); tr,rc,sd=[],[],[]
        for g in (0.5,0.7,0.9,1.0,1.2,1.5,2.0):
            true=g*raw_m/0.75; got=[]
            for _ in range(6):
                try:
                    v=fitf(synth(true,nraw,len(D),om,dt,TR0,A0,B0,rng),[A0,true,B0,TR0])[1]
                    if np.isfinite(v): got.append(v)
                except Exception: pass
            if len(got)>2: tr.append(true); rc.append(np.mean(got)); sd.append(np.std(got,ddof=1))
        tr,rc,sd=np.array(tr),np.array(rc),np.array(sd)
        slope=float(np.polyfit(tr,rc,1)[0]); inside=bool(rc.min()<=raw_m<=rc.max())
        vm=float(np.interp(raw_m,rc,tr)); em=float(np.interp(raw_m,rc,sd))/max(slope,1e-9)
        vb=float(np.interp(raw_b,rc,tr)) if rc.min()<=raw_b<=rc.max() else float("nan")
        s0,s0e=fit_S0(D,dt,nu,raw_m)
        vals=[v for v in (vm,vb,s0) if np.isfinite(v)]
        out[M]=dict(tau_model=vm,err=em,tau_block=vb,tau_S0=s0,S0err=s0e,slope=slope,
                    rec=rec,Ltau=rec/vm,inside=inside,tau_r=float(p[3]),raw_m=float(raw_m))
        print(f"| {M} | {len(D)} | {rec:.0f} | **{rec/vm:.0f}** | **{slope:.3f}** | {vm:.0f} ± {em:.0f} | "
              f"{vb:.0f} | {s0:.0f} ± {s0e:.0f} | {max(vals)/min(vals):.2f}x | {'yes' if inside else 'NO'} |")
    json.dump({str(k):v for k,v in out.items()},open(os.path.join(HERE,"261005_topup_results.json"),"w"),indent=1)
    print("\n| set | estimator | b ± σ | σ from 1 | a ± σ |")
    print("|---|---|---|---|---|")
    B={}
    for lab,ms in (("all five",[10,20,50,100,200]),("light {10,20,50}",[10,20,50]),("heavy {50,100,200}",[50,100,200])):
        for est,kv in (("modelled","tau_model"),("block","tau_block"),("S(0)","tau_S0")):
            t=[out[m][kv] for m in ms]
            if any(not np.isfinite(v) for v in t): continue
            e=[max(out[m]["err"],0.02*out[m][kv]) for m in ms]
            b,db,a,da=powerlaw(ms,t,e); B[(lab,est)]=(b,db)
            print(f"| {lab} | {est} | {b:.3f} ± {db:.3f} | {abs(b-1)/db:.1f} | {a:.2f} ± {da:.2f} |")
    print("\n### VERDICT RULE, applied as written")
    c1=all(out[m]["Ltau"]>=60 for m in out); c2=all(out[m]["slope"]>=0.6 for m in out)
    print(f"  every mass L/tau >= 60 : {'PASS' if c1 else 'FAIL'}  " +
          ", ".join(f"M{m}={out[m]['Ltau']:.0f}" for m in sorted(out)))
    print(f"  every mass slope >= 0.6: {'PASS' if c2 else 'FAIL'}  " +
          ", ".join(f"M{m}={out[m]['slope']:.3f}" for m in sorted(out)))
    ok3=True
    for est in ("modelled","block","S(0)"):
        kl,kh=("light {10,20,50}",est),("heavy {50,100,200}",est)
        if kl in B and kh in B:
            bl,dl=B[kl]; bh,dh=B[kh]; s=abs(bh-bl)/math.sqrt(dl*dl+dh*dh)
            print(f"  subsets agree within 2σ ({est}): {'PASS' if s<2 else 'FAIL'}  "
                  f"light {bl:.3f}±{dl:.3f} vs heavy {bh:.3f}±{dh:.3f} -> {s:.1f}σ")
            ok3 = ok3 and s < 2
    print(f"\n  ==> b IS {'A RESULT' if (c1 and c2 and ok3) else 'NOT QUOTABLE -- report the failing diagnostic and stop'}")

if __name__ == "__main__":
    main()
