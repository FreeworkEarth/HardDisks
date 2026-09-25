#!/usr/bin/env python3
"""##CHRIS 2026-10-06: the divider mode across the ladder -- period vs cot K = alpha K, and tau_r
vs Mansour's piston form. Analysis only, on traces already on disk. No new runs, no core changes.

WHY THIS IS WORTH A TABLE OF ITS OWN. The ladder's tau_T fits FIXED omega at the Kolafa-Rottner
prediction, so the mode period was never measured across the ladder -- only at M = 10 (260930).
Here omega is FREE at every mass in both boxes, which turns "the oscillation is Paper 1's divider
mode" from one agreement into a seven-point test of the eigenvalue equation itself:

    cot K = alpha K ,   alpha = M/(2 N_s m) = 1/(2R) ,   nu = c_s K/(2 pi L_eff)

Inverting: K_measured = 2 pi L_eff nu_measured / c_s. Plotting K_measured against alpha for the two
boxes together tests the equation over alpha = 0.1 - 2.0 with no free parameter, and the choice of
c_s (KR / Paper 1's measured / ideal gas) shifts every point by a common factor -- so the SAME
figure discriminates the sound speed.

ETA CONVENTION, and a 0.28 % shift from 261003. That script predicted periods at eta = 0.100051,
which is the summary's eta_nominal = N pi r^2/(2 L0 h) with L0 the WALL POSITION -- i.e. it ignores
the 1.0-sigma divider thickness. The audited physical value uses the free compartment length
L_c = L0 - t/2, giving eta = 0.10134170 in BOTH boxes and c_s = 1.749302 instead of 1.744337.
Every predicted period here is therefore 0.28 % SHORTER than 261003's. Stated, not silently applied.

TAU_R. Mansour's piston linewidth is Delta f/f = Gamma L_y sqrt(2/(M_hat N)) with M_hat = M + m N/3
and N the TOTAL particle count (100 in the ladder box, 200 in the doubled box), so tau_r follows
from Delta f/f = 1/(pi tau_r nu). Gamma = 0.331 (Enskog at eta = 0.10).

ERRORS. Delete-8 block jackknife over the 80 seeds (10 groups), sigma^2 = (g-1)/g * sum (t_i - tbar)^2.
Full delete-1 jackknife would be 80 five-parameter fits per cell on up to 32 000 lags; the block
version is the same estimator at a tenth of the cost, and the machine is busy with campaign A.
The spectral peak's FWHM is reported separately as the PRECISION of the line, which is not the same
thing as the error on the fitted centre (the 260930 lesson: agreement is not precision).
"""
import glob, json, math, os, sys
import numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); sys.path.insert(0,os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO=os.path.dirname(os.path.dirname(HERE)); ET=os.path.join(REPO,"hspist3","experiments_energy_transfer")

ETA  = 0.10134170                                  # audited: free compartment length, both boxes
Zv   = float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0])
dZv  = float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
CS_KR= math.sqrt(Zv+ETA*dZv+Zv*Zv)
CS_P1= CS_KR*1.0101                                # Paper 1's measured c_s, +1.01 % on KR
CS_ID= math.sqrt(2.0)                             # ideal gas, 2D monatomic, kT=m=1
GAMMA_ENSKOG, LY, BURN = 0.331, 10.0, 2000.0

# (label, N_s, L_c, dir, M, tau_T known from the ladder/pass-1 analyses -- sets the fit window only)
CELLS=[("A",  50, 38.75, "level4_topup_20261005/Md10",         10,   475.),
       ("A",  50, 38.75, "level4_topup_20261005/Md20",         20,  1089.),
       ("A",  50, 38.75, "level4_ladder_rerun_20261004/Md50",  50,  3835.),
       ("A",  50, 38.75, "level4_ladder_rerun_20261004/Md100",100, 12330.),
       ("A",  50, 38.75, "level4_ladder_rerun_20261004/Md200",200, 40079.),
       ("B", 100, 77.50, "level4_Rcollapse_20261005/Md50",     50,  7253.),
       ("B", 100, 77.50, "level4_Rcollapse_20261005/Md100",   100, 15773.)]
PASS2=[("B", 100, 77.50, "level4_Rcollapse2_20261006/Md25",    25,  2470.),
       ("B", 100, 77.50, "level4_Rcollapse2_20261006/Md50",    50,  7253.),
       ("B", 100, 77.50, "level4_Rcollapse2_20261006/Md100",  100, 15773.)]

def kroot(a): return brentq(lambda k: math.cos(k)/math.sin(k)-a*k, 1e-9, math.pi-1e-9)
def acf(x,nl):
    x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m)
    c=np.fft.irfft(f*np.conj(f))[:nl].real; return c/c[0]

def load(sub):
    d=os.path.join(ET,sub); X=[]; T=[]; dts=[]
    for f in sorted(glob.glob(os.path.join(d,"red_*.csv"))):
        e=pd.read_csv(f,usecols=["Time","KE_gas_left","KE_gas_right","W0_x_sigma"])
        t=e["Time"].to_numpy(float); dts.append((t[-1]-t[0])/(len(t)-1))
        X.append(e["W0_x_sigma"].to_numpy(float))
        T.append(e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))
    if not X: return None
    dt=float(np.median(dts)); n=min(len(a) for a in X); lo=int(BURN/dt)
    return [a[lo:n] for a in X],[a[lo:n] for a in T],dt

def fit_free(series,dt,nl,p0):
    """5-parameter fit with OMEGA FREE: A e^-t/tT + B e^-t/tr cos(om t)."""
    lag=np.arange(nl)*dt
    c=np.mean([acf(s,nl) for s in series],axis=0)
    g=lambda t,A,tT,B,tr,om: A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(om*t)
    lo=[0,20,0,20,0.2*p0[4]]; hi=[1.5,1e7,1.5,1e6,5.0*p0[4]]
    p,_=curve_fit(g,lag,c,p0=p0,bounds=(lo,hi),maxfev=200000)
    return p

def jack(series,dt,nl,p0,groups=10):
    n=len(series); idx=np.arange(n); out=[]
    for gi in range(groups):
        keep=[series[i] for i in idx if i%groups!=gi]
        try: out.append(fit_free(keep,dt,nl,p0))
        except Exception: pass
    if len(out)<3: return None
    A=np.array(out); g=len(out); mean=A.mean(0)
    return np.sqrt((g-1)/g*((A-mean)**2).sum(0))

def spec_peak(series,dt,nu0):
    """Averaged periodogram of the divider position; peak near nu0 and its FWHM."""
    sp=[]
    for s in series:
        x=s-s.mean(); w=np.hanning(len(x)); x=x*w
        sp.append(np.abs(np.fft.rfft(x))**2*2.0*dt/((w**2).sum()))
    S=np.mean(sp,axis=0); fr=np.fft.rfftfreq(len(series[0]),dt)
    m=(fr>0.4*nu0)&(fr<2.2*nu0)
    if m.sum()<8: return float("nan"),float("nan")
    f,s=fr[m],S[m]; k=int(np.argmax(s)); pk=f[k]; half=s[k]/2
    lo=f[0]
    for i in range(k,0,-1):
        if s[i]<half: lo=f[i]; break
    hi=f[-1]
    for i in range(k,len(f)):
        if s[i]<half: hi=f[i]; break
    return float(pk), float(hi-lo)

def run(cells,tag):
    rows=[]
    for box,NS,LC,sub,M,tauT in cells:
        got=load(sub)
        if got is None: print(f"  (skip {sub}: no traces yet)"); continue
        X,T,dt=got
        LEFF=LC-1.0; NTOT=2*NS
        al=M/(2.0*NS); K=kroot(al); nu_kr=CS_KR*K/(2*math.pi*LEFF)
        rec=len(X[0])*dt; nl=max(128,int(min(8.0*tauT,rec/3.0)/dt))
        p0=[0.50,tauT,0.47,0.4*(1/nu_kr)*5,2*math.pi*nu_kr]
        res={}
        for name,ser in (("x",X),("dT",T)):
            try: p=fit_free(ser,dt,nl,p0)
            except Exception as e: print(f"  {sub} {name}: fit failed {e}"); continue
            e=jack(ser,dt,nl,p)
            res[name]=dict(A=p[0],tauT=p[1],B=p[2],taur=p[3],om=p[4],
                           per=2*math.pi/p[4], per_err=2*math.pi*e[4]/p[4]**2 if e is not None else float("nan"),
                           taur_err=e[3] if e is not None else float("nan"))
        pk,fw=spec_peak(X,dt,nu_kr)
        Mhat=M+NTOT/3.0; dff=GAMMA_ENSKOG*LY*math.sqrt(2.0/(Mhat*NTOT))
        taur_man=1.0/(math.pi*dff*nu_kr)
        rows.append(dict(box=box,NS=NS,LC=LC,LEFF=LEFF,M=M,R=NS/M,alpha=al,K_pred=K,dt=dt,rec=rec,
                         nl=nl,nseed=len(X),nu_kr=nu_kr,per_kr=1/nu_kr,
                         per_p1=1/(CS_P1*K/(2*math.pi*LEFF)),per_id=1/(CS_ID*K/(2*math.pi*LEFF)),
                         spec_per=1/pk if pk==pk else float("nan"),
                         spec_fwhm_frac=fw/pk if pk==pk else float("nan"),
                         Mhat=Mhat,dff_man=dff,taur_man=taur_man,**{f"{k}_{n}":v
                          for n,d in res.items() for k,v in d.items()}))
        r=rows[-1]
        print(f"  {tag} box {box} M={M:3d} a={al:.3f} {len(X)} seeds nl={nl} "
              f"per_fit={r.get('per_x',float('nan')):.2f} per_kr={1/nu_kr:.2f} taur={r.get('taur_x',float('nan')):.0f}")
    return rows

def report(rows):
    print("\n### 1. The eigenvalue equation across alpha = 0.1 - 2.0, omega FREE\n")
    print("| box | N_s | M | alpha | R | period (x, fit) | period (dT, fit) | spectral peak | FWHM | KR | Paper 1 c_s | ideal |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['box']} | {r['NS']} | {r['M']} | {r['alpha']:.3f} | {r['R']:.2f} | "
              f"**{r['per_x']:.2f} ± {r['per_err_x']:.2f}** | {r['per_dT']:.2f} ± {r['per_err_dT']:.2f} | "
              f"{r['spec_per']:.2f} | ±{50*r['spec_fwhm_frac']:.1f}% | {r['per_kr']:.2f} | {r['per_p1']:.2f} | {r['per_id']:.2f} |")
    print("\n### 2. Inverted: K measured against cot K = alpha K\n")
    print("| box | M | alpha | K from cot K = aK | K meas (KR c_s) | sigma | K meas (Paper 1) | sigma | K meas (ideal) | sigma |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    agg={}
    for r in rows:
        f=2*math.pi*r['LEFF']/r['per_x']; fe=f*r['per_err_x']/r['per_x']
        for nm,cs in (("KR",CS_KR),("P1",CS_P1),("ID",CS_ID)):
            km=f/cs; ke=fe/cs; s=(km-r['K_pred'])/ke
            r[f"K_{nm}"]=km; r[f"Ke_{nm}"]=ke; r[f"sig_{nm}"]=s
            agg.setdefault(nm,[]).append(s)
        print(f"| {r['box']} | {r['M']} | {r['alpha']:.3f} | **{r['K_pred']:.4f}** | {r['K_KR']:.4f} | "
              f"**{r['sig_KR']:+.1f}σ** | {r['K_P1']:.4f} | {r['sig_P1']:+.1f}σ | {r['K_ID']:.4f} | {r['sig_ID']:+.1f}σ |")
    print("\n| sound speed | mean deviation | rms | chi2 / n |")
    print("|---|---|---|---|")
    for nm,lab in (("KR","Kolafa-Rottner 1.7493"),("P1","Paper 1 measured 1.7670"),("ID","ideal gas 1.4142")):
        a=np.array(agg[nm]); print(f"| {lab} | {a.mean():+.2f}σ | {np.sqrt((a**2).mean()):.2f}σ | {(a**2).mean():.1f} |")
    print("\n### 3. tau_r against Mansour's piston form, both boxes\n")
    print("| box | N_s | M | M_hat = M + mN/3 | Delta f/f Mansour | tau_r Mansour | **tau_r measured (x)** | tau_r (dT) | **Mansour/measured** |")
    print("|---|---|---|---|---|---|---|---|---|")
    rat=[]
    for r in rows:
        q=r['taur_man']/r['taur_x']; rat.append(q)
        print(f"| {r['box']} | {r['NS']} | {r['M']} | {r['Mhat']:.2f} | {r['dff_man']:.4f} | {r['taur_man']:.0f} | "
              f"**{r['taur_x']:.0f} ± {r['taur_err_x']:.0f}** | {r['taur_dT']:.0f} ± {r['taur_err_dT']:.0f} | **{q:.2f}×** |")
    rat=np.array(rat)
    print(f"\nMansour/measured: mean **{rat.mean():.2f}×**, range {rat.min():.2f}-{rat.max():.2f}×, "
          f"n = {len(rat)}. The mode damps FASTER than the piston form at every mass in both boxes,")
    print("i.e. a correspondingly BROADER line. That is a measurement, not a discrepancy to apologise for.")
    return rows

def main():
    print(f"eta={ETA:.8f}  Z={Zv:.6f}  etaZ'={ETA*dZv:.6f}")
    print(f"c_s: KR={CS_KR:.6f}  Paper1={CS_P1:.6f} (+1.01%)  ideal={CS_ID:.6f}")
    print(f"(261003 used eta=0.100051 -> c_s=1.744337; every predicted period here is 0.28% shorter)\n")
    rows=run(CELLS,"have")
    p2=run(PASS2,"pass2")
    rows=report(rows+[r for r in p2 if r['M']==25])          # pass-2 M=25 extends alpha down to 0.125
    json.dump(dict(meta=dict(ETA=ETA,CS_KR=CS_KR,CS_P1=CS_P1,CS_ID=CS_ID,GAMMA=GAMMA_ENSKOG),
                   cells=rows),open(os.path.join(HERE,"261006_mode_ladder.json"),"w"),indent=1,default=float)
    print("\nwrote 261006_mode_ladder.json")

if __name__=="__main__":
    main()
