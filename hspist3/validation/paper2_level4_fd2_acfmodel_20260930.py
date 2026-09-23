import glob, math, os, sys, numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
sys.path.insert(0,os.path.abspath('hspist3/validation')); sys.path.insert(0,os.path.abspath('hspist3'))
import plot_speed_of_sound_edmd as sos
P="hspist3/experiments_energy_transfer/level4_equilibrium_20260929/Md10"; NS=50; dt=5.0
eta=50*math.pi*0.25/(39.25*10)
Z=float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0]); dZ=float(sos.dZ_kolafa_rottner_2006(np.array([eta]))[0])
CS=math.sqrt(Z+eta*dZ+Z*Z); K=brentq(lambda k: math.cos(k)/math.sin(k)-0.1*k,1e-9,math.pi-1e-9)
LEFF=38.75-1.0; NU=CS*K/(2*math.pi*LEFF); OM=2*math.pi*NU
print(f"predicted mode: K={K:.6f} cs={CS:.6f} Leff={LEFF} nu={NU:.6f} period={1/NU:.2f} omega={OM:.6f}")
fs=sorted(glob.glob(f"{P}/red_*.csv")); Ds=[];Xs=[]
for f in fs:
    e=pd.read_csv(f)
    Ds.append((e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))/NS)
    Xs.append(e["W0_x_sigma"].to_numpy(float))
n=min(len(a) for a in Ds); lo=int(2000/dt)
def acf1(x,nl):
    x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m); c=np.fft.irfft(f*np.conj(f))[:nl].real; return c/c[0]
NL=int(1500/dt)+1
def meanacf(S,idx=None):
    idx=range(len(S)) if idx is None else idx
    return np.mean([acf1(S[i][lo:n],NL) for i in idx],axis=0)
lag=np.arange(NL)*dt
def model(t,A,tT,B,tr): return A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(OM*t)
def fit(c):
    try:
        p,_=curve_fit(model,lag,c,p0=[0.5,500.,0.5,200.],
                      bounds=([0,20,0,20],[1.5,1e4,1.5,5e3]),maxfev=40000)
        return p
    except Exception: return np.full(4,np.nan)
print("\n=== MODEL FIT, omega FIXED from cot K = alpha K ===")
res={}
for nm,S in (("T1-T2",Ds),("divider x",Xs)):
    c=meanacf(S); p=fit(c)
    js=np.array([fit(meanacf(S,[j for j in range(len(S)) if j!=i])) for i in range(len(S))])
    k=len(js); err=np.sqrt((k-1)/k*np.nansum((js-np.nanmean(js,axis=0))**2,axis=0))
    res[nm]=(p,err)
    print(f"{nm:10s}: A={p[0]:.3f}+-{err[0]:.3f}  tau_T={p[1]:.0f}+-{err[1]:.0f}  B={p[2]:.3f}+-{err[2]:.3f}  tau_r={p[3]:.0f}+-{err[3]:.0f}")
    rms=np.sqrt(np.mean((c-model(lag,*p))**2)); print(f"            fit rms residual = {rms:.4f}")
# omega free, as an independent measurement of the period
def model2(t,A,tT,B,tr,w): return A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(w*t)
print("\n=== same fit with omega FREE (independent period measurement) ===")
for nm,S in (("T1-T2",Ds),("divider x",Xs)):
    c=meanacf(S)
    p,_=curve_fit(model2,lag,c,p0=[0.5,500.,0.5,200.,OM],bounds=([0,20,0,20,0.03],[1.5,1e4,1.5,5e3,0.12]),maxfev=60000)
    js=[]
    for i in range(len(S)):
        try:
            q,_=curve_fit(model2,lag,meanacf(S,[j for j in range(len(S)) if j!=i]),p0=p,
                          bounds=([0,20,0,20,0.03],[1.5,1e4,1.5,5e3,0.12]),maxfev=60000); js.append(q)
        except Exception: pass
    js=np.array(js); k=len(js); e=np.sqrt((k-1)/k*((js-js.mean(0))**2).sum(0))
    print(f"{nm:10s}: period = {2*math.pi/p[4]:.2f} +- {2*math.pi/p[4]**2*e[4]:.2f}   tau_T={p[1]:.0f}+-{e[1]:.0f}  tau_r={p[3]:.0f}+-{e[3]:.0f}")
np.save('/tmp/l4fd2_res.npy',np.array([res[k][0] for k in ("T1-T2","divider x")]))
