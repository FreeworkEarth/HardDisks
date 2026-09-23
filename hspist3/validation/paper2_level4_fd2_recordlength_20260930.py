import math, numpy as np
from scipy.optimize import curve_fit
rng=np.random.default_rng(47); dt=5.0
OM=0.066025; A0,B0,TR0=0.47,0.50,190.0; NSD=20
def synth(tT,nraw):
    a=math.exp(-dt/tT); s=np.zeros(nraw); s[0]=rng.normal(); e=rng.normal(size=nraw)*math.sqrt(1-a*a)
    for j in range(1,nraw): s[j]=a*s[j-1]+e[j]
    rho=math.exp(-dt/TR0); c1=2*rho*math.cos(OM*dt); c2=-rho*rho
    v=np.zeros(nraw); w=rng.normal(size=nraw)
    for j in range(2,nraw): v[j]=c1*v[j-1]+c2*v[j-2]+w[j]
    v/=(v[200:].std() or 1.0)
    return math.sqrt(A0)*s+math.sqrt(B0)*v
def fitset(series,NL):
    lag=np.arange(NL)*dt; cs=[]
    for x in series:
        x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m); c=np.fft.irfft(f*np.conj(f))[:NL].real; cs.append(c/c[0])
    c=np.mean(cs,axis=0)
    def model(t,A,tT,B,tr): return A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(OM*t)
    try:
        p,_=curve_fit(model,lag,c,p0=[0.5,500.,0.5,200.],bounds=([0,20,0,20],[1.5,1e4,1.5,5e3]),maxfev=40000)
        return p[1]
    except Exception: return np.nan
print("MODELLED estimator (OU + damped cosine, omega fixed) vs record length")
print("| record [sigma] | L/tau | seeds | rec@504 | rec@618 | slope | separation |")
print("|---|---|---|---|---|---|---|")
for mult,ns,ntr in ((1,20,20),(2,20,20),(4,20,20),(4,80,12)):
    nraw=int(8200*mult/dt); NL=int(min(1500*mult,3000)/dt)+1
    out=[]
    for true in (504.,618.):
        g=[]
        for _ in range(ntr):
            r=fitset([synth(true,nraw) for _ in range(ns)],NL)
            if np.isfinite(r): g.append(r)
        g=np.array(g); out.append((g.mean(),g.std(ddof=1)))
    (m1,s1),(m2,s2)=out; sl=(m2-m1)/114.0; sep=(m2-m1)/math.sqrt((s1**2+s2**2)/2)
    print(f"| {8200*mult} | {8200*mult/504:.0f} | {ns} | {m1:.0f}+-{s1:.0f} | {m2:.0f}+-{s2:.0f} | {sl:.2f} | {sep:.1f} |")
