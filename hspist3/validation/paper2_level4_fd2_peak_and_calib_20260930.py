import glob, math, os, sys, numpy as np, pandas as pd
from scipy.optimize import curve_fit, brentq
sys.path.insert(0,os.path.abspath('hspist3/validation')); sys.path.insert(0,os.path.abspath('hspist3'))
import plot_speed_of_sound_edmd as sos
P="hspist3/experiments_energy_transfer/level4_equilibrium_20260929/Md10"; NS=50; dt=5.0
eta=50*math.pi*0.25/(39.25*10)
Z=float(sos.Z_kolafa_rottner_2006(np.array([eta]))[0]); dZ=float(sos.dZ_kolafa_rottner_2006(np.array([eta]))[0])
CS=math.sqrt(Z+eta*dZ+Z*Z); K=brentq(lambda k: math.cos(k)/math.sin(k)-0.1*k,1e-9,math.pi-1e-9)
LEFF=37.75; NU=CS*K/(2*math.pi*LEFF); OM=2*math.pi*NU
fs=sorted(glob.glob(f"{P}/red_*.csv")); Ds=[];Xs=[]
for f in fs:
    e=pd.read_csv(f)
    Ds.append((e["KE_gas_left"].to_numpy(float)-e["KE_gas_right"].to_numpy(float))/NS)
    Xs.append(e["W0_x_sigma"].to_numpy(float))
n=min(len(a) for a in Ds); lo=int(2000/dt)
print("=== SPECTRAL PEAK, with width as the error ===")
def psd(s):
    s=s[lo:n]; s=(s-s.mean())*np.hanning(len(s))
    return np.fft.rfftfreq(len(s),dt), np.abs(np.fft.rfft(s))**2
fr,_=psd(Ds[0])
for nm,S in (("T1-T2",Ds),("divider x",Xs)):
    p=np.mean([psd(a)[1] for a in S],axis=0)
    band=(fr>0.004)&(fr<0.020); i=np.argmax(np.where(band,p,-1))
    # parabolic interpolation of the peak
    d=0.5*(p[i-1]-p[i+1])/(p[i-1]-2*p[i]+p[i+1]); fpk=fr[i]+d*(fr[1]-fr[0])
    half=p[i]/2; j=i
    while j>0 and p[j]>half: j-=1
    k2=i
    while k2<len(p)-1 and p[k2]>half: k2+=1
    fwhm=fr[k2]-fr[j]; sig=fwhm/2.355
    print(f"{nm:10s}: f_peak={fpk:.6f}  period={1/fpk:.2f}  FWHM={fwhm:.6f} -> period err ~ {sig/fpk**2:.2f}")
print(f"\npredicted (KR)            : {1/NU:.2f}")
print(f"predicted (ideal)         : {2*math.pi*LEFF/(math.sqrt(2.0)*K):.2f}")
print(f"predicted (Paper1 cs +1.01%): {1/(CS*1.0101*K/(2*math.pi*LEFF)):.2f}")
NL=int(1500/dt)+1; lag=np.arange(NL)*dt
def acf1(x):
    x=x-x.mean(); m=len(x); f=np.fft.rfft(x,2*m); c=np.fft.irfft(f*np.conj(f))[:NL].real; return c/c[0]
def model(t,A,tT,B,tr): return A*np.exp(-t/tT)+B*np.exp(-t/tr)*np.cos(OM*t)
def fitset(series):
    c=np.mean([acf1(s) for s in series],axis=0)
    try:
        p,_=curve_fit(model,lag,c,p0=[0.5,500.,0.5,200.],bounds=([0,20,0,20],[1.5,1e4,1.5,5e3]),maxfev=40000)
        return p[1]
    except Exception: return np.nan
print("\n=== BIAS CALIBRATION of the MODELLED estimator (OU + damped cosine) ===")
rng=np.random.default_rng(31); nraw=n-lo; NSD=len(fs)
A0,B0,TR0=0.47,0.50,190.0
def synth(tT,nraw):
    a=math.exp(-dt/tT); s=np.zeros(nraw); s[0]=rng.normal(); e=rng.normal(size=nraw)*math.sqrt(1-a*a)
    for j in range(1,nraw): s[j]=a*s[j-1]+e[j]
    rho=math.exp(-dt/TR0); c1=2*rho*math.cos(OM*dt); c2=-rho*rho
    v=np.zeros(nraw); w=rng.normal(size=nraw)
    for j in range(2,nraw): v[j]=c1*v[j-1]+c2*v[j-2]+w[j]
    v/= (v[200:].std() or 1.0)
    return math.sqrt(A0)*s+math.sqrt(B0)*v
print("| true tau_T | recovered (mean of 25) | spread |")
print("|---|---|---|")
TRUE=[250,350,450,504,550,618,700,900]; REC=[]
for t in TRUE:
    g=[]
    for _ in range(25):
        r=fitset([synth(float(t),nraw) for _ in range(NSD)])
        if np.isfinite(r): g.append(r)
    g=np.array(g); REC.append(g.mean()); print(f"| {t} | {g.mean():.0f} | {g.std(ddof=1):.0f} |")
TRUE=np.array(TRUE,float); REC=np.array(REC)
sl=np.polyfit(TRUE,REC,1)[0]; print(f"\nresponse slope = {sl:.3f}  (block estimator was 0.335)")
for nm,mv,ev in (("T1-T2",311.,34.),("divider x",316.,38.)):
    inv=np.interp(mv,REC,TRUE); print(f"{nm:10s}: recovered {mv:.0f} -> true tau_T = {inv:.0f} +- {ev/sl:.0f}")
    for lbl,pr in (("hard disk 503",503.),("ideal 618",618.)):
        cr=np.interp(pr,TRUE,REC); sd=np.interp(pr,TRUE,[0]*0+list(np.zeros(len(TRUE)))) if False else None
        print(f"     vs {lbl}: control recovers {cr:.0f}, measured {mv:.0f} -> {(mv-cr)/ev:+.1f} sigma")
