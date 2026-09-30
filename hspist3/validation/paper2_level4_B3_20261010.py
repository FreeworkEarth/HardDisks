#!/usr/bin/env python3
"""##CHRIS 2026-10-10: B3 -- the out-of-equilibrium tau_T test, applied exactly as pre-registered in
261007 section 3B (committed eb7387e before any record existed). Analysis only.
  D(t) = (T_1 - T_2) N k / W_in, seed-mean; the divider mode (period 186 in the compressed box,
  section 4.9) removed by a running mean over ONE period; fit D_0 exp(-t/tau) on t > 3 tau_r;
  verdict: tau within 2 sigma (combined) of the ladder's 40 079 +- 1708."""
import glob, json, math, os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
HERE=os.path.dirname(os.path.abspath(__file__)); REPO=os.path.dirname(os.path.dirname(HERE))
P=os.path.join(REPO,"hspist3","experiments_energy_transfer","level4_B3_20261010","B3")
NS, M = 50, 200.0
TAU_T, TAU_T_ERR = 40079.0, 1708.0          # ladder, M = 200
TAU_R, PERIOD    = 2090.0, 186.2            # mode damping; period after the push (261007 4.9)
D0_REF, D0_ERR   = 0.224, 0.049             # 32-seed cell at ~3 tau_r, consistency line only
fs=sorted(glob.glob(os.path.join(P,"red_*.csv")))
S=[pd.read_csv(f,low_memory=False) for f in fs]
n=min(len(d) for d in S); t=S[0]["Time"].to_numpy(float)[:n]; dt=t[1]-t[0]
W=np.array([d["PistonWork"].to_numpy(float)[n-1] for d in S])
D=np.array([((d["KE_gas_right"].to_numpy(float)[:n]-d["KE_gas_left"].to_numpy(float)[:n])/NS)*NS/W[k] for k,d in enumerate(S)])
# per-seed ledger, exact: W_in - dKE_R - dKE_L - KE_div
led=0.0
for k,d in enumerate(S):
    dKl=d["KE_gas_left"].to_numpy(float)[:n]-d["KE_gas_left"].to_numpy(float)[0]
    dKr=d["KE_gas_right"].to_numpy(float)[:n]-d["KE_gas_right"].to_numpy(float)[0]
    KEd=0.5*M*d["W0_v"].to_numpy(float)[:n]**2; KEd-=KEd[0]
    led=max(led, float(np.abs(d["PistonWork"].to_numpy(float)[:n]-(dKl+dKr+KEd)).max()))
win=max(3,int(round(PERIOD/dt)) | 1)          # one mode period, odd
ker=np.ones(win)/win
Dm=np.convolve(D.mean(0),ker,mode="same")
Dse=D.std(0,ddof=1)/math.sqrt(len(S))
m=(t>3*TAU_R)&(t<t[-1]-PERIOD)                # fit window: after the mechanical stage, inside the smoothed range
f=lambda tt,D0,tau: D0*np.exp(-tt/tau)
p,cov=curve_fit(f,t[m],Dm[m],p0=[0.2,TAU_T],sigma=np.maximum(Dse[m],1e-6),absolute_sigma=True,maxfev=60000)
# jackknife over seeds for an error that respects the time correlation
J=[]
for k in range(len(S)):
    Dk=np.delete(D,k,axis=0).mean(0); Dk=np.convolve(Dk,ker,mode="same")
    try: J.append(curve_fit(f,t[m],Dk[m],p0=p,maxfev=60000)[0])
    except Exception: pass
J=np.array(J); g=len(J); jk=np.sqrt((g-1)/g*((J-J.mean(0))**2).sum(0))
tau,D0=p[1],p[0]; etau,eD0=jk[1],jk[0]
z=abs(tau-TAU_T)/math.hypot(etau,TAU_T_ERR)
print(f"B3: {len(S)} seeds, record {t[-1]:.0f} sigma, dt = {dt:.1f}, running-mean window {win} samples = {win*dt:.0f} sigma")
print(f"W_in = {W.mean():.3f} +- {W.std(ddof=1)/math.sqrt(len(W)):.3f}   per-seed ledger max |residual| = {led:.3e} kT")
print(f"fit window t in [{t[m][0]:.0f}, {t[m][-1]:.0f}]  ({m.sum()} points)")
print(f"\n| quantity | fitted | reference | sigma |")
print(f"|---|---|---|---|")
print(f"| **tau_T** | **{tau:.0f} +- {etau:.0f}** | 40 079 +- 1708 (ladder) | **{z:.2f}** |")
print(f"| tau_T vs Gruber-Piasecki | | 10 070 | {abs(tau-10070)/etau:.1f} |")
print(f"| D_0 | {D0:.4f} +- {eD0:.4f} | 0.224 +- 0.049 (32-seed cell, ~3 tau_r) | {abs(D0-D0_REF)/math.hypot(eD0,D0_ERR):.1f} (consistency, not a gate) |")
for tt in (6300,10000,20000,40000,80000,120000):
    i=min(int(np.searchsorted(t,tt)),n-1); print(f"   D({t[i]:6.0f}) = {Dm[i]:+.4f} +- {Dse[i]:.4f}")
print(f"\n### VERDICT: fitted tau_T is {z:.2f} sigma from the ladder's 40 079  ->  **{'CONFIRMED out of equilibrium' if z<2 else 'NOT within 2 sigma -- reported as measured'}**")
json.dump(dict(n=len(S),tau=tau,tau_err=etau,D0=D0,D0_err=eD0,z=z,W=float(W.mean()),ledger_max=led,window=[float(t[m][0]),float(t[m][-1])]),
          open(os.path.join(HERE,"261010_B3_results.json"),"w"),indent=1)
