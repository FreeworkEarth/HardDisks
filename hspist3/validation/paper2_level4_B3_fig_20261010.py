#!/usr/bin/env python3
"""##CHRIS 2026-10-10: figure only for B3 -- reads the same traces and the fit result JSON written by
paper2_level4_B3_20261010.py (the pre-registered fit is not touched here)."""
import glob, json, math, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); REPO=os.path.dirname(os.path.dirname(HERE))
P=os.path.join(REPO,"hspist3","experiments_energy_transfer","level4_B3_20261010","B3")
OUT=os.path.join(REPO,"0000_PLAN_OVERALL","paper2_energytransfer","experiments","final")
BLUE,RED,BLACK,GREY="#1f4e9c","#c0392b","#000000","#7f7f7f"
J=json.load(open(os.path.join(HERE,"261010_B3_results.json")))
NS=50; PERIOD=186.2
S=[pd.read_csv(f,low_memory=False) for f in sorted(glob.glob(os.path.join(P,"red_*.csv")))]
n=min(len(d) for d in S); t=S[0]["Time"].to_numpy(float)[:n]; dt=t[1]-t[0]
W=np.array([d["PistonWork"].to_numpy(float)[n-1] for d in S])
D=np.array([(d["KE_gas_right"].to_numpy(float)[:n]-d["KE_gas_left"].to_numpy(float)[:n])/W[k] for k,d in enumerate(S)])
win=max(3,int(round(PERIOD/dt))|1); ker=np.ones(win)/win
Dm=np.convolve(D.mean(0),ker,mode="same"); Dse=np.convolve(D.std(0,ddof=1),ker,mode="same")/math.sqrt(len(S))
fig,ax=plt.subplots(figsize=(8.2,4.6))
ax.fill_between(t,Dm-Dse,Dm+Dse,color=BLUE,alpha=.18,lw=0,label=r"80-seed mean $\pm$ sem (one-period running mean)")
ax.plot(t,Dm,color=BLUE,lw=1.3)
tt=np.linspace(J["window"][0],t[-1],300)
ax.plot(tt,J["D0"]*np.exp(-tt/J["tau"]),color=BLACK,lw=2,label=rf"fit: $\tau_T={J['tau']:.0f}\pm{J['tau_err']:.0f}$")
ax.plot(tt,J["D0"]*np.exp(-tt/40079.),color=RED,lw=1.5,ls="--",label=r"ladder $\tau_T=40\,079$ (same $D_0$)")
ax.plot(tt,J["D0"]*np.exp(-tt/10070.),color=GREY,lw=1.3,ls=":",label=r"Gruber–Piasecki $\tau_T=10\,070$")
ax.axvspan(0,J["window"][0],color=GREY,alpha=.12,lw=0,label=r"mechanical stage ($<3\tau_r$), not fitted")
ax.axhline(0,color=GREY,lw=.7)
ax.set_xlabel(r"$t$ after the push  [$\sigma$-time]"); ax.set_ylabel(r"$D=(T_1-T_2)Nk/W_{\rm in}$")
ax.set_title(r"B3: the thermal stage watched — $u=1.0$, $M_d=200$, 80 seeds, binary 05215ea",fontsize=10)
ax.set_xlim(0,t[-1]); ax.legend(fontsize=8,loc="upper right"); ax.grid(alpha=.25)
fig.tight_layout()
for ext in ("png","pdf"):
    fp=os.path.join(OUT,f"261010_p2_B3_tauT.{ext}"); fig.savefig(fp,dpi=180)
print(os.path.join(OUT,"261010_p2_B3_tauT.png"))
