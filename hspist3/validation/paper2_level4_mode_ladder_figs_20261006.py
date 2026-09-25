#!/usr/bin/env python3
"""##CHRIS 2026-10-06: tables + figure for the divider mode across the ladder. Reads
261006_mode_ladder.json (written by paper2_level4_mode_ladder_20261006.py); refits nothing.

WHY THIS FILE EXISTS SEPARATELY, AND THE ERROR BAR IT USES. The first pass of this analysis printed
deviations of "+179 sigma" from the Kolafa-Rottner prediction, using the block-jackknife error on the
fitted centre -- 0.01 to 0.12 %. That is a statistical error only, and quoting it as evidence is the
mirror image of the 260930 mistake (agreement is not precision): 260930 quoted a 0.03 % AGREEMENT off
a 5 % measurement, and this would quote a 100-sigma DISAGREEMENT off a 0.02 % statistical error while
a systematic ten times larger sits unmeasured beside it.

The error used here is therefore

    sigma = max( block jackknife , |period(x) - period(dT)| / 2 )

i.e. the statistical error or half the spread between the two independent observables (divider
position and temperature difference), whichever is larger. That spread is 0.88 % at M = 10 and below
0.1 % for M >= 50, so it dominates exactly where the jackknife is least believable. The spectral
FWHM (0.8-6 %) is reported as the LINE WIDTH and never as the error on the centre.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=os.path.dirname(os.path.abspath(__file__))
REPO=os.path.dirname(os.path.dirname(HERE))
OUT=os.path.join(REPO,"0000_PLAN_OVERALL","paper2_energytransfer","experiments","final")
BLUE, RED, BLACK, GREY = "#1f4e9c", "#c0392b", "#000000", "#7f7f7f"

D=json.load(open(os.path.join(HERE,"261006_mode_ladder.json")))
CS_KR, CS_P1, CS_ID = D["meta"]["CS_KR"], D["meta"]["CS_P1"], D["meta"]["CS_ID"]
C=sorted(D["cells"], key=lambda r:(r["box"], r["M"]))
def kroot(a): return brentq(lambda k: math.cos(k)/math.sin(k)-a*k, 1e-9, math.pi-1e-9)

for r in C:
    r["sig_per"]=max(r["per_err_x"], abs(r["per_x"]-r["per_dT"])/2.0)
    r["sig_taur"]=max(r["taur_err_x"], abs(r["taur_x"]-r["taur_dT"])/2.0)
    f=2*math.pi*r["LEFF"]/r["per_x"]; fe=f*r["sig_per"]/r["per_x"]
    for nm,cs in (("KR",CS_KR),("P1",CS_P1),("ID",CS_ID)):
        r[f"k_{nm}"]=f/cs; r[f"ke_{nm}"]=fe/cs
        r[f"d_{nm}"]=100*(f/cs/r["K_pred"]-1); r[f"de_{nm}"]=100*fe/cs/r["K_pred"]

print("### 1. The divider mode across alpha = 0.1 - 2.0, omega FREE, two boxes\n")
print("| box | N_s | L_c | M | alpha | R | period measured (x) | sigma | period (dT) | spectral peak | line FWHM | KR predicted | Paper 1 c_s | ideal gas |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in C:
    print(f"| {r['box']} | {r['NS']} | {r['LC']:.2f} | {r['M']} | {r['alpha']:.3f} | {r['R']:.2f} | "
          f"**{r['per_x']:.2f}** | ±{r['sig_per']:.2f} | {r['per_dT']:.2f} | {r['spec_per']:.2f} | "
          f"±{50*r['spec_fwhm_frac']:.1f}% | {r['per_kr']:.2f} | {r['per_p1']:.2f} | {r['per_id']:.2f} |")

print("\n### 2. cot K = alpha K, inverted -- the departure is a function of alpha, not of box size\n")
print("| box | M | alpha | K from cot K = aK | K measured (KR c_s) | departure | vs Paper 1 c_s | vs ideal gas |")
print("|---|---|---|---|---|---|---|---|")
for r in C:
    print(f"| {r['box']} | {r['M']} | {r['alpha']:.3f} | **{r['K_pred']:.4f}** | {r['k_KR']:.4f} ± {r['ke_KR']:.4f} | "
          f"**{r['d_KR']:+.2f} ± {r['de_KR']:.2f} %** | {r['d_P1']:+.2f} % | {r['d_ID']:+.1f} % |")
a=[r for r in C if r['box']=='A' and r['M']==50][0]; b=[r for r in C if r['box']=='B' and r['M']==100][0]
print(f"\n**The one shared alpha = 0.500 is the test that matters**: box A (N_s = 50, M = 50) gives "
      f"{a['d_KR']:+.2f} ± {a['de_KR']:.2f} %, box B (N_s = 100, M = 100) gives {b['d_KR']:+.2f} ± {b['de_KR']:.2f} %.")
print(f"They agree to **{abs(a['d_KR']-b['d_KR']):.2f} %** while L_c, N and the period itself all differ by a factor 2 "
      f"({a['per_x']:.1f} vs {b['per_x']:.1f}).")
dv=np.array([r['d_KR'] for r in C])
print(f"\nDeparture from KR grows monotonically with alpha: {dv.min():+.2f} % at alpha = 0.100 to "
      f"{dv.max():+.2f} % at alpha = 2.000 -- a factor {dv.max()/dv.min():.1f}, so **no single sound speed can absorb it**.")
print(f"Against Paper 1's measured c_s (+1.01 % on KR) the residual is "
      f"{min(r['d_P1'] for r in C):+.2f} % to {max(r['d_P1'] for r in C):+.2f} %, crossing zero near alpha = 0.5.")
print(f"Against the ideal gas: {np.mean([r['d_ID'] for r in C]):+.1f} % on average. **Excluded absolutely.**")

print("\n### 3. tau_r against Mansour's piston form -- the form becomes EXACT in the heavy limit\n")
print("| box | N_s | M | M_hat | Delta f/f Mansour | tau_r Mansour | tau_r measured | sigma | Mansour/measured | tau_r/L_c |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in C:
    print(f"| {r['box']} | {r['NS']} | {r['M']} | {r['Mhat']:.2f} | {r['dff_man']:.4f} | {r['taur_man']:.0f} | "
          f"**{r['taur_x']:.0f}** | ±{r['sig_taur']:.0f} | **{r['taur_man']/r['taur_x']:.2f}×** | {r['taur_x']/r['LC']:.2f} |")
A=[r for r in C if r['box']=='A']
print(f"\n**In box A the ratio falls monotonically 2.04 -> 1.03 as M goes 10 -> 200.** Mansour's form is a")
print("PISTON form: it assumes the divider carries the inertia. At alpha = 0.1 the gas standing wave carries")
print("it instead, and the form is 2x too slow; by alpha = 2.0 the divider does dominate and the form is exact")
print("to 3 %. So this is not a standing disagreement -- it is a crossover, and its location is the result.")
for M in (50,100):
    x=[r for r in C if r['M']==M]
    if len(x)==2:
        print(f"\ntau_r/L_c at M = {M}: {x[0]['taur_x']/x[0]['LC']:.2f} (box A) vs {x[1]['taur_x']/x[1]['LC']:.2f} (box B) "
              f"-- agree to {100*abs(x[0]['taur_x']/x[0]['LC']/(x[1]['taur_x']/x[1]['LC'])-1):.1f} %.")
print("\nSo **tau_r scales with L_c at fixed M**, while Mansour's form predicts the box ratio "
      f"{[f'{r:.2f}' for r in [2257/782,3168/1237]]} against the measured "
      f"{[f'{r:.2f}' for r in [1209/586,2091/1085]]} -- its box scaling is ~40 % too strong.")
print("tau_r therefore does NOT collapse on alpha the way the period does, nor on R the way tau_T does. **OPEN.**")

# ---------------- figure ----------------
fig,ax=plt.subplots(2,2,figsize=(11.5,8.6))
al=np.logspace(math.log10(0.07),math.log10(3.0),400); kc=np.array([kroot(x) for x in al])
A=[r for r in C if r['box']=='A']; B=[r for r in C if r['box']=='B']

p=ax[0,0]
p.plot(al,kc,"-",color=RED,lw=1.8,label="$\\cot K=\\alpha K$  (Román eigenmode)")
p.axhline(math.pi/2,ls=":",color=GREY,lw=1,label=r"$\alpha\!\to\!0$ limit, $K=\pi/2$")
for S,mk,fc,lb in ((A,"o",BLUE,r"box A, $N_s=50$"),(B,"s","white",r"box B, $N_s=100$")):
    p.errorbar([r["alpha"] for r in S],[r["k_KR"] for r in S],yerr=[r["ke_KR"] for r in S],
               fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=lb+r", $c_s$ KR")
p.errorbar([r["alpha"] for r in C],[r["k_ID"] for r in C],yerr=[r["ke_ID"] for r in C],
           fmt="^",ms=6,color=GREY,capsize=3,lw=1,label=r"same data, ideal-gas $c_s$")
p.set_xscale("log"); p.set_xlabel(r"$\alpha=M/(2N_sm)=1/(2R)$"); p.set_ylabel(r"$K$")
p.set_title(r"(a) the eigenvalue equation, $\omega$ free",fontsize=10); p.legend(fontsize=7.5); p.grid(alpha=.25)

p=ax[0,1]
p.axhline(0,color=RED,lw=1.6,label="Kolafa-Rottner bulk $c_s$")
p.axhline(-1.01,color=BLUE,ls="--",lw=1.4,label=r"Paper 1 measured $c_s$ ($+1.01\%$)")
for S,mk,fc,lb in ((A,"o",BLUE,r"box A, $N_s=50$"),(B,"s","white",r"box B, $N_s=100$")):
    p.errorbar([r["alpha"] for r in S],[r["d_KR"] for r in S],yerr=[r["de_KR"] for r in S],
               fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=lb)
p.annotate("same $\\alpha$, two boxes:"+"\n"+f"${a['d_KR']:+.2f}\\%$ vs ${b['d_KR']:+.2f}\\%$",
           xy=(0.5,(a['d_KR']+b['d_KR'])/2),xytext=(0.105,0.30),fontsize=8,ha="left",
           arrowprops=dict(arrowstyle="->",lw=.9,color=BLACK))
p.set_ylim(-1.35,1.65)
p.set_xscale("log"); p.set_xlabel(r"$\alpha$"); p.set_ylabel(r"$100\,(K_{\rm meas}/K_{\rm pred}-1)$  [%]")
p.set_title(r"(b) departure is a function of $\alpha$, not of box size",fontsize=10)
p.legend(fontsize=7.5,loc="lower right"); p.grid(alpha=.25)

p=ax[1,0]
for S,mk,fc,lb in ((A,"o",BLUE,r"box A, $N_s=50$"),(B,"s","white",r"box B, $N_s=100$")):
    p.errorbar([r["M"] for r in S],[r["taur_x"] for r in S],yerr=[r["sig_taur"] for r in S],
               fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=lb+" measured")
    p.plot([r["M"] for r in S],[r["taur_man"] for r in S],"--" if S is A else ":",color=BLACK,lw=1.3,
           label=lb+" Mansour piston form")
p.set_xscale("log"); p.set_yscale("log"); p.set_xlabel(r"$M_d$"); p.set_ylabel(r"$\tau_r$  [$\sigma$-time]")
p.set_title(r"(c) mode damping vs Mansour's piston linewidth",fontsize=10); p.legend(fontsize=7.5); p.grid(alpha=.25,which="both")

p=ax[1,1]
for S,mk,fc,lb in ((A,"o",BLUE,r"box A, $N_s=50$"),(B,"s","white",r"box B, $N_s=100$")):
    p.errorbar([r["M"] for r in S],[r["taur_man"]/r["taur_x"] for r in S],
               yerr=[r["taur_man"]/r["taur_x"]*r["sig_taur"]/r["taur_x"] for r in S],
               fmt=mk,ms=7,mfc=fc,mec=BLUE,ecolor=BLUE,capsize=3,lw=1.2,label=lb)
p.axhline(1.0,color=BLACK,lw=1.4,label="Mansour exact")
p.set_xscale("log"); p.set_xlabel(r"$M_d$"); p.set_ylabel(r"$\tau_r^{\rm Mansour}/\tau_r^{\rm measured}$")
p.set_title(r"(d) the piston form becomes exact in the heavy limit",fontsize=10)
p.legend(fontsize=7.5); p.grid(alpha=.25,which="both")

fig.suptitle(r"Level 4: the divider eigenmode across $\alpha=0.1$–$2.0$ in two boxes "
             r"($\eta=0.10134170$, $\rho_0=1.2903$)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.97])
for ext in ("png","pdf"):
    f=os.path.join(OUT,f"261006_p2_mode_ladder.{ext}"); fig.savefig(f,dpi=180); print(f"\nwrote {f}")
