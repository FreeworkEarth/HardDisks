#!/usr/bin/env python3
"""##CHRIS 2026-10-08: the demonstration run. Applies the pre-registered rules of
261007_paper2_level4_Rcollapse2_demo.md section 3 exactly as written (committed 1b9b517, before any
record existed). Analysis only."""
import glob,csv,json,math,os,sys
import numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); REPO=os.path.dirname(os.path.dirname(HERE))
P=os.path.join(REPO,"hspist3","experiments_energy_transfer","level4_demo_20261007")
NS, CS, LC, TAUR, TAUT, TAUGP = 50, 1.749302, 38.75, 2090., 40079., 10070.
W_KR, W_ID, DT_KR, DT_ID = 7.0715, 5.5556, 0.14143, 0.11111
X_FINAL_SHIFT = -1.9375                      # divider displacement: half the 3.875 piston travel

def cell(tag):
    S=[]; 
    for f in sorted(glob.glob(os.path.join(P,tag,"summary_*.csv"))):
        S.append(list(csv.DictReader(open(f)))[0])
    D=[]
    for f in sorted(glob.glob(os.path.join(P,tag,"red_*.csv"))):
        D.append(pd.read_csv(f))
    return S,D

def arr(S,k): return np.array([float(r[k]) for r in S])
def sem(a): return a.std(ddof=1)/math.sqrt(len(a))

out={}
for tag in ("B1long","B2pilot"):
    S,D=cell(tag)
    if not D: print(f"  {tag}: no data"); continue
    u=arr(S,"piston_right_v_step_constant").mean(); W=arr(S,"W_in_max"); stop=arr(S,"piston_stop_t_rel")
    n=min(len(d) for d in D)
    t=D[0]["Time"].values[:n]
    dT=np.array([(d["KE_gas_right"].values[:n]-d["KE_gas_left"].values[:n])/NS for d in D])
    x =np.array([ d["W0_x_sigma"].values[:n]-d["W0_x_sigma"].values[0] for d in D])
    m=len(D); noise=0.199/math.sqrt(m)
    mech=(t>10*TAUR)                                       # mechanical equilibrium reached
    dTm=dT[:,mech].mean(axis=1); xm=x[:,mech].mean(axis=1)
    Dfrac=dTm.mean()*NS/W.mean()
    out[tag]=dict(u=float(u),n=m,W=float(W.mean()),Wsem=float(sem(W)),stop=float(stop.mean()),
                  dT_mech=float(dTm.mean()),dT_mech_sem=float(sem(dTm)),noise=float(noise),
                  x_mech=float(xm.mean()),x_mech_sem=float(sem(xm)),D=float(Dfrac),
                  D_sem=float(sem(dTm)*NS/W.mean()),Ma=float(u/CS),record=float(t[-1]))
    r=out[tag]
    print(f"\n=== {tag}: u = {u:.2f}, Ma = {u/CS:.4f}, {m} seeds, record {t[-1]:.0f} ===")
    print(f"  W_in            = {W.mean():.4f} ± {sem(W):.4f}   (KR isentropic {W_KR}, ideal {W_ID})")
    print(f"  push ends at      {stop.mean():.3f} sigma   (transit L_c/c_s = {LC/CS:.2f})")
    print(f"  T1-T2 at mech eq  {dTm.mean():+.5f} ± {sem(dTm):.5f}   (seed noise floor {noise:.4f})")
    print(f"  **D = (T1-T2) N/W_in = {Dfrac:+.4f} ± {sem(dTm)*NS/W.mean():.4f}**   predicted Ma^2 = {(u/CS)**2:.4f}")
    print(f"  divider dx        {xm.mean():+.4f} ± {sem(xm):.4f}   (predicted {X_FINAL_SHIFT})")
    # decay of |dT| across the record
    mean=dT.mean(axis=0); se=dT.std(axis=0,ddof=1)/math.sqrt(m)
    for frac in (0.0,0.5,1.0,2.0,5.0,10.0,20.0):
        i=int(np.searchsorted(t,stop.mean()+frac*TAUR)); i=min(i,n-1)
        print(f"     t = push+{frac:4.1f} tau_r = {t[i]:8.0f}:  T1-T2 = {mean[i]:+.4f} ± {se[i]:.4f}")

print("\n\n### The pre-registered rules, applied")
b=out.get("B1long")
if b:
    print(f"(i)  dKE_right/W_in and dKE_left/W_in at push end: reported separately from B1zoom.")
    print(f"(ii) T1-T2 at mechanical equilibrium = {b['dT_mech']:+.5f} ± {b['dT_mech_sem']:.5f}, "
          f"noise floor {b['noise']:.4f} -> {'BELOW' if abs(b['dT_mech'])<b['noise'] else 'ABOVE'} the floor.")
    print(f"     Fit rule: window is where |T1-T2| > 3x noise = {3*b['noise']:.4f}.")
    print(f"(iii) divider dx = {b['x_mech']:+.4f} ± {b['x_mech_sem']:.4f} vs predicted {X_FINAL_SHIFT} "
          f"({abs(b['x_mech']-X_FINAL_SHIFT)/b['x_mech_sem']:.1f} sigma)")
if "B2pilot" in out and b:
    p2=out["B2pilot"]
    print(f"(ii') ORDERING: D(u=1.0) = {p2['D']:+.4f} ± {p2['D_sem']:.4f} vs D(u=0.2) = {b['D']:+.4f} ± {b['D_sem']:.4f}"
          f"  -> {'HOLDS' if p2['D']>b['D'] else 'FAILS'}")
json.dump(out,open(os.path.join(HERE,"261007_demo_results.json"),"w"),indent=1)
print(f"\nwrote 261007_demo_results.json")
