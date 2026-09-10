#!/usr/bin/env python3
"""##CHRIS 2026-09-08 -- TASK 5 part D: structural boundary and stationarity.
Reads the block CSVs (headerless, 12 cols) from runs/ and the eta=0.69 N=900
equilibration ladder. Produces: ladder table (Z_pair per seed vs equilibration
time, with the equilibration trace), per-trajectory drift tests over the 30
measurement blocks, and stationarity panels for eta = 0.60, 0.69, 0.72.
Usage: analyze_pressure_partD.py CAMPAIGN_DIR -> CAMPAIGN_DIR/analysis/"""
import sys, os, glob, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=sys.argv[1]; OUT=os.path.join(R,"analysis"); os.makedirs(OUT,exist_ok=True)
C="eta N seed idx t0 t1 T Zp Zx Zy p6g p6l".split()
def read(f):
    rows=[dict(zip(C,[float(v) for v in p])) for p in (l.strip().split(",") for l in open(f)) if len(p)==12 and p[0]!="eta"]
    return sorted(rows,key=lambda r:r["idx"])
def drift(vals):
    """slope per block and its t-value, plus first/second-half means"""
    y=np.array(vals); x=np.arange(len(y),dtype=float)
    if len(y)<4: return float("nan"),float("nan"),float("nan"),float("nan")
    A=np.vstack([np.ones_like(x),x]).T; coef,res,_,_=np.linalg.lstsq(A,y,rcond=None)
    resid=y-A@coef; s2=resid@resid/max(1,len(y)-2); cov=s2*np.linalg.inv(A.T@A)
    return coef[1], coef[1]/math.sqrt(cov[1,1]) if cov[1,1]>0 else float("nan"), y[:len(y)//2].mean(), y[len(y)//2:].mean()
md=["## D. Structural boundary: equilibration ladder and stationarity\n"]

# ---- ladder ----
L=f"{R}/ladder_0p69_N900"
md.append("### D1. Ladder at eta = 0.69, N = 900, chunk 0.3: one trajectory per seed, three windows (seeds = production seeds k=0,1,2)\n")
md.append("| equil | seed | equil blocks | Z_pair (meas, 30 blk) | sem | Z_wall | T | psi6_g | psi6_l | Z_pair last-5 equil blk | health |\n|---|---|---|---|---|---|---|---|---|---|---|")
ladder={}
for f in sorted(glob.glob(f"{L}/blk_eq*_s*.csv")):
    b=os.path.basename(f); eq=int(b.split("eq")[1].split("_")[0]); k=int(b.split("_s")[1].split(".")[0])
    rows=read(f); eqb=[r for r in rows if r["idx"]<0]; mb=[r for r in rows if r["idx"]>=0]
    if len(mb)<30: continue   # not finished
    tj=f"{L}/traj_eq{eq}_s{k}.csv"
    h="—"
    if os.path.exists(tj):
        t=open(tj).read().strip().split(","); h=str(int(float(t[22]))+int(float(t[23]))+int(float(t[24]))+int(float(t[25])))
    Zp=np.array([r["Zp"] for r in mb]); Zw=np.array([0.5*(r["Zx"]+r["Zy"]) for r in mb])
    last5=np.mean([r["Zp"] for r in eqb[-5:]]) if len(eqb)>=5 else float("nan")
    ladder[(eq,k)]=(Zp.mean(),Zp.std(ddof=1)/math.sqrt(len(Zp)))
    md.append(f"| {eq} | k={k} | {len(eqb)} | {Zp.mean():.4f} | {Zp.std(ddof=1)/math.sqrt(len(Zp)):.4f} | {Zw.mean():.4f} | {np.mean([r['T'] for r in mb]):.3f} | "
              f"{np.mean([r['p6g'] for r in mb]):.3f} | {np.mean([r['p6l'] for r in mb]):.3f} | {last5:.4f} | {h} |")
md.append("\nZ_pair per seed and window (window = [t_eq, t_eq+600]; blank = not finished):\n\n| seed | window 400–1000 | window 1600–2200 | window 6400–7000 |\n|---|---|---|---|")
for k in (0,1,2):
    md.append(f"| k={k} | "+" | ".join(f"{ladder[(eq,k)][0]:.4f} ± {ladder[(eq,k)][1]:.4f}" if (eq,k) in ladder else "" for eq in (400,1600,6400))+" |")
# ---- block-identity check: are the three "rungs" one trajectory per seed? ----
# The rungs share seed and chunk and the runner is deterministic, so if the
# equilibration/measurement blocks at equal t_start are bit-identical, the ladder
# is one 7000-unit trajectory per seed read at three windows -- not three
# differently equilibrated systems.
def blocks_by_t(f):
    out={}
    for l in open(f):
        p=l.strip().split(",")
        if len(p)!=12 or p[0]=="eta": continue
        out[p[4]]=tuple(p[6:12])          # key t_start -> (T,Zp,Zx,Zy,p6g,p6l) as written
    return out
ident_rows=[]; all_identical=True
for k in (0,1,2):
    ref=f"{L}/blk_eq6400_s{k}.csv"
    if not os.path.exists(ref): continue
    B6=blocks_by_t(ref)
    for eq in (400,1600):
        f=f"{L}/blk_eq{eq}_s{k}.csv"
        if not os.path.exists(f): continue
        Bq=blocks_by_t(f); common=[t for t in Bq if t in B6]
        same=sum(1 for t in common if Bq[t]==B6[t]); diff=len(common)-same
        if diff: all_identical=False
        ident_rows.append(f"| k={k} | eq{eq} vs eq6400 | {len(Bq)} | {len(common)} | {same} | {diff} |")
md.append("\nBlock-identity check (blocks compared at equal t_start; a block is 'identical' if T, Z_pair, Z_wall_x, Z_wall_y, psi6_global, psi6_local match to the written precision):\n")
md.append("| seed | comparison | blocks in shorter run | blocks at common t_start | identical | different |\n|---|---|---|---|---|---|")
md.extend(ident_rows)
if all_identical and ident_rows:
    md.append("\n**Reading:** every block of the eq = 400 and eq = 1600 runs is bit-identical to the block at the same t_start in the eq = 6400 run of the same seed. "
              "The three rungs are therefore **one trajectory per seed, read at three windows** (t ∈ [400,1000], [1600,2200], [6400,7000] of a 7000-unit run), "
              "not three differently equilibrated systems. What the ladder shows is that at eta = 0.69, N = 900, Z_pair is **stationary over t = 400–7000 in three independent trajectories** "
              "(seed-to-seed sd ≈ 0.02; the k = 1 trajectory moves up by ≈ 0.065 ≈ 3 block-σ between its first and second window and stays there). "
              "It does not show that a differently equilibrated route to 0.69 lands on the same Z; that test was skipped (below).\n")
else:
    md.append("\n**Reading:** blocks differ between rungs; the rungs are not a single trajectory. Treat the rows above as independent equilibrations.\n")
md.append("\nSaved-configuration restart (equilibrated eta=0.65 compressed to 0.69): **skipped** -- the EDMD core has no state save/load; adding one is a core change, excluded during a campaign.\n")

# ---- drift tests on production cells ----
md.append("\n### D2. Stationarity of the 30 measurement blocks (dense cells): linear drift of Z_pair and both psi6\n")
md.append("| eta | N | seed | Z_pair slope/blk | t | first-half | second-half | psi6_g slope t | psi6_l slope t |\n|---|---|---|---|---|---|---|---|---|")
byc={}
for f in sorted(glob.glob(f"{R}/runs/*_blk_*.csv")):
    rows=[r for r in read(f) if r["idx"]>=0]
    if len(rows)<30: continue
    e=rows[0]["eta"]; N=int(rows[0]["N"]); s=int(rows[0]["seed"])
    if e<0.6: continue
    sZ,tZ,h1,h2=drift([r["Zp"] for r in rows]); _,tg,_,_=drift([r["p6g"] for r in rows]); _,tl,_,_=drift([r["p6l"] for r in rows])
    byc.setdefault((e,N),[]).append((s,rows))
    md.append(f"| {e:.3f} | {N} | {s} | {sZ:+.5f} | {tZ:+.2f} | {h1:.4f} | {h2:.4f} | {tg:+.2f} | {tl:+.2f} |")
md.append("\n|t| > ~2.5 would indicate a drift over the window; the table is what it is.\n")

# ---- stationarity panels ----
def panel(cells, title, fname, with_ladder=None):
    fig,axs=plt.subplots(6,1,figsize=(10,12),sharex=True)
    lab=["T","Z_pair","Z_wall_x","Z_wall_y","|<psi6>|","<|psi6_i|>"]; keys=["T","Zp","Zx","Zy","p6g","p6l"]
    for (s,rows) in cells:
        x=[r["idx"] for r in rows]
        for a,k in zip(axs,keys): a.plot(x,[r[k] for r in rows],"-",lw=0.9,alpha=0.85,label=f"seed {s}")
    if with_ladder:
        for f in sorted(glob.glob(with_ladder)):
            rows=read(f); x=[r["idx"] for r in rows]
            for a,k in zip(axs,keys): a.plot(x,[r[k] for r in rows],"-",lw=0.8,alpha=0.6)
    for a,l in zip(axs,lab): a.set_ylabel(l,fontsize=8); a.grid(True,ls=":",alpha=0.5); a.axvline(0,color="k",lw=0.8,ls="--")
    axs[0].set_title(title,fontsize=10); axs[0].legend(fontsize=6,ncol=4); axs[-1].set_xlabel("block index  (negative = equilibration, block_dt = 20)")
    fig.tight_layout(); fig.savefig(f"{OUT}/{fname}.pdf"); fig.savefig(f"{OUT}/{fname}.png",dpi=150); plt.close(fig)
for e,N,fn,lad in ((0.6,900,"stationarity_eta0p60_N900",None),
                   (0.69,900,"stationarity_eta0p69_N900_ladder",f"{L}/blk_eq6400_s*.csv"),
                   (0.72,900,"stationarity_eta0p72_N900",None)):
    cells=byc.get((e,N),[])
    if cells or lad: panel(cells,f"eta={e} N={N}: measurement blocks"+(" + eq=6400 ladder (equilibration shown at negative index)" if lad else ""),fn,lad)
open(f"{OUT}/tables_D.md","w").write("\n".join(md)+"\n"); print("wrote tables_D.md + stationarity panels")
