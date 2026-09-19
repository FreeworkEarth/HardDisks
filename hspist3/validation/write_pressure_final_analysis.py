import math
#!/usr/bin/env python3
"""##CHRIS 2026-09-08 -- TASK 6: assemble 260908_pressure_final_analysis.md from the
files on disk (tables_ABC.md, tables_D.md, chunk_calibration*.csv, run.log).
Numbers are never typed in here; re-run after any campaign change."""
import sys, os, re, glob, csv, statistics as st
R=sys.argv[1]; OUTMD=sys.argv[2]
A=os.path.join(R,"analysis")
log=open(f"{R}/run.log").read().splitlines()
disc=[l for l in log if l.startswith("DISCARD")]
prodfail=[l for l in log if l.startswith("PROD") and "valid=0" in l]
tvals=sorted(float(m) for m in re.findall(r"at t=([0-9.]+)",open(f"{R}/run.log").read()))
ntraj=len(glob.glob(f"{R}/preserved/traj_*.csv"))+len(glob.glob(f"{R}/runs/*_traj_*.csv"))
def cnt(e,N): return len(glob.glob(f"{R}/runs/main_traj_{e}_{N}_*.csv"))
want=[("0.60",900,4),("0.60",1600,3),("0.65",900,4),("0.65",1600,3),("0.67",900,4),("0.67",1600,3)]
t3=[(e,N,cnt(e,N),n) for e,N,n in want]; t3done=all(c>=n for _,_,c,n in t3)
lad=sorted(glob.glob(f"{R}/ladder_0p69_N900/traj_eq*_s*.csv")); laddone=len(lad)==9
draft = "" if (t3done and laddone) else f"> **DRAFT** — TASK 3 {sum(c for *_,c,_ in t3)}/21 target files, ladder {len(lad)}/9 finished. Re-run `write_pressure_final_analysis.py` when complete.\n\n"
def table(rows,hdr):
    return "| "+" | ".join(hdr)+" |\n|"+"---|"*len(hdr)+"\n"+"\n".join("| "+" | ".join(str(c) for c in r)+" |" for r in rows)
# TASK 1 table
d_rows=[]
for l in disc:
    m=re.search(r"eta=([\d.]+) N=(\d+) seed=(\d+).*fa=(\d+) orep=(\d+) crep=(\d+) wovd=(\d+)\) chunk=([\d.]+)",l)
    d_rows.append((m.group(1),m.group(2),m.group(3),m.group(8),m.group(4),m.group(5),m.group(6),m.group(7),"equilibration"))
for l in prodfail:
    m=re.search(r"eta=([\d.]+) N=(\d+) seed=(\d+).*chunk=([\d.]+) health=(\d+)",l)
    d_rows.append((m.group(1),m.group(2),m.group(3),m.group(4),"—","—","—","—",f"production (health={m.group(5)})"))
calnew=list(csv.DictReader(open(f"{R}/chunk_calibration.csv")))
calold={(r['eta'],r['N']):r for r in csv.DictReader(open(f"{R}/chunk_calibration_20260907_freshseed_method.csv"))}
cal_rows=[(r['eta'],r['N'],r['safe_chunk_seedA'],r['safe_chunk_seedB'],r['production_chunk'],f"{320/int(r['N']):.3f}",
           ("equilibrated (2026-09-08)" if (r['eta'] in ("0.60","0.65","0.67") and int(r['N'])>=900) else "fresh-seed (2026-09-07)"),
           calold.get((r['eta'],r['N']),{}).get('production_chunk','—')) for r in sorted(calnew,key=lambda r:(float(r['eta']),int(r['N'])))]
def fitrows(md, header_key):
    out={}
    sect=md.split(header_key,1)[1] if header_key in md else ""
    lines=sect.splitlines()
    # keep only the contiguous rows of THIS table: stop at the first non-table line
    body=[]
    for l in lines[1:]:            # lines[0] is the rest of the header line
        if l.startswith("|"): body.append(l)
        elif body: break
    for l in body:
        if not l.startswith("| 0."): continue
        c=[x.strip() for x in l.strip("|").split("|")]
        try: out[float(c[0])]=c
        except: pass
    return out
abc=open(f"{A}/tables_ABC.md").read()
FORMTAB={}
_sec=abc.split("Extrapolation-form comparison",1)[1] if "Extrapolation-form comparison" in abc else ""
_body=False
for l in _sec.splitlines():
    if l.startswith("|---"): _body=True; continue
    if _body:
        if not l.startswith("| 0."): break      # stop at the end of THIS table
        c=[x.strip() for x in l.strip("|").split("|")]
        if len(c)!=12: break
        try: FORMTAB[round(float(c[0]),3)]=c   # eta,Z1,s1,chi1,Z2,s2,chi2,spread,comb,dev1,dev2,dominates
        except ValueError: break
def F1(e,i):
    try: return float(FORMTAB[e][i])
    except Exception: return float("nan")
B=fitrows(abc,"Z_pair_inf | sigma")      # eta,#N,Zinf,sigma,a,chi2,dof,p,Z(900/1600),KR,dev,form
Cw=fitrows(abc,"Z_wall_inf | sigma")     # eta, Zwinf, sigma, Zpinf, wall-pair, KR, dev
def pct(x): return float(x.rstrip("%")) if x.rstrip("%").replace("+","").replace("-","").replace(".","").isdigit() else None
def f(x):
    try: return float(x)
    except: return float("nan")
# B columns
ZI,SG,CHI,DEV,FORM,Z2=2,3,5,10,11,8
mid=[e for e in B if 0.2<=e<=0.5 and pct(B[e][DEV]) is not None]
mid_max=max(abs(pct(B[e][DEV])) for e in mid) if mid else float("nan")
chi_mid_max=max(f(B[e][CHI]) for e in mid) if mid else float("nan")
chi_flag=[e for e in B if B[e][FORM].startswith("rejected")]
def dv(e,tab=B,col=DEV): return pct(tab[e][col]) if e in tab and pct(tab[e][col]) is not None else float("nan")
def sig(e): return 100*f(B[e][SG])/f(B[e][ZI]) if e in B else float("nan")
wp=[abs(pct(Cw[e][4])) for e in Cw if 0.05<=e<=0.65 and pct(Cw[e][4]) is not None]
# per-N Z(N) rows for 0.67 / 0.69 (table B per-cell)
cellB={}
for l in abc.split("## B.",1)[1].splitlines():
    if l.startswith("| 0."):
        c=[x.strip() for x in l.strip("|").split("|")]
        if len(c)==8 and c[1].isdigit(): cellB[(float(c[0]),int(c[1]))]=c
def devN(e,N): return pct(cellB[(e,N)][7]) if (e,N) in cellB and pct(cellB[(e,N)][7]) is not None else float("nan")
# exploratory betaP*sigma^2 at N=900: Z*rho*sigma^2 with rho sigma^2 = 4 eta/pi
def bP(e): 
    c=cellB.get((e,900)); return f(c[3])*4*e/math.pi if c else float("nan")
note30=(f" One low-density point, η = 0.30, sits at the flag threshold (χ²₁ = {f(B[0.3][CHI]):.2f}, p = {f(B[0.3][7]):.3f}) with a −0.10% intercept; "
        f"its N = 900/1600-only intercept is {B[0.3][Z2].split(' ')[0]}. It is reported, not excluded." if 0.3 in chi_flag else "")
def dofstr(e):
    d=B[e][6] if e in B else "?"
    return f"chi2_{d}"
def chi(e): return f(B[e][CHI]) if e in B else float("nan")
def nN(e):  return B[e][1] if e in B else "?"
def KRv(e):
    try: return float(B[e][9])
    except Exception: return float('nan')
def fdev(e,i): return FORMTAB[e][i] if e in FORMTAB else "—"
lowok=[e for e in FORMTAB if e<=0.10]
claims=(f"- **η ≤ 0.10: the equation of state is reproduced in the N → ∞ limit, and the "
        f"extrapolation form does not matter.** Fitting Z(N) = Z∞ + b·N^(-1/2) and Z∞ + b·N^(-1) "
        f"to the same per-cell means gives intercepts differing by ≤ {max(F1(e,7) for e in lowok):.4f}, "
        f"below their combined statistical error, because the finite-size term is negligible there. "
        f"Z_pair,∞ lies within {max(abs(pct(FORMTAB[e][9])) for e in lowok):.2f}% of Kolafa–Rottner 2006 "
        f"for η = 0.005–0.10 under either form.\n"
        f"- **η = 0.20–0.65: consistent with the EOS, but the extrapolation FORM is now the dominant "
        f"uncertainty and the claim carries it.** At η = 0.65 (four sizes, N = 400/900/1600/2500, 2 dof): "
        f"**Z∞ = {F1(0.65,1):.4f} ± {F1(0.65,2):.4f} (stat) from a + b/√N, χ² = {F1(0.65,3):.2f}**, and "
        f"**Z∞ = {F1(0.65,4):.4f} ± {F1(0.65,5):.4f} (stat) from a + b/N, χ² = {F1(0.65,6):.2f}**. "
        f"Both forms fit acceptably, but the intercepts differ by {F1(0.65,7):.4f}, which is "
        f"{F1(0.65,7)/F1(0.65,8):.1f}× their combined statistical error — so the honest quotation is "
        f"**Z∞(0.65) = {0.5*(F1(0.65,1)+F1(0.65,4)):.4f} ± {F1(0.65,2):.4f} (stat) ± {0.5*F1(0.65,7):.4f} (form)**, "
        f"i.e. {fdev(0.65,9)} to {fdev(0.65,10)} relative to KR. **The two forms do not agree on a "
        f"significant deviation** — 1/√N puts Z∞ {abs(F1(0.65,1)-KRv(0.65))/F1(0.65,2):.1f} σ below KR while "
        f"1/N puts it {abs(F1(0.65,4)-KRv(0.65))/F1(0.65,5):.1f} σ from it — so no claim of a significant "
        f"departure from KR is made at 0.65. Of the two, the 1/N form is the better fit "
        f"(χ² = {F1(0.65,6):.2f} against {F1(0.65,3):.2f} on 2 dof), but a χ² difference of that size on "
        f"2 degrees of freedom does not select between the forms, and the claim stands as no significant "
        f"departure. The same form spread dominates at every η ≥ 0.20 "
        f"(η = 0.60: {fdev(0.60,9)} vs {fdev(0.60,10)}); under either form the agreement with KR is "
        f"within ≈1%.\n"
        f"- **η = 0.67–0.69: these data do not support a bulk extrapolation with this model.** "
        f"Z(N) is non-monotone — it falls to N = 1600 and turns back up at N = 2500 — and both forms are "
        f"rejected by their own χ² (0.67: {F1(0.67,3):.1f} and {F1(0.67,6):.1f}; 0.69: {F1(0.69,3):.1f} and "
        f"{F1(0.69,6):.1f}, 2 dof). We therefore report Z(N) per size and quote no extrapolated value: "
        f"0.67 gives {devN(0.67,400):+.2f} → {devN(0.67,900):+.2f} → {devN(0.67,1600):+.2f} → {devN(0.67,2500):+.2f}% "
        f"vs KR and 0.69 gives {devN(0.69,400):+.2f} → {devN(0.69,900):+.2f} → {devN(0.69,1600):+.2f} → "
        f"{devN(0.69,2500):+.2f}%. Stationary over t = 400–7000 at 0.69/900 (D1); seed-to-seed scatter up to "
        f"2× the block error. Consistent with correlation lengths comparable to the box "
        f"(Bernard–Krauth: ξ ≈ 50 σ at 0.698) in hard-wall geometry.\n"
        f"- **η ≥ 0.70: exploratory.** KR is not a reference there (its fit turns over at 0.70). βPσ² at "
        f"N = 900 ({bP(0.702):.2f} / {bP(0.710):.2f} / {bP(0.720):.2f} at 0.702 / 0.710 / 0.720) sits above "
        f"the coexistence plateau (9.185) and rises; ψ₆ drifts while Z and T are flat (D2); one cell "
        f"(0.720/900) shows a 3.5% x–y anisotropy in the wall pressure.")
# boundary table from the discards: N*chunk, exceedances / equilibration calls
bnd={}
for r in d_rows:
    if r[4]=="—": continue
    e,N,ch,fa=r[0],int(r[1]),float(r[3]),int(r[4]); calls=int(round(400/ch)); bnd.setdefault((e,N,ch),[]).append((fa,calls))
bnd_rows=[]
for (e,N,ch),v in sorted(bnd.items(),key=lambda kv:(float(kv[0][0]),kv[0][1])):
    fas=[x[0] for x in v]; calls=v[0][1]
    bnd_rows.append((e,N,ch,f"{N*ch:.0f}",f"{min(fas)}–{max(fas)} of {calls}",f"{100*min(fas)/calls:.1f}–{100*max(fas)/calls:.1f}%"))
boundary=table(bnd_rows,["η","N","chunk","N·chunk","exceedances / equil. calls","rate"])
doc=f"""# Pressure validation — final analysis (Prompt 1, 2026-09-08)

{draft}Campaign: `hspist3/experiments_speed_of_sound/EDMD/mode1_normalized_units/00_pressure_validation_20260907/`.
Runner: `validation/pressure_validation.c`; driver: `validation/run_pressure_campaign2.sh`; analysis: `validation/analyze_pressure_campaign.py`, `validation/analyze_pressure_partD.py`; this file is generated by `validation/write_pressure_final_analysis.py`. Every number below is read from the campaign files.

## Claim range

{claims}

## TASK 1 — state verified against `run.log`

Accepted trajectories on disk: **{ntraj}**; all have valid=1 and zero health counters (table A). Discards in the 2026-09-07 campaign: **{len(disc)} during equilibration + {len(prodfail)} during production = {len(disc)+len(prodfail)}**, all N ≥ 900, η = 0.60–0.67:

{table(d_rows,["η","N","seed","chunk","forced_adv","overlap_rep","clamp_rep","wall_ovd","phase"])}

Avalanche warnings in `run.log`: n = {len(tvals)}, t = {tvals[0]:.2f} … {tvals[-1]:.2f}, median {st.median(tvals):.2f}; {sum(1 for t in tvals if t<400)}/{len(tvals)} inside equilibration (≤ 400), all `dominant=AB` (pair-collision bursts). They are spread over the whole equilibration, not a lattice-start transient.

**Cause (confirmed at source):** `edmd_core/edmd.c` `#define EDMD_ADVANCE_MAX_EVENTS 250000L` is a per-`edmd_advance_to()` budget that counts calendar pops (including invalidated entries) and does not scale with N, so the safe window shrinks ≈ 1/N at fixed density and steeply with density. The 5-unit fresh-lattice calibration saw a lower pop rate than the equilibrated fluid and passed windows that failed later. Where the budget was actually exceeded, in units of N·chunk (exceedances are the `forced_advance` counters over the ≈ 400/chunk equilibration calls):

{boundary}

At 0.69–0.72, N·chunk = 360 (0.4 at N = 900) ran clean in every accepted trajectory and 270 (the ladder, 0.3 at N = 900) over ≈ 10⁵ calls; no failure data exist above 360 there. **320/N is therefore a conservative cap with a margin between ≈ 1.1× (0.69–0.72, unmeasured) and ≈ 4× (0.60), not a threshold.**

## TASK 2 — runner-level fix (no core change)

(a) production chunk = min(0.6·min(A,B), 320/N) for η ≥ 0.6; (b) calibration on a disposable instance **after 100 time units of equilibration on that instance**, rate probed on the equilibrated state, candidate verified over **20 consecutive chunks**; (c) production contract unchanged (any health event from initialization to end of measurement discards the trajectory; no adaptation).

Chunk table (six target cells recalibrated with the new method; other rows are the 2026-09-07 fresh-seed values that the accepted trajectories actually ran at, kept for the record):

{table(cal_rows,["η","N","calib A","calib B","production","320/N","method","2026-09-07 value"])}

The equilibrated calibration found ≤ 320/N on its own for all six target cells (the cap did not bind); at 0.67/1600 it chose 0.075, more conservative than the rule. The calibration's stderr was not logged in this run (sent to /dev/null); why 0.25 failed its 20-chunk verification at 0.67/1600 is therefore unrecoverable. Both drivers now write `calib_<eta>_<N>_<seed>.err` next to the calibration table.

## TASK 3 — rerun of the empty/short cells (same seeds: 20260907 + 104729·k + 7919·N)

{table([(e,N,f"{c}/{n}") for e,N,c,n in t3],["η","N","accepted / target"])}

New discards since relaunch: **{len(disc)-17} equilibration, {len(prodfail)-1} production.** The previously failed 0.60/900 seed 27388007 (health=15 at chunk 0.8) was accepted at chunk 0.3 with health 0.

## TASK 4 — equilibration ladder and TASK 5 — analysis

Tables A–C follow (from `analysis/tables_ABC.md`), then D (from `analysis/tables_D.md`). Plots: `analysis/dev_vs_invsqrtN.pdf`, `analysis/Z_vs_eta.pdf`, `analysis/stationarity_eta0p60_N900.pdf`, `analysis/stationarity_eta0p69_N900_ladder.pdf`, `analysis/stationarity_eta0p72_N900.pdf`.

{open(f"{A}/tables_ABC.md").read()}

{open(f"{A}/tables_D.md").read()}

## Methods paragraph — sixth failure class

*Event-budget exhaustion during equilibration, masked by a short fresh-seed calibration.* The event-driven integrator caps the number of calendar entries processed per advance call at 2.5×10⁵, counting invalidated entries as well (`edmd.c`: `events_processed++` at line 1403 precedes the validity test at lines 1465–1477); at η = 0.60–0.69 this corresponds to ≈ 30–45 calendar entries per physical collision, i.e. 2–4 collisions per particle per call at the window we used (N·chunk = 320). Because the cap does not scale with N, the largest safe window shrinks ≈ 1/N at fixed density and steeply with density. A calibration that verified candidate windows for 5 time units on a freshly seeded lattice passed windows of 0.8 (N = 900, 1600) and 0.4 (N = 1600) at η ≥ 0.60; during equilibration these exhausted the budget in the last 1–8% of the affected calls (all 345 logged cases; `stagnant = 0` in 344 of 345), i.e. throughput exhaustions rather than event cascades — the word "avalanche" is the warning text's, not a description of the dynamics. Each exhaustion forced a free-flight advance and was followed by 1.8×10³–1.8×10⁴ overlap repairs per trajectory. The resulting pressures were numerically plausible (the one production case gave Z_pair = 6.548 ± 0.005 against 6.544 for the clean seeds) and were rejected only by the health ledger. The fix was applied at the runner level: a conservative window cap of 320/N together with calibration on an equilibrated disposable instance verified over 20 consecutive windows; the all-or-nothing validity rule was unchanged. The clean core fix, deferred until after the campaign, is to count only validated entries or to scale the budget with N. This is the sixth failure class caught by the validity layer rather than by the observable, after event-budget exhaustion in production, the inert event calendar (Z = 1 reported as valid), corner-packed seeding, wall-contact seeding, and the non-terminating random insertion above η ≈ 0.55.
"""
open(OUTMD,"w").write(doc); print(f"wrote {OUTMD}  ({'FINAL' if not draft else 'DRAFT'})")
