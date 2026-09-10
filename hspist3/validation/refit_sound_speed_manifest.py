#!/usr/bin/env python3
"""##CHRIS: TASK C -- trajectory-level fit-input manifest and refit of c_s(eta)
for the route-A campaign, on ELIGIBLE trajectories only.

Joins, per (eta, M, run):
  runs.csv (nu, peak flags, quality)  x  leaf run.log (run -> seed, health counters)
  x  speed_of_sound_psi6.csv (seed)   x  GPT structural flags (seed).
Eligibility (all must hold):
  frequency_quality_pass == 1, peak_on_search_boundary == 0,
  forced_advance == overlap_repair == clamp_repair == 0,
  wall_overdue == 0 OR wall_overdue is the documented t=0 seeder class
  (260909_wall_overdue_and_temperature_resolution_CC.md: all 662 events at
  t=0, RIGHT wall, gap=0; none during hold or measurement).
Fit: Roman relation, K the fundamental root of cot K = alpha K (bisection,
alpha = M/(2 N_side)), x = K/(2 pi L_eff), L_eff = L0 - 2r, c_s = slope of
nu vs x through the origin (weights 1/sem^2 of the per-mass mean nu).
EOS columns: Kolafa-Rottner 2006 (quoted only for eta <= 0.69) and Liu 2021
global, both mapped through the adiabatic c_s^2 = (kT/m)[Z + eta Z' + Z^2].
Originals are never modified; unmatched joins are recorded, not guessed.
"""
import csv, glob, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

CAMP = sys.argv[1]; GPT = sys.argv[2]; OUT = sys.argv[3]
os.makedirs(OUT, exist_ok=True)
R, N_SIDE, KR_MAX = 0.5, 50, 0.69

gpt = {}
for r in csv.DictReader(open(f"{GPT}/trajectories.csv")):
    if "campaign_r25_psi6_20260823" in r["source"]:
        gpt[int(r["seed"])] = r

man, summ, unmatched = [], [], []
for leaf in sorted(glob.glob(f"{CAMP}/eta_*")):
    if not os.path.isdir(leaf) or not os.path.exists(f"{leaf}/run.log"): continue
    tag = os.path.basename(leaf)
    # match the analysis directory by eta (the analysis names it from the
    # computed eta to 6 decimals; dense leaves are named to 3), tolerance 5e-4
    pe = os.path.join(leaf, "speed_of_sound_psi6.csv")
    leaf_eta = float(next(csv.DictReader(open(pe)))["eta"]) if os.path.exists(pe) else float(tag.replace("eta_","").replace("p","."))
    ana = [d for d in glob.glob(f"{CAMP}/analysis/analysis/eta_*_L0_*/speed_of_sound_runs.csv")
           if abs(float(re.search(r"eta_([\dp]+)_L0", d).group(1).replace("p",".")) - leaf_eta) < 5e-4]
    if len(ana) != 1: unmatched.append((tag, "runs.csv", f"{len(ana)} candidates for eta={leaf_eta}")); continue
    runs = list(csv.DictReader(open(ana[0])))
    seed_of, health = {}, {}
    for ln in open(f"{leaf}/run.log"):
        m = re.search(r"Running: L0 = ([\d.]+), M = (\d+)\*m, run = (\d+), seed = (\d+)", ln)
        if m: seed_of[(int(m.group(2)), int(m.group(3)))] = int(m.group(4))
        h = re.search(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): forced_advance=(\d+) wall_clamp_repairs=(\d+) overlap_repairs=(\d+) wall_overdue=(\d+)", ln)
        if h: health[(int(h.group(1)), int(h.group(2)))] = tuple(int(h.group(i)) for i in (4,5,6,7))
    psi = {int(r["seed"]): r for r in csv.DictReader(open(f"{leaf}/speed_of_sound_psi6.csv"))} if os.path.exists(f"{leaf}/speed_of_sound_psi6.csv") else {}
    L0 = float(runs[0]["L0"]); eta = N_SIDE*2*math.pi*R*R/(2*L0*10.0)
    per_mass = {}
    for r in runs:
        M, run = int(float(r["wall_mass_factor"])), int(float(r["run"]))
        seed = seed_of.get((M, run))
        fa, cr, orp, wo = health.get((M, run), (0,0,0,0))
        g = gpt.get(seed) if seed is not None else None
        p = psi.get(seed) if seed is not None else None
        if seed is None: unmatched.append((tag, f"seed for M={M} run={run}", "run.log"))
        if g is None and seed is not None: unmatched.append((tag, f"GPT flags seed={seed}", "trajectories.csv"))
        qual = r["frequency_quality_pass"] in ("1","True","true")
        bnd  = r["peak_on_search_boundary"] in ("1","True","true")
        wo_class = "none" if wo == 0 else ("t0_seeder_right_wall" if (fa==0 and cr==0 and orp==0) else "unclassified")
        reasons = []
        if not qual: reasons.append("quality_fail:"+r["frequency_quality_reasons"])
        if bnd: reasons.append("peak_on_search_boundary")
        if fa or cr or orp: reasons.append(f"health fa={fa} crep={cr} orep={orp}")
        if wo_class == "unclassified": reasons.append(f"wall_overdue={wo} with other health")
        elig = not reasons
        nu = float(r["nu"]) if r["nu"] not in ("", "nan") else float("nan")
        man.append(dict(eta=f"{eta:.6f}", L0=L0, M=M, run=run, seed=seed, nu=r["nu"], peak_snr=r["peak_snr"],
                        peak_on_search_boundary=int(bnd), frequency_quality_pass=int(qual),
                        forced_advance=fa, clamp_repair=cr, overlap_repair=orp, wall_overdue=wo, wall_overdue_class=wo_class,
                        psi6_global_hold=(p or {}).get("psi6_global_hold",""), psi6_global_end=(p or {}).get("psi6_global_end",""),
                        gpt_health_reason=(g or {}).get("health_reason",""), gpt_eligible_structure=(g or {}).get("eligible_structure",""),
                        eligible=int(elig), exclusion_reasons=";".join(reasons)))
        if elig and nu == nu: per_mass.setdefault(M, []).append(nu)
    # ---- refit ----
    xs, ys, ws, n_used = [], [], [], 0
    L_eff = L0 - 2*R
    for M, nus in sorted(per_mass.items()):
        if len(nus) < 3: continue
        K = sos.k_root_bisect(M/(2.0*N_SIDE))
        xs.append(K/(2*math.pi*L_eff)); ys.append(float(np.mean(nus)))
        ws.append(max(float(np.std(nus, ddof=1))/math.sqrt(len(nus)), 1e-12)); n_used += len(nus)
    if len(xs) < 3: unmatched.append((tag, "fit", f"only {len(xs)} masses eligible")); continue
    x, y, s = np.array(xs), np.array(ys), np.array(ws)
    cs, cs_err = sos.weighted_linreg(x, y, s, force_zero_intercept=True)[:2]
    cf, cf_err, b0, b0_err = sos.weighted_linreg(x, y, s, force_zero_intercept=False)[:4]
    resid = y - cs*x; r2 = 1 - float(((resid)**2).sum()/((y-y.mean())**2).sum())
    a = np.array([eta]); Zk = sos.Z_kolafa_rottner_2006(a); Zl = sos.Z_liu_global(a); dZl = sos.dZ_liu_global(a)
    h = 1e-5; dZk = (sos.Z_kolafa_rottner_2006(a+h)-sos.Z_kolafa_rottner_2006(a-h))/(2*h)
    cs_kr = float(sos.cs_adiabatic_2d_monatomic(Zk, dZk, a, kbt=1, m=1)[0]) if eta <= KR_MAX else float("nan")
    cs_liu = float(sos.cs_adiabatic_2d_monatomic(Zl, dZl, a, kbt=1, m=1)[0])
    n_all = len(runs); n_el = sum(1 for m_ in man if m_["eta"] == f"{eta:.6f}" and m_["eligible"])
    summ.append(dict(eta=f"{eta:.6f}", L0=L0, L_eff=L_eff, c_s=f"{cs:.5f}", c_s_err=f"{cs_err:.5f}",
                     c_s_free_intercept=f"{cf:.5f}", intercept=f"{b0:.5g}", fit_r2=f"{r2:.5f}",
                     n_masses=len(xs), n_runs_all=n_all, n_runs_eligible=n_el, n_runs_in_fit=n_used,
                     cs_KR2006=f"{cs_kr:.5f}" if cs_kr==cs_kr else "", dev_KR=f"{100*(cs-cs_kr)/cs_kr:+.2f}%" if cs_kr==cs_kr else "",
                     cs_Liu2021=f"{cs_liu:.5f}", dev_Liu=f"{100*(cs-cs_liu)/cs_liu:+.2f}%"))

with open(f"{OUT}/fit_input_manifest.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(man[0].keys())); w.writeheader(); w.writerows(man)
with open(f"{OUT}/refit_cs_vs_eta.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
with open(f"{OUT}/unmatched_joins.csv","w",newline="") as fh:
    w = csv.writer(fh); w.writerow(["leaf","item","detail"]); w.writerows(unmatched)
print(f"manifest rows {len(man)}  eligible {sum(m['eligible'] for m in man)}  refit eta {len(summ)}  unmatched {len(unmatched)}")
print("exclusion reasons:", {k: sum(1 for m in man if k in m['exclusion_reasons']) for k in ("quality_fail","peak_on_search_boundary","health fa","wall_overdue")})
print("wall_overdue classes:", {c: sum(1 for m in man if m['wall_overdue_class']==c) for c in ("none","t0_seeder_right_wall","unclassified")})
