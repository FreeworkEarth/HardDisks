#!/usr/bin/env python3
"""Independent verification of D57-3 task 2 (the bistable run). Read-only on the repo; writes only into this script's folder.
Run with python3 -B (no __pycache__ written into the repo's helper folders)."""
import sys, os, math, glob, json
sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from scipy.optimize import brentq, curve_fit

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
for p in (os.path.join(HS, "validation"), HS):
    sys.path.insert(0, p)
import resched_testT_followup_261007 as FU
import tests_20260913 as T

OUTD = os.path.dirname(os.path.abspath(__file__))
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016", "P2", "merged", "stage2", "gen3")
CELL, M, RUN = "P2_N0900_eta0.780000_G2_pin0", 1800, 9
D = os.path.join(P2, CELL, f"m_{M}")
CCB = os.path.join(HS, "experiments_gen3_d57_261016", "bistable")
res = {}


def kroot_own(alpha):
    # root K in (0, pi/2) of cot K = alpha K
    return brentq(lambda k: math.cos(k) / math.sin(k) - alpha * k, 1e-9, math.pi / 2 - 1e-12, xtol=1e-15, rtol=1e-15, maxiter=500)


def g_model(t, A, tT, B, tr, om, ph):
    return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(om * t + ph)


def fit_acf_own(a, dt, nu_ref):
    # my own transcription of the item-4a model, bounds and start (to cross-check FU.fit_acf)
    c = a / a[0]; lag = np.arange(len(c)) * dt; om0 = 2 * math.pi * nu_ref; P = 1 / nu_ref
    p, _ = curve_fit(g_model, lag, c, p0=[0.3, 2 * P, 0.7, 5 * P, om0, 0.0],
                     bounds=([0, 2 * dt, 0, 2 * dt, 0.2 * om0, -math.pi / 2], [1.5, 1e7, 1.5, 1e7, 5 * om0, math.pi / 2]), maxfev=200000)
    return p


def seed_avg_fit(acfs, dt, nu_ref, fitter=FU.fit_acf):
    nl = min(len(a) for a in acfs)
    m = np.mean(np.vstack([np.asarray(a[:nl], float) for a in acfs]), axis=0)
    p = fitter(m, dt, nu_ref)
    rms = float(np.sqrt(np.mean((m / m[0] - g_model(np.arange(nl) * dt, *p)) ** 2)))
    return p, nl, rms


def e1_mass(d, exclude=(), nu_ref_over="used"):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run")
    nud = pd.read_csv(os.path.join(d, "nud_runs.csv")); hl = pd.read_csv(os.path.join(d, "runs_health.csv"))
    z = np.load(os.path.join(d, "acf_runs.npz"))
    clean = set(hl.loc[(hl.gen3_clean == 1) & (hl.gen3_records == 1), "run"].astype(int))
    ok = set(nud.loc[nud.fit_ok == 1, "run"].astype(int))
    used = sorted((set(red.run.astype(int)) & clean & ok) - set(exclude))
    dt = float(red.dt.iloc[0])
    if nu_ref_over == "used":
        nu_ref = float(red.set_index("run").loc[used, "nu"].mean())
    else:  # all folder runs
        nu_ref = float(red.nu.mean())
    nd = nud.set_index("run").loc[used, "nu_d"].to_numpy(float)
    acfs = [np.asarray(z[f"run{r}"], float) for r in used]
    p, nl, rms = seed_avg_fit(acfs, dt, nu_ref)
    nu0 = math.sqrt(nd.mean() ** 2 + (1 / (2 * math.pi * p[3])) ** 2)
    jk, jkp = [], []
    for i in range(len(used)):
        keep = [k for k in range(len(used)) if k != i]
        pi_, _, _ = seed_avg_fit([acfs[k] for k in keep], dt, nu_ref)
        jkp.append(pi_)
        jk.append(math.sqrt(nd[keep].mean() ** 2 + (1 / (2 * math.pi * pi_[3])) ** 2))
    jk = np.array(jk); n = len(jk)
    se = math.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2))
    # stall check: start values A=0.3, tau_T=2/nu_ref, B=0.7, tau_r=5/nu_ref
    def stalled(q):
        return bool(abs(q[0] - 0.3) < 1e-3 and abs(q[2] - 0.7) < 1e-3)
    return dict(used=len(used), used_runs=used, nu_ref=nu_ref, nl=nl, nud=float(nd.mean()), A=float(p[0]), tau_T=float(p[1]), B=float(p[2]),
                tau_r=float(p[3]), omega=float(p[4]), phi=float(p[5]), rms=rms, tau_r_start=5 / nu_ref, stalled=stalled(p),
                jk_stalls=int(sum(stalled(q) for q in jkp)), jk_tau_r=[float(q[3]) for q in jkp],
                Q=math.pi * nd.mean() * p[3], nu0=nu0, nu0_se=se, jk_nu0=jk.tolist(), p=p)


def main():
    R = pd.read_csv(os.path.join(D, "p2_runs.csv")).sort_values("run").reset_index(drop=True)
    F = pd.read_csv(os.path.join(D, "p2_frames.csv"))
    NUD = pd.read_csv(os.path.join(D, "nud_runs.csv")).set_index("run")
    CC = pd.read_csv(os.path.join(CCB, "bistable_runs.csv")).set_index("run")
    alpha = M / 900.0
    L0s = R.L0.unique(); assert len(L0s) == 1, L0s
    L0 = float(L0s[0])
    K = kroot_own(alpha); Kt = T.k_root(alpha)
    x = K / (2 * math.pi * (L0 - 1.025))
    res["kroot"] = dict(own=K, T=Kt, reldiff=abs(K - Kt) / Kt, L0=L0, x=x)

    # (a) per-run c
    rows = []
    for r in R.itertuples():
        nd = float(NUD.loc[r.run, "nu_d"]); tr = float(r.fit_tau_r)
        c = math.sqrt(nd ** 2 + (1 / (2 * math.pi * tr)) ** 2) / x
        rows.append(dict(run=int(r.run), seed=int(r.seed), nu_d_nud=nd, nu_d_p2=float(r.nu_d), tau_r=tr, fit_A=float(r.fit_A), fit_B=float(r.fit_B),
                         fit_rms=float(r.fit_rms), stall=int(r.stall), fit_ok_p2=int(r.fit_ok), c_mine=c, c_cc=float(CC.loc[r.run, "c_run"])))
    A = pd.DataFrame(rows); A["c_rel"] = (A.c_mine - A.c_cc).abs() / A.c_cc
    A["nud_vs_p2_rel"] = (A.nu_d_nud - A.nu_d_p2).abs() / A.nu_d_nud
    A["stall_like"] = ((A.fit_A - 0.3).abs() < 1e-3) & ((A.fit_B - 0.7).abs() < 1e-3)
    A.to_csv(os.path.join(OUTD, "a_per_run_c.csv"), index=False, float_format="%.10g")
    sib = A[A.run != RUN]
    res["a"] = dict(run9_c=float(A.loc[A.run == RUN, "c_mine"].iloc[0]), sib_median=float(sib.c_mine.median()), sib_min=float(sib.c_mine.min()),
                    sib_max=float(sib.c_mine.max()), max_rel_vs_cc=float(A.c_rel.max()), max_nud_vs_p2_rel=float(A.nud_vs_p2_rel.max()),
                    stall_like_runs=A.loc[A.stall_like, "run"].tolist(), stall_col=A.loc[A.stall == 1, "run"].tolist(),
                    rms_max=float(A.fit_rms.max()), run9_nud=float(A.loc[A.run == RUN, "nu_d_nud"].iloc[0]), run9_tau_r=float(A.loc[A.run == RUN, "tau_r"].iloc[0]),
                    sib_nud_median=float(sib.nu_d_nud.median()), sib_nud_min=float(sib.nu_d_nud.min()), sib_nud_max=float(sib.nu_d_nud.max()),
                    sib_tau_median=float(sib.tau_r.median()), sib_tau_min=float(sib.tau_r.min()), sib_tau_max=float(sib.tau_r.max()),
                    # c with tau_r -> infinity (the damping term's share)
                    run9_c_nud_only=float(A.loc[A.run == RUN, "nu_d_nud"].iloc[0]) / x)

    # (b) strain per compartment
    rec = F[F.phase == 1]
    eb = rec.groupby("run")[["e_0", "e_1"]].agg(["mean", lambda v: v.std(ddof=1), "count"])
    eb.columns = ["e0_mean", "e0_sd", "n0", "e1_mean", "e1_sd", "n1"]
    eb = eb.join(R.set_index("run")[["hold_e_0", "hold_e_1", "end_e_0", "end_e_1", "start_e_0", "start_e_1", "psi6_rec_0", "psi6_rec_1",
                                    "orient_rec_0", "orient_rec_1", "end_line_0", "end_line_1"]])
    # last record snapshot and last hold snapshot from p2_frames, vs p2_runs hold_e / end_e
    lastrec = rec.sort_values("t").groupby("run").tail(1).set_index("run")
    hold = F[F.phase == 0]; lasthold = hold.sort_values("t").groupby("run").tail(1).set_index("run")
    eb["frames_lastrec_e0"] = lastrec.e_0; eb["frames_lastrec_e1"] = lastrec.e_1
    eb["frames_lasthold_e0"] = lasthold.e_0; eb["frames_lasthold_e1"] = lasthold.e_1
    eb["n_rec"] = rec.groupby("run").size(); eb["n_hold"] = hold.groupby("run").size()
    # psi6 record means and circular orientation from the psi6_t files (independent of p2_runs)
    for r in eb.index:
        ps = pd.read_csv(glob.glob(os.path.join(D, f"psi6_t_L0_*_run{r}.csv*"))[0]); pr = ps[ps.phase == "record"]
        for k in (0, 1):
            eb.loc[r, f"psi6_rec_{k}_mine"] = pr[f"psi6_{k}"].mean()
            th = np.deg2rad(pr[f"psi6_orient_deg_{k}"].to_numpy(float)) * 6
            ang = math.degrees(math.atan2(np.sin(th).mean(), np.cos(th).mean())) / 6  # in (-30, 30]
            eb.loc[r, f"orient_dev_{k}_mine"] = ang  # pin 0: deviation from 0 deg
            ow = (float(eb.loc[r, f"orient_rec_{k}"]) + 30) % 60 - 30
            eb.loc[r, f"orient_dev_{k}_from_p2"] = ow
        eb.loc[r, "n_psi_rec"] = len(pr); eb.loc[r, "t_rec_min"] = pr.t_sigma.min()
    eb.to_csv(os.path.join(OUTD, "b_strain_psi6.csv"), float_format="%.10g")
    me = eb.loc[RUN]; sb = eb.drop(index=RUN)
    res["b"] = dict(run9=dict(e_rec_left=me.e0_mean, e_rec_right=me.e1_mean, sd_left=me.e0_sd, sd_right=me.e1_sd, n_rec=int(me.n0),
                              hold_e_0=me.hold_e_0, hold_e_1=me.hold_e_1, end_e_0=me.end_e_0, end_e_1=me.end_e_1,
                              frames_lasthold=(me.frames_lasthold_e0, me.frames_lasthold_e1), frames_lastrec=(me.frames_lastrec_e0, me.frames_lastrec_e1),
                              psi6_rec_0=me.psi6_rec_0, psi6_rec_1=me.psi6_rec_1, psi6_rec_0_mine=me.psi6_rec_0_mine, psi6_rec_1_mine=me.psi6_rec_1_mine,
                              orient_rec_0=me.orient_rec_0, orient_rec_1=me.orient_rec_1, orient_dev_0_mine=me.orient_dev_0_mine,
                              orient_dev_1_mine=me.orient_dev_1_mine, end_line=(me.end_line_0, me.end_line_1)),
                    sib=dict(left_median=sb.e0_mean.median(), left_min=sb.e0_mean.min(), left_max=sb.e0_mean.max(),
                             right_median=sb.e1_mean.median(), right_min=sb.e1_mean.min(), right_max=sb.e1_mean.max(),
                             hold0=(sb.hold_e_0.median(), sb.hold_e_0.min(), sb.hold_e_0.max()), hold1=(sb.hold_e_1.median(), sb.hold_e_1.min(), sb.hold_e_1.max()),
                             end0=(sb.end_e_0.median(), sb.end_e_0.min(), sb.end_e_0.max()), end1=(sb.end_e_1.median(), sb.end_e_1.min(), sb.end_e_1.max()),
                             psi0=(sb.psi6_rec_0.median(), sb.psi6_rec_0.min(), sb.psi6_rec_0.max()), psi1=(sb.psi6_rec_1.median(), sb.psi6_rec_1.min(), sb.psi6_rec_1.max()),
                             od0=(sb.orient_dev_0_mine.median(), sb.orient_dev_0_mine.min(), sb.orient_dev_0_mine.max()),
                             od1=(sb.orient_dev_1_mine.median(), sb.orient_dev_1_mine.min(), sb.orient_dev_1_mine.max()),
                             sd0=(sb.e0_sd.median(), sb.e0_sd.min(), sb.e0_sd.max()), sd1=(sb.e1_sd.median(), sb.e1_sd.min(), sb.e1_sd.max())),
                    cmp_cc=dict(max_abs_e_rec=float(max((eb.e0_mean - CC.e_rec_0).abs().max(), (eb.e1_mean - CC.e_rec_1).abs().max())),
                                max_abs_e_sd=float(max((eb.e0_sd - CC.e_rec_sd_0).abs().max(), (eb.e1_sd - CC.e_rec_sd_1).abs().max())),
                                max_abs_psi_mine_vs_p2=float(max((eb.psi6_rec_0_mine - eb.psi6_rec_0).abs().max(), (eb.psi6_rec_1_mine - eb.psi6_rec_1).abs().max())),
                                max_abs_od_mine_vs_cc=float(max((eb.orient_dev_0_mine - CC.orient_dev_0).abs().max(), (eb.orient_dev_1_mine - CC.orient_dev_1).abs().max())),
                                max_abs_od_p2wrap_vs_cc=float(max((eb.orient_dev_0_from_p2 - CC.orient_dev_0).abs().max(), (eb.orient_dev_1_from_p2 - CC.orient_dev_1).abs().max())),
                                max_abs_lasthold_vs_p2=float(max((eb.frames_lasthold_e0 - eb.hold_e_0).abs().max(), (eb.frames_lasthold_e1 - eb.hold_e_1).abs().max())),
                                max_abs_lastrec_vs_p2=float(max((eb.frames_lastrec_e0 - eb.end_e_0).abs().max(), (eb.frames_lastrec_e1 - eb.end_e_1).abs().max()))))
    # run 9 left trajectory facts
    f9 = F[F.run == RUN].sort_values("t")
    up = f9[f9.e_0 > 0.01]
    sibF = F[(F.run != RUN) & (F.frame > 0)]
    res["b"]["run9_traj"] = dict(first_t_above_1pct=float(up.t.iloc[0]), n_pos=int((f9.e_0 > 0).sum()), n_snap=len(f9),
                                 dt_snap=float(f9.t.diff().median()), sib_pos=int((sibF.e_0 > 0).sum()), sib_n=len(sibF),
                                 first_frames=f9[["frame", "t", "phase", "e_0", "e_1", "theta0_0", "ordered_0"]].head(6).to_dict("records"))

    # (c) E1 with all 10 and without run 9
    full = e1_mass(D)
    less_used = e1_mass(D, exclude=(RUN,), nu_ref_over="used")
    less_all = e1_mass(D, exclude=(RUN,), nu_ref_over="all")
    # cross-check FU.fit_acf vs my own transcription on the full seed average
    z = np.load(os.path.join(D, "acf_runs.npz"))
    acfs = [np.asarray(z[f"run{r}"], float) for r in full["used_runs"]]
    p_own, _, _ = seed_avg_fit(acfs, full["p"][0] * 0 + float(pd.read_csv(os.path.join(D, "red_nu.csv")).dt.iloc[0]), full["nu_ref"], fitter=fit_acf_own)
    amf = pd.read_csv(os.path.join(D, "acf_mean_fit.csv")).iloc[0]
    res["c"] = {}
    for lab, pm in (("all10", full), ("no9_nuref_used", less_used), ("no9_nuref_allfolder", less_all)):
        d_ = {k: v for k, v in pm.items() if k not in ("p",)}
        d_["c"] = pm["nu0"] / x; d_["c_se"] = pm["nu0_se"] / x
        res["c"][lab] = d_
    res["c"]["fit_own_vs_FU_maxrel"] = float(np.max(np.abs(p_own - full["p"]) / np.abs(full["p"])))
    res["c"]["koa_acf_mean_fit"] = dict(tau_r=float(amf.tau_r), A=float(amf.A), B=float(amf.B), nu_ref=float(amf.nu_ref), nl=int(amf.nl))
    res["c"]["nud_runs_nu_ref"] = float(NUD.nu_ref.iloc[0])

    # cell E1 over the 5 masses (extra): own per_mass for each mass, stage-1 e1 rule
    mrows = []
    for dd in sorted(glob.glob(os.path.join(P2, CELL, "m_*"))):
        MM = int(dd.split("_")[-1]); a_ = MM / 900.0
        L0m = float(pd.read_csv(os.path.join(dd, "p2_runs.csv")).L0.iloc[0]); xm = kroot_own(a_) / (2 * math.pi * (L0m - 1.025))
        pm = e1_mass(dd, nu_ref_over="all")
        mrows.append(dict(M=MM, alpha=a_, used=pm["used"], c=pm["nu0"] / xm, sig_c=pm["nu0_se"] / xm, Q=pm["Q"], tau_r=pm["tau_r"], stalled=pm["stalled"]))
    MR = pd.DataFrame(mrows)

    def e1(cs, ss):
        cs = np.asarray(cs); ss = np.asarray(ss); w = 1 / ss ** 2; m = (w * cs).sum() / w.sum(); ch = (w * (cs - m) ** 2).sum(); dof = len(cs) - 1
        return dict(cs=float(m), sig=float(w.sum() ** -0.5 * max(1.0, math.sqrt(ch / dof))), chi2=float(ch), dof=dof)
    cA = MR.c.to_numpy().copy(); sA = MR.sig_c.to_numpy().copy()
    i2 = int(np.where(MR.M.to_numpy() == M)[0][0])
    cB = cA.copy(); sB = sA.copy(); cB[i2] = less_used["nu0"] / x; sB[i2] = less_used["nu0_se"] / x
    res["c"]["cell_masses"] = MR.to_dict("records")
    res["c"]["cell_e1_all10"] = e1(cA, sA); res["c"]["cell_e1_no9"] = e1(cB, sB)
    PM = pd.read_csv(os.path.join(HS, "experiments_gen3_p2_261016", "report_261016", "tables", "p2_masses.csv"))
    PM = PM[PM.cell == CELL][["M", "c", "sig_c", "used", "tau_r"]]
    res["c"]["cc_p2_masses"] = PM.to_dict("records")

    # (d) the scan
    sc = []
    for dd in sorted(glob.glob(os.path.join(P2, "*", "m_*"))):
        cid = dd.split(os.sep)[-2]; MM = int(dd.split("_")[-1])
        FF = pd.read_csv(os.path.join(dd, "p2_frames.csv")); FF = FF[FF.phase == 1]
        for r_, g in FF.groupby("run"):
            for k in (0, 1):
                v = g[f"e_{k}"].to_numpy(float)
                sc.append(dict(cell=cid, M=MM, run=int(r_), side=k, n=len(v), e_rec=v.mean(), sd_t=v.std(ddof=1)))
    S = pd.DataFrame(sc)
    S["med"] = S.groupby(["cell", "side"]).e_rec.transform("median")
    S["thermal"] = S.groupby(["cell", "side"]).sd_t.transform("median")
    S["nruns_cellside"] = S.groupby(["cell", "side"]).e_rec.transform("size")
    S["dev"] = S.e_rec - S.med; S["dev_sd"] = S.dev / S.thermal
    S.to_csv(os.path.join(OUTD, "d_scan.csv"), index=False, float_format="%.10g")
    far = S[S.dev.abs() > 0.01]
    srt = S.dev_sd.abs().sort_values(ascending=False)
    top = S.loc[srt.index[:5], ["cell", "M", "run", "side", "e_rec", "med", "thermal", "dev", "dev_sd"]]
    srt_pct = S.dev.abs().sort_values(ascending=False)
    top_pct = S.loc[srt_pct.index[:3], ["cell", "M", "run", "side", "e_rec", "med", "dev", "dev_sd"]]
    CS = pd.read_csv(os.path.join(CCB, "bistable_scan.csv"))
    mg = S.merge(CS, on=["cell", "M", "run", "side"], suffixes=("", "_cc"))
    res["d"] = dict(n=len(S), n_cells=S.cell.nunique(), nruns_per_cellside=sorted(S.nruns_cellside.unique().tolist()), n_rec_snaps=sorted(S.n.unique().tolist()),
                    n_far=len(far), far=far[["cell", "M", "run", "side", "e_rec", "med", "thermal", "dev", "dev_sd"]].to_dict("records"),
                    top5_sd=top.to_dict("records"), top3_pct=top_pct.to_dict("records"),
                    merged=len(mg), max_abs_diff=dict(e_rec=float((mg.e_rec - mg.e_rec_cc).abs().max()), sd_t=float((mg.sd_t - mg.sd_t_cc).abs().max()),
                                                     med=float((mg.med - mg.median_cell_side).abs().max()), thermal=float((mg.thermal - mg.thermal_cell_side).abs().max()),
                                                     dev_sd=float((mg.dev_sd - mg.dev_in_thermal_sd).abs().max())),
                    # thermal sd variants (robustness): median of sd within the (cell, side, M) group instead of the whole cell-side
                    )
    # robustness: per (cell, M, side) median and thermal sd
    S["med_M"] = S.groupby(["cell", "M", "side"]).e_rec.transform("median"); S["th_M"] = S.groupby(["cell", "M", "side"]).sd_t.transform("median")
    S["dev_sd_M"] = (S.e_rec - S.med_M) / S.th_M
    srtM = S.dev_sd_M.abs().sort_values(ascending=False)
    res["d"]["per_mass_variant_top3"] = S.loc[srtM.index[:3], ["cell", "M", "run", "side", "e_rec", "med_M", "th_M", "dev_sd_M"]].to_dict("records")
    res["d"]["per_mass_variant_n_far"] = int(((S.e_rec - S.med_M).abs() > 0.01).sum())

    with open(os.path.join(OUTD, "verify_results.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else (list(o) if isinstance(o, (tuple, np.ndarray)) else str(o)))
    print(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else (list(o) if isinstance(o, (tuple, np.ndarray)) else str(o))))


if __name__ == "__main__":
    main()
