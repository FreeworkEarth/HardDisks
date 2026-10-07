#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.13; third plan-author decision of 2026-10-07, items 1, 4 and 5): INFORMATION ONLY, on the
existing data (Test T, the campaign anchor, the replay). No verdict: Test T's verdict (FAIL, sec. 4.4.12) stands, and the
registered estimators and results are not changed by anything here.
  item 1   the tail probability of Test T's two largest |z| under the null, by simulation: 1e6 Gaussian draws of the nine per-mass
           means of both policies around their pooled mean with the observed SEs, the registered estimators (resched_testT_261007
           .estimators) applied to every draw, numpy default_rng(20261011); as a check, 1e5 relabelings of the trajectories within
           each mass pool (default_rng(20261012)), and the formula for eleven independent normal z.
  item 4a  a refined frequency per trajectory from acf_runs.npz (the mean-removed position ACF to 20 predicted periods, reduce_B.py):
           a damped cosine with free phase plus the slow mode, C(t)/C(0) = A e^{-t/tau_T} + B e^{-t/tau_r} cos(omega t + phi),
           least squares over the whole stored ACF (bounds of paper1_confinement_results_261004._fit_acf, phi in [-pi/2, pi/2],
           start omega = 2 pi x the mass's pooled argmax mean, the same for both policies); nu_d = omega/2 pi. The phase is free
           because the ACF of a noise-driven damped oscillator is e^{-t/tau}[cos(w_d t) + sin(w_d t)/(w_d tau)]: a cosine without
           phase would pull omega down by ~1/(4Q^2). The other option of the decision, the mean zero-crossing period over the first
           10 periods, is reported by its crossing count only: the slow mode A e^{-t/tau_T} keeps the ACF positive after about one
           period at the light masses. Printed: per-seed SD of both estimators per mass; the eleven numbers with the refined
           estimator next to the argmax one.
  item 4b  M = 300: minimal minus legacy in the first and the second 50 seeds (seed-list order), each with z; the same-seed correlation.
  item 4c  per-mass c_s,M = nu_M/x_M: campaign anchor (25 seeds, legacy, 279282b), replay (25, minimal, 73fc07f), Test T minimal,
           Test T legacy; z campaign vs Test T legacy per mass and for the registered c_s; each mass's share of the c_s difference
           (the unweighted slope is sum_M w_M c_s,M with w_M = x_M^2 / sum x^2).
  item 5   per mass of this cell: the gas-inertia fraction (2/3) N_s m / M, the exact root K of cot K = alpha K that the registered
           x_M = K/(2 pi L_eff) uses, the effective-mass root (alpha + 1/3)^(-1/2) and their ratio, the weights w_M, the damping Q
           from the fit of the seed-averaged ACF (Test T legacy, same model as 4a), the spectral-peak shift of a simple damped
           oscillator -1/(4 Q^2), and the Test T legacy deficits c_s,M / c_s,weighted - 1 with the argmax, nu_d and
           nu_0 = sqrt(nu_d^2 + (1/(2 pi tau_r))^2) frequencies (each against its own weighted slope).
usage (from hspist3/): python3 validation/resched_testT_followup_261007.py
"""
import contextlib, io, math, os, sys, warnings
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm, chi2 as CHI2
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import tests_20260913 as T
import resched_testT_261007 as TT
from paper1_populate_cs_err_20261002 import slope_with_errors

CID, NDRAW, NPERM, RNG_G, RNG_P = TT.CID, 1_000_000, 100_000, 20261011, 20261012
POL = ("minimal", "legacy")


def _g(t, A, tT, B, tr, om, ph):
    return A * np.exp(-t / tT) + B * np.exp(-t / tr) * np.cos(om * t + ph)


def fit_acf(a, dt, nu_ref):
    c = a / a[0]; lag = np.arange(len(c)) * dt; om0 = 2 * math.pi * nu_ref; P = 1 / nu_ref
    p, _ = curve_fit(_g, lag, c, p0=[0.3, 2 * P, 0.7, 5 * P, om0, 0.0],
                     bounds=([0, 2 * dt, 0, 2 * dt, 0.2 * om0, -math.pi / 2], [1.5, 1e7, 1.5, 1e7, 5 * om0, math.pi / 2]), maxfev=200000)
    return p


def crossings(a, dt, nu_ref):
    c = a / a[0]; n = min(len(c), int(10 / (nu_ref * dt)) + 1)
    return int(np.sum(np.signbit(c[1:n]) != np.signbit(c[:n - 1])))


def load():
    with contextlib.redirect_stdout(io.StringIO()):
        c = {x["cid"]: x for x in R.cells()}[CID]; co = dict(c); R.method_B(co)
    Ms, Ns, x = c["Ms"], c["Ns"], np.asarray(co["x"], float)
    D = {}
    for pol in POL:
        for M in Ms:
            d = os.path.join(TT.LOC, TT.REL_T, pol, CID, f"m_{M}")
            red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
            z = np.load(os.path.join(d, "acf_runs.npz"))
            D[(pol, M)] = dict(red=red, acf=[z[f"run{r}"].astype(float) for r in red["run"]], dt=float(red["dt"].iloc[0]))
    rb = os.path.join(R.REL_B, CID); camp, rep = {}, {}
    for out, root in ((camp, HS), (rep, os.path.join(HS, "experiments_resched_gate_261005"))):
        for M in Ms:
            d = os.path.join(root, rb, f"m_{M}"); red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
            z = np.load(os.path.join(d, "acf_runs.npz"))
            out[M] = dict(nu=red["nu"].to_numpy(float), acf=[z[f"run{r}"].astype(float) for r in red["run"]], dt=float(red["dt"].iloc[0]))
    return c, co, Ms, Ns, x, D, camp, rep


def zrows(nu_min, nu_leg, Ms, x, Ns):
    out = {}
    for name, a, sa, b, sb in TT.eleven(nu_min, nu_leg, Ms, x, Ns):
        s = math.hypot(sa, sb); out[name] = (a - b, s, b, (a - b) / s)
    return out


def tail_stats(Z, z1, z2, zs, chi_obs, nmass):
    A = np.sort(np.abs(Z), axis=1)[:, ::-1]
    return np.array([(A[:, 0] >= z1).sum(), ((A[:, 0] >= z1) & (A[:, 1] >= z2)).sum(), (A[:, 0] >= zs).sum(),
                     ((Z[:, 2:2 + nmass] ** 2).sum(1) >= chi_obs).sum()], float)


def item1(nu, Ms, x, Ns):
    print("## Item 1 -- Test T's largest |z| under the null, by simulation\n")
    obs = zrows(nu["minimal"], nu["legacy"], Ms, x, Ns); names = list(obs)
    zo = np.array([obs[k][3] for k in names]); order = np.argsort(-np.abs(zo)); z1, z2 = abs(zo[order[0]]), abs(zo[order[1]])
    zs = norm.isf(TT.FW / (2 * TT.NNUM)); chi_obs = float(sum(obs[f"nu M={M}"][3] ** 2 for M in Ms))
    print(f"observed: largest |z| = {z1:.4f} ({names[order[0]]}), second = {z2:.4f} ({names[order[1]]}); nine-mass chi2 = {chi_obs:.2f}; z* = {zs:.4f}\n")
    mm = np.array([nu["minimal"][M].mean() for M in Ms]); ml = np.array([nu["legacy"][M].mean() for M in Ms])
    sm = np.array([nu["minimal"][M].std(ddof=1) / math.sqrt(len(nu["minimal"][M])) for M in Ms])
    sl = np.array([nu["legacy"][M].std(ddof=1) / math.sqrt(len(nu["legacy"][M])) for M in Ms])
    mu = 0.5 * (mm + ml)

    def zs_of(a_m, s_m, a_l, s_l):
        cm, cme, km, kme = TT.estimators(a_m, s_m, Ms, x, Ns); cl, cle, kl, kle = TT.estimators(a_l, s_l, Ms, x, Ns)
        return np.column_stack([(km - kl) / np.hypot(kme, kle), (cm - cl) / np.hypot(cme, cle), (a_m - a_l) / np.hypot(s_m, s_l)])

    rng = np.random.default_rng(RNG_G); acc = np.zeros(4)
    for _ in range(NDRAW // 100000):
        a_m = mu + sm * rng.standard_normal((100000, len(Ms))); a_l = mu + sl * rng.standard_normal((100000, len(Ms)))
        acc += tail_stats(zs_of(a_m, sm, a_l, sl), z1, z2, zs, chi_obs, len(Ms))
    pg = acc / NDRAW
    rng = np.random.default_rng(RNG_P); accp = np.zeros(4); CH = 10000
    for _ in range(NPERM // CH):
        a_m = np.empty((CH, len(Ms))); a_l = np.empty_like(a_m); s_m = np.empty_like(a_m); s_l = np.empty_like(a_m)
        for j, M in enumerate(Ms):
            pool = np.r_[nu["minimal"][M], nu["legacy"][M]]; na = len(nu["minimal"][M])
            v = pool[np.argsort(rng.random((CH, len(pool))), axis=1)]; pa, pb = v[:, :na], v[:, na:]
            a_m[:, j], a_l[:, j] = pa.mean(1), pb.mean(1)
            s_m[:, j], s_l[:, j] = pa.std(1, ddof=1) / math.sqrt(na), pb.std(1, ddof=1) / math.sqrt(len(pool) - na)
        accp += tail_stats(zs_of(a_m, s_m, a_l, s_l), z1, z2, zs, chi_obs, len(Ms))
    pp = accp / NPERM
    p1, p2 = 2 * norm.sf(z1), 2 * norm.sf(z2); n = TT.NNUM
    ind = [1 - (1 - p1) ** n, 1 - (1 - p1) ** n - n * p1 * (1 - p2) ** (n - 1), 1 - (1 - 2 * norm.sf(zs)) ** n, CHI2.sf(chi_obs, len(Ms))]
    print(f"| event under the null | Gaussian, {NDRAW:.0e} draws, registered estimators | relabelings, {NPERM:.0e} | eleven independent normal z |\n|---|---|---|---|")
    labs = [f"largest abs(z) >= {z1:.2f}", f"largest abs(z) >= {z1:.2f} AND second largest >= {z2:.2f}",
            f"largest abs(z) >= z* = {zs:.4f} (the design's family-wise rate)", f"nine-mass chi2 >= {chi_obs:.2f}"]
    for lab, a, b, c in zip(labs, pg, pp, ind):
        print(f"| {lab} | {a:.5f} (1 in {1 / a:.0f}) | {b:.5f} (1 in {1 / b:.0f}) | {c:.5f} (1 in {1 / c:.0f}) |" if a > 0 and b > 0 else f"| {lab} | {a:.5f} | {b:.5f} | {c:.5f} |")
    print("\nThe joint event is defined from the observed values after seeing them (post hoc); its probability is not a test level.")


def item4a(D, nu, Ms, x, Ns, camp, rep):
    print("\n## Item 4a -- refined frequency per trajectory (damped cosine with free phase plus slow mode) next to the argmax estimator\n")
    warnings.filterwarnings("ignore")
    ref = {M: 0.5 * (nu["minimal"][M].mean() + nu["legacy"][M].mean()) for M in Ms}
    nd, n0, fails, cross, fit_mean = {p: {} for p in POL}, {p: {} for p in POL}, {}, {}, {}
    for pol in POL:
        for M in Ms:
            e = D[(pol, M)]; v = []; ff = 0; cc = []
            for a in e["acf"]:
                cc.append(crossings(a, e["dt"], ref[M]))
                try: v.append(fit_acf(a, e["dt"], ref[M])[4] / (2 * math.pi))
                except Exception: v.append(np.nan); ff += 1
            v = np.array(v); nd[pol][M] = v[np.isfinite(v)]; fails[(pol, M)] = ff; cross[(pol, M)] = (int(np.median(cc)), min(cc), max(cc))
            p = fit_acf(np.mean([a[:min(len(b) for b in e["acf"])] for a in e["acf"]], axis=0), e["dt"], ref[M])   # seed-averaged ACF
            fit_mean[(pol, M)] = dict(nud=p[4] / (2 * math.pi), tau_r=p[3], A=p[0], B=p[2], tau_T=p[1], phi=p[5])
            n0[pol][M] = math.sqrt(nd[pol][M].mean() ** 2 + (1 / (2 * math.pi * p[3])) ** 2)
    for src in (camp, rep):
        for M in Ms:
            v = []
            for a in src[M]["acf"]:
                try: v.append(fit_acf(a, src[M]["dt"], ref[M])[4] / (2 * math.pi))
                except Exception: v.append(np.nan)
            v = np.array(v); src[M]["nud"] = v[np.isfinite(v)]; src[M]["fail"] = int((~np.isfinite(v)).sum())
    print(f"campaign and replay nu_d fit failures: {sum(camp[M]['fail'] for M in Ms)} and {sum(rep[M]['fail'] for M in Ms)} of {sum(len(camp[M]['acf']) for M in Ms)} each\n")
    print("| policy | M | trajectories | fit failures | per-seed SD, argmax [% of nu] | per-seed SD, nu_d [% of nu] | SD ratio argmax/nu_d | "
          "zero crossings in the first 10 periods: median (min-max), expected 20 |\n|---|---|---|---|---|---|---|---|")
    for pol in POL:
        for M in Ms:
            a = nu[pol][M]; b = nd[pol][M]
            print(f"| {pol} | {M} | {len(a)} | {fails[(pol, M)]} | {100 * a.std(ddof=1) / ref[M]:.3f} | {100 * b.std(ddof=1) / ref[M]:.3f} | "
                  f"{a.std(ddof=1) / b.std(ddof=1):.2f} | {cross[(pol, M)][0]} ({cross[(pol, M)][1]}-{cross[(pol, M)][2]}) |")
    za = zrows(nu["minimal"], nu["legacy"], Ms, x, Ns); zd = zrows(nd["minimal"], nd["legacy"], Ms, x, Ns)
    print("\n### The eleven numbers, minimal minus legacy, with both estimators (registered estimators of k_S^dyn and c_s applied to each)\n")
    print("| number | argmax: relative [%] | 95 % interval [%] | z | nu_d: relative [%] | 95 % interval [%] | z |\n|---|---|---|---|---|---|---|")
    for k in za:
        (d1, s1, b1, z1), (d2, s2, b2, z2) = za[k], zd[k]
        print(f"| {k} | {100 * d1 / b1:+.3f} | [{100 * (d1 - 1.96 * s1) / b1:+.3f}, {100 * (d1 + 1.96 * s1) / b1:+.3f}] | {z1:+.2f} | "
              f"{100 * d2 / b2:+.3f} | [{100 * (d2 - 1.96 * s2) / b2:+.3f}, {100 * (d2 + 1.96 * s2) / b2:+.3f}] | {z2:+.2f} |")
    ca = sum(za[f"nu M={M}"][3] ** 2 for M in Ms); cd = sum(zd[f"nu M={M}"][3] ** 2 for M in Ms)
    print(f"\nnine-mass chi2: argmax {ca:.2f} (nominal p {CHI2.sf(ca, len(Ms)):.3f}); nu_d {cd:.2f} (nominal p {CHI2.sf(cd, len(Ms)):.3f})")
    return nd, n0, fit_mean


def item4b(D, nu, nd):
    print("\n## Item 4b -- M = 300: the two halves of the seed list, and the same-seed correlation\n")
    print("| estimator | seeds (r) | minimal mean | legacy mean | relative difference [%] | z |\n|---|---|---|---|---|---|")
    for lab, src in (("argmax (registered)", nu), ("nu_d (4a)", nd)):
        for lo, hi in ((0, 50), (50, 100), (0, 100)):
            a, b = src["minimal"][300][lo:hi], src["legacy"][300][lo:hi]
            d = a.mean() - b.mean(); s = math.hypot(a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b)))
            print(f"| {lab} | {lo}-{hi - 1} | {a.mean():.7f} | {b.mean():.7f} | {100 * d / b.mean():+.3f} | {d / s:+.2f} |")
    for lab, src in (("argmax", nu), ("nu_d", nd)):
        r = np.corrcoef(src["minimal"][300], src["legacy"][300])[0, 1]; n = len(src["minimal"][300])
        print(f"\nsame-seed correlation of nu between the policies, M = 300, {lab}: r = {r:+.3f} (n = {n}; |r| > {1.96 / math.sqrt(n - 3):.3f} "
              f"would be outside the 95 % range of r = 0, Fisher z)")


def item4c(co, Ms, Ns, x, nu, nd, camp, rep):
    print("\n## Item 4c -- per-mass implied sound speed: campaign anchor, replay, Test T\n")
    sets = [("campaign (legacy, 279282b)", {M: camp[M]["nu"] for M in Ms}), ("replay (minimal, 73fc07f)", {M: rep[M]["nu"] for M in Ms}),
            ("Test T minimal", nu["minimal"]), ("Test T legacy", nu["legacy"])]
    cs = {lab: ({M: v[M].mean() / xm for M, xm in zip(Ms, x)}, {M: v[M].std(ddof=1) / math.sqrt(len(v[M])) / xm for M, xm in zip(Ms, x)},
                {M: len(v[M]) for M in Ms}) for lab, v in sets}
    w = x * x / (x * x).sum()
    print("| M | weight w_M of the unweighted slope [%] | " + " | ".join(f"{lab}: c_s,M +- SE (n)" for lab, _ in sets) +
          " | z campaign - Test T legacy | share of the c_s difference [%] |\n|---|---|" + "---|" * len(sets) + "---|---|")
    c0, cl = cs[sets[0][0]], cs[sets[3][0]]
    dcs = sum(wj * (c0[0][M] - cl[0][M]) for wj, M in zip(w, Ms)); zz = []
    for wj, M in zip(w, Ms):
        z = (c0[0][M] - cl[0][M]) / math.hypot(c0[1][M], cl[1][M]); zz.append(z)
        print(f"| {M} | {100 * wj:.1f} | " + " | ".join(f"{cs[lab][0][M]:.5f} +- {cs[lab][1][M]:.5f} ({cs[lab][2][M]})" for lab, _ in sets) +
              f" | {z:+.2f} | {100 * wj * (c0[0][M] - cl[0][M]) / dcs:+.1f} |")
    zmax = max(abs(v) for v in zz); chz = float(np.sum(np.square(zz)))
    print(f"\ncampaign vs Test T legacy (both the legacy path, different seeds): largest abs(z) over the nine masses {zmax:.2f}, "
          f"P(largest >= that | nine independent normal z) = {1 - (1 - 2 * norm.sf(zmax)) ** len(Ms):.4f}; nine-mass chi2 {chz:.2f} "
          f"(nominal p {CHI2.sf(chz, len(Ms)):.4f})")
    print("\n| data | registered c_s (unweighted slope) | c_s_err (unscaled) | chi2_red | c_s_err_scaled | weighted c_s +- SE | unweighted / weighted - 1 [%] |\n|---|---|---|---|---|---|---|")
    reg = {}
    for lab, v in sets:
        y = np.array([v[M].mean() for M in Ms]); sy = np.array([v[M].std(ddof=1) / math.sqrt(len(v[M])) for M in Ms])
        s, err, errs, chi2r = slope_with_errors(x, y, sy); wt = 1 / sy ** 2; sw = float((wt * x * y).sum() / (wt * x * x).sum())
        reg[lab] = (s, err, errs)
        print(f"| {lab} | {s:.5f} | {err:.5f} | {chi2r:.2f} | {errs:.5f} | {sw:.5f} +- {1 / math.sqrt(float((wt * x * x).sum())):.5f} | {100 * (s / sw - 1):+.3f} |")
    (s0, e0, es0), (s1, e1, es1) = reg[sets[0][0]], reg[sets[3][0]]
    print(f"\ncampaign minus Test T legacy, registered c_s: {s0 - s1:+.5f}; z = {(s0 - s1) / math.hypot(e0, e1):+.2f} with the unscaled errors, "
          f"{(s0 - s1) / math.hypot(es0, es1):+.2f} with the scaled ones; sum over masses of w_M x (c_s,M difference) = {dcs:+.5f} (identity check)")
    print("\n### The same per-mass comparison with the refined frequency nu_d (item 4a fit on the campaign's and the replay's own ACFs)\n")
    print("| M | campaign nu_d (25) | Test T legacy nu_d | z campaign - Test T legacy, nu_d | z, argmax (above) | replay nu_d (25) | Test T minimal nu_d | z replay - Test T minimal, nu_d |\n|---|---|---|---|---|---|---|---|")
    zc, zr = [], []
    for M, za_ in zip(Ms, zz):
        a, b, r_, t_ = camp[M]["nud"], nd["legacy"][M], rep[M]["nud"], nd["minimal"][M]
        z1 = (a.mean() - b.mean()) / math.hypot(a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b)))
        z2 = (r_.mean() - t_.mean()) / math.hypot(r_.std(ddof=1) / math.sqrt(len(r_)), t_.std(ddof=1) / math.sqrt(len(t_))); zc.append(z1); zr.append(z2)
        print(f"| {M} | {a.mean():.7f} | {b.mean():.7f} | {z1:+.2f} | {za_:+.2f} | {r_.mean():.7f} | {t_.mean():.7f} | {z2:+.2f} |")
    for lab, v in (("campaign vs Test T legacy, nu_d", zc), ("replay vs Test T minimal, nu_d", zr)):
        m = max(abs(q) for q in v); ch = float(np.sum(np.square(v)))
        print(f"{lab}: largest abs(z) {m:.2f} (P(largest >= that | nine independent normal z) = {1 - (1 - 2 * norm.sf(m)) ** len(Ms):.4f}); "
              f"nine-mass chi2 {ch:.2f} (nominal p {CHI2.sf(ch, len(Ms)):.4f})")
    return cs


def item5(c, co, Ms, Ns, x, nu, nd, n0, fit_mean):
    print("\n## Item 5 -- the gas inertia in the registered model, and the light-mass deficit of the plain-fluid baseline (Test T legacy)\n")
    w = x * x / (x * x).sum(); damp = {}
    pcsv = os.path.join(T.PLOTS, "261004_p1_confinement_damping.csv")
    if os.path.exists(pcsv):
        dd = pd.read_csv(pcsv); dd = dd[dd["cell"] == CID]; damp = {int(r.M): float(r.Q) for r in dd.itertuples()}

    def deficits(y, per_seed):
        """c_s,M / weighted slope - 1, with the SE of the mean from the per-seed values that the frequency is built from."""
        y = np.array([y[M] for M in Ms]); sy = np.array([per_seed[M].std(ddof=1) / math.sqrt(len(per_seed[M])) for M in Ms])
        wt = 1 / sy ** 2; sw = float((wt * x * y).sum() / (wt * x * x).sum())
        return y / x / sw - 1, sy / x / sw, sw
    da, sa, swa = deficits({M: nu["legacy"][M].mean() for M in Ms}, nu["legacy"])
    dn, sn, swn = deficits({M: nd["legacy"][M].mean() for M in Ms}, nd["legacy"])
    d0, s0_, sw0 = deficits(n0["legacy"], nd["legacy"])        # nu_0 = f(mean nu_d, tau_r of the averaged ACF): the nu_d seed SE
    print("| M | alpha = M/(2 N_s m) | gas inertia (2/3) N_s m / M [%] | K, cot K = alpha K (registered) | K_eff = (alpha + 1/3)^(-1/2) | K_eff/K - 1 [%] | "
          "w_M [%] | Q, campaign (sec. 13 model) | Q, Test T legacy (4a model, averaged ACF) | -1/(4Q^2) [%] (Test T Q) | "
          "deficit, argmax [%] | deficit, nu_d [%] | deficit, nu_0 [%] |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for j, M in enumerate(Ms):
        al = M / (2.0 * c["Ns"]); K = T.k_root(al); Ke = (al + 1 / 3) ** -0.5; fm = fit_mean[("legacy", M)]
        Q = math.pi * fm["nud"] * fm["tau_r"]
        print(f"| {M} | {al:g} | {100 * (2 / 3) * c['Ns'] / M:.1f} | {K:.5f} | {Ke:.5f} | {100 * (Ke / K - 1):+.3f} | {100 * w[j]:.1f} | "
              f"{damp.get(M, float('nan')):.1f} | {Q:.1f} | {-100 / (4 * Q * Q):+.3f} | {100 * da[j]:+.3f} +- {100 * sa[j]:.3f} | "
              f"{100 * dn[j]:+.3f} +- {100 * sn[j]:.3f} | {100 * d0[j]:+.3f} +- {100 * s0_[j]:.3f} |")
    print(f"\nweighted slopes (Test T legacy): argmax {swa:.5f}, nu_d {swn:.5f}, nu_0 {sw0:.5f}; deficit = c_s,M / (that slope) - 1, SE from the seeds; "
          "nu_0 = sqrt(nu_d^2 + (1/(2 pi tau_r))^2) with tau_r from the fit of the mass's seed-averaged ACF")
    for lab, d, s in (("argmax", da, sa), ("nu_d", dn, sn), ("nu_0", d0, s0_)):
        ch = float(((d / s) ** 2).sum())
        print(f"single-c_s chi2 at the weighted slope, {lab}: {ch:.1f} (8 dof, p {CHI2.sf(ch, 8):.3g})")


def main():
    c, co, Ms, Ns, x, D, camp, rep = load()
    nu = {pol: {M: D[(pol, M)]["red"]["nu"].to_numpy(float) for M in Ms} for pol in POL}
    print("# Test T follow-up (261012 sec. 4.4.13, third decision of 2026-10-07, items 1, 4, 5) -- information only, no verdict\n")
    item1(nu, Ms, x, Ns)
    nd, n0, fit_mean = item4a(D, nu, Ms, x, Ns, camp, rep)
    item4b(D, nu, nd)
    item4c(co, Ms, Ns, x, nu, nd, camp, rep)
    item5(c, co, Ms, Ns, x, nu, nd, n0, fit_mean)


if __name__ == "__main__":
    main()
