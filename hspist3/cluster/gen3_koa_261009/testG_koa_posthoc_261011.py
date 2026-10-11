#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 27, task B; 261012 sec. 4.7.56): POST HOC, INFORMATION -- tables on the fetched data of TEST
G ON KOA. They change no verdict: TEST G (KOA): FAIL (DIFFERENCE) at M = 50 stands (sec. 4.7.54, decision 27). Committed before its first
run on the data, with a self-test (--selftest) on synthetic data with known answers. Reads only; writes nothing.
DATA: hspist3/experiments_gen3_gate_koa_261009/testG/<gen2|gen3>/epi8_H_H10_L10/m_<M>/ (pi/8, N = 100, N_s = 50, L0 = H = 10; the cell of
Test T and of methods sec. 16), M = 50, 300, 1500, 2000: red_nu.csv (per trajectory: run, seed, the registered argmax nu, nu_pred, dt, n),
acf_runs.npz (per trajectory the mean-removed ACF of the divider displacement in sigma, unbiased normalization, to 20 predicted periods;
reduce_B.py) and run.log (the seeding lines and the [EDMD-ENERGY] lines of HD_KE_TRACE=1).
THE ESTIMATORS ARE THE PROJECT'S, imported, not rewritten:
  nu_arg   the registered argmax (red_nu.csv's nu);
  nu_d     validation/resched_testT_followup_261007.fit_acf (sec. 4.4.13 item 4a: C(t)/C(0) = A e^{-t/tau_T} + B e^{-t/tau_r} cos(omega t + phi),
           free phase) on each trajectory's ACF, started at 2 pi x the mass's pooled argmax mean of both engines (item 4a pooled the two
           policies the same way); nu_d = omega / 2 pi;
  tau_r, Q the same fit on the mass's seed-averaged ACF (truncated to the shortest); Q = pi nu_d(averaged fit) tau_r (item 5);
  nu_0     methods sec. 16 (item 5): sqrt((mean nu_d)^2 + (1/(2 pi tau_r))^2);
  x_M      tests_20260913.x_of(M, L0) = K(alpha) / (2 pi L_eff), cot K = alpha K, alpha = M / (2 N_s m) (the registered x).
SEs: delete-one-group jackknife over the trajectories in run order (groups of 10 at n = 400, 5 at n = 100), every quantity recomputed in
each replicate (the averaged-ACF fit included).
B1  per engine and mass: mean nu_arg, mean nu_d, tau_r, Q, nu_0, each with its SE; c_k = nu_k / x_M; delta_M = c_M / c_ref - 1, c_ref the
    inverse-variance weighted mean of c_M over M >= 300 (300, 1500, 2000) of the same engine and estimator (= sec. 16's weighted slope
    through the origin, restricted to M >= 300 as decision 27 says); SE by the delta method, the reference's own error included (the
    masses are independent samples); delta_50 per engine with SE, its z against 0 and against sec. 16's -0.424 +- 0.078 % (nu_0, Test T
    legacy, nine-mass reference), and gen3 - gen2 per estimator.
B2  the argmax bins: nu_arg = k df, df = 1 / (n dt) (no zero padding: tests_20260913._spectrum); per mass the counts per bin in both arms,
    the bin width in % of nu, each arm's SD in bins, the mean difference in bins, and whether the two arms' counts agree bin by bin.
B3  pairing by start: per mass the seeds common to both arms, whether each pair's seeding lines (the [KE-AUDIT] lines and "[EDMD-ENERGY]
    0 state loaded" up to E_tot) are byte-identical, then the paired difference gen3 - gen2 per seed for nu_arg and nu_d: mean, SE
    (SD / sqrt(pairs)), the same-seed correlation, and the unpaired SE for comparison.
B4  the divider's kinetic temperature against equipartition at fixed total energy over f = 2N + 1 = 201 quadratic degrees of freedom,
    <M V^2> = 2 E_tot / f:
    (i)  exact, one instant per trajectory: "[EDMD-ENERGY] 2 end of record" gives E_div = M V^2 / 2 and E_tot; mean of 2 E_div with its SE
         over trajectories, per engine and mass and pooled over masses;
    (ii) over the record, a proxy: the mean square of the divider's average velocity over one sample interval, M <(dx/dt)^2> =
         M 2 (C(0) - C(k dt)) / (k dt)^2 from the fetched ACF (k = 1 and 2), per trajectory; biased low against <M V^2> by the velocity's
         decorrelation within k dt (an exact oscillation alone gives 1 - (omega k dt)^2 / 12, about 0.3 % at k = 1 here); the same estimator in
         both arms, so its gen3 - gen2 difference is the informative number.
B5  per engine and mass the [EDMD-HEALTH] lines in the fetched logs (printed only when a gen2 counter is nonzero).
usage (from hspist3/): python3 cluster/gen3_koa_261009/testG_koa_posthoc_261011.py [--selftest] [<data folder>]
"""
import math, os, re, sys, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import resched_testT_followup_261007 as FU                      # fit_acf: the sec. 4.4.13 item 4a model
MASSES, ENG, CELL, L0, NS, NGAS = (50, 300, 1500, 2000), ("gen2", "gen3"), "epi8_H_H10_L10", 10.0, 50, 100
F_DOF = 2 * NGAS + 1
SEC16 = (-0.424, 0.078)                                           # methods sec. 16, nu_0 deficit at M = 50, Test T legacy [%]
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------- the pieces (each checked by the self-test)
def groups_of(n):
    g = 10 if n >= 200 else 5
    return [np.arange(i, min(i + g, n)) for i in range(0, n, g)]


def jackknife(stat, n):
    """delete-one-group jackknife of a vector statistic: (estimate on all, SE, replicates that failed)"""
    full = np.atleast_1d(np.asarray(stat(np.arange(n)), float))
    reps = np.array([np.atleast_1d(stat(np.setdiff1d(np.arange(n), g))) for g in groups_of(n)], float)
    ok = np.all(np.isfinite(reps), axis=1); r = reps[ok]; m = len(r)
    return full, np.sqrt((m - 1) / m * ((r - r.mean(axis=0)) ** 2).sum(axis=0)), int((~ok).sum())


def avg_fit(acfs, dt, nu_ref):
    """(nu_d, tau_r) of the 4a fit of the averaged ACF"""
    try: p = FU.fit_acf(np.mean(acfs, axis=0), dt, nu_ref)
    except Exception: return math.nan, math.nan
    return p[4] / (2 * math.pi), p[3]


def nu_d_each(acfs, dt, nu_ref):
    v = []
    for a in acfs:
        try: v.append(FU.fit_acf(a, dt, nu_ref)[4] / (2 * math.pi))
        except Exception: v.append(math.nan)
    return np.array(v)


def block_stats(nu_arg, nud, acfs, dt, nu_ref):
    """jackknifed vector: mean nu_arg, mean nu_d, tau_r, Q, nu_0"""
    ok = np.isfinite(nud)

    def stat(idx):
        nd, tr = avg_fit(acfs[idx], dt, nu_ref); md = float(nud[idx[ok[idx]]].mean())
        return [float(nu_arg[idx].mean()), md, tr, math.pi * nd * tr, math.sqrt(md ** 2 + (1 / (2 * math.pi * tr)) ** 2)]
    return jackknife(stat, len(nu_arg))


def deltas(c, s, ref_masses):
    """delta_M = c_M / c_ref - 1 and its SE (delta method), c_ref the inverse-variance weighted mean over ref_masses"""
    w = {M: 1 / s[M] ** 2 for M in ref_masses}; W = sum(w.values())
    cref = sum(w[M] * c[M] for M in ref_masses) / W
    out = {}
    for M in c:
        g = {j: -(c[M] / cref ** 2) * (w[j] / W) for j in ref_masses}
        g[M] = g.get(M, 0.0) + 1 / cref
        out[M] = (c[M] / cref - 1, math.sqrt(sum(g[j] ** 2 * s[j] ** 2 for j in g)))
    return out, cref


def run_sections(text):
    """per ##RUN section: run -> seed, the seeding lines (the [KE-AUDIT] lines and '[EDMD-ENERGY] 0 state loaded' up to E_tot),
    E_div and E_tot of '[EDMD-ENERGY] 2 end of record', and the [EDMD-HEALTH] line count"""
    out = {}
    for sec in re.split(r"^##RUN ", text, flags=re.M)[1:]:
        h = re.match(r"\S+ \S+ run (\d+) seed (\d+)", sec)
        if not h: continue
        start = re.findall(r"^\[KE-AUDIT\].*$", sec, re.M) + re.findall(r"^\[EDMD-ENERGY\] 0 state loaded.*?E_tot=\S+", sec, re.M)
        m = re.search(r"^\[EDMD-ENERGY\] 2 end of record\s+t=\S+ E_gas=\S+ E_div=(\S+) E_tot=(\S+)", sec, re.M)
        out[int(h.group(1))] = dict(seed=int(h.group(2)), start="\n".join(start), nstart=len(start),
                                    ediv=float(m.group(1)) if m else math.nan, etot=float(m.group(2)) if m else math.nan,
                                    health=len(re.findall(r"^\[EDMD-HEALTH\]", sec, re.M)))
    return out


def hist_bins(nu, df):
    k = np.rint(np.asarray(nu) / df).astype(int)
    return {int(b): int((k == b).sum()) for b in np.unique(k)}


def fd_proxy(acfs, dt, M, k):
    """per trajectory M 2 (C(0) - C(k dt)) / (k dt)^2: M times the mean square of the average velocity over k dt"""
    return M * 2 * (acfs[:, 0].astype(float) - acfs[:, k].astype(float)) / (k * dt) ** 2


# ---------------------------------------------------------------- the self-test (synthetic data with known answers)
def selftest():
    rng = np.random.default_rng(20261011); ok = {}
    # 1. the 4a fit recovers a known damped cosine with a slow mode (the M = 50 sampling: dt 0.4333, 20 periods)
    dt, nu, tr = 0.4333, 0.0684, 60.0; t = np.arange(int(20 / (nu * dt)) + 1) * dt
    a = 0.3 * np.exp(-t / 30.0) + 0.7 * np.exp(-t / tr) * np.cos(2 * math.pi * nu * t + 0.1)
    p = FU.fit_acf(a + 1e-4 * rng.standard_normal(len(t)), dt, nu * 1.01)
    ok["fit_acf recovers nu_d (0.1 %) and tau_r (2 %) of a known ACF"] = abs(p[4] / (2 * math.pi) / nu - 1) < 1e-3 and abs(p[3] / tr - 1) < 0.02
    # 2. the jackknife: of a mean, ~ SD/sqrt(n) (40 groups); a failed replicate is counted, not used
    x = rng.standard_normal(400); full, se, nf = jackknife(lambda i: x[i].mean(), 400)
    ok["jackknife SE of a mean ~ SD/sqrt(n) (25 %)"] = abs(full[0] - x.mean()) < 1e-12 and abs(se[0] / (x.std(ddof=1) / 20) - 1) < 0.25 and nf == 0
    _, se1, nf = jackknife(lambda i: math.nan if 0 not in i else x[i].mean(), 400)     # only the replicate without group 0 fails
    ok["jackknife: a replicate that fails is counted and left out"] = nf == 1 and np.isfinite(se1[0])
    # 3. delta: c_50 1 % below a flat reference of three masses; the SE of a mass outside the reference
    c = {50: 0.99, 300: 1.0, 1500: 1.0, 2000: 1.0}; s = {50: 0.002, 300: 0.001, 1500: 0.001, 2000: 0.002}
    d, cref = deltas(c, s, (300, 1500, 2000)); W = 2 / 0.001 ** 2 + 1 / 0.002 ** 2
    ok["delta_50 = -1 % against a flat reference; SE by the delta method"] = abs(cref - 1) < 1e-12 and abs(d[50][0] + 0.01) < 1e-12 \
        and abs(d[50][1] - math.sqrt(0.002 ** 2 + 0.99 ** 2 / W)) < 1e-12
    # 4. the bin histogram on a known grid
    df = 3.4e-4; nus = np.array([213, 213, 214, 212, 213]) * df
    ok["argmax bins counted on their grid"] = hist_bins(nus, df) == {212: 1, 213: 3, 214: 1}
    # 5. the log parser: the seeding lines without the engine's resched tag, the end-of-record energies, the health lines
    sec = ("##RUN 2026-10-11 01:47:55 run {r} seed {s} node cn (3 s)\n[KE-AUDIT] 1 after velocity draw        N=100  KE_tot=91.4\n"
           "[EDMD-ENERGY] 0 state loaded (hold)      t=0 E_gas=99.9 E_div=0 E_tot=99.9 resched={e}\n{h}"
           "[EDMD-ENERGY] 2 end of record            t=1 E_gas=99 E_div={d} E_tot=100 resched={e}\n")
    a2 = run_sections(sec.format(r=7, s=42, e="minimal", h="[EDMD-HEALTH] forced_advance=1\n", d=0.5))
    a3 = run_sections(sec.format(r=7, s=42, e="gen3", h="", d=1.25))
    ok["log parser: seeding lines equal across engines, energies, health lines"] = a2[7]["start"] == a3[7]["start"] and a2[7]["nstart"] == 2 \
        and (a2[7]["seed"], a2[7]["ediv"], a2[7]["etot"], a2[7]["health"], a3[7]["ediv"], a3[7]["health"]) == (42, 0.5, 100.0, 1, 1.25, 0)
    # 6. the finite-difference proxy on an exact oscillation (reduce_B's ACF): M <(dx/dt)^2> = M A^2 omega^2/2 x 2(1 - cos w dt)/(w dt)^2
    A, om, Mt = 0.3, 2 * math.pi * 0.0684, 50.0; tt = np.arange(200000) * dt; y = A * np.cos(om * tt + 0.4); y -= y.mean()
    f = np.fft.rfft(y, 2 * len(y)); ac = np.fft.irfft(f * np.conj(f))[:4] / np.arange(len(y), len(y) - 4, -1)
    want = Mt * A * A * om * om / 2 * 2 * (1 - math.cos(om * dt)) / (om * dt) ** 2
    ok["finite-difference proxy equals its exact value on an oscillation (1e-3)"] = abs(fd_proxy(ac[None, :], dt, Mt, 1)[0] / want - 1) < 1e-3
    # 7. equipartition at fixed energy, not circular: canonical momenta with unequal masses (one heavy of 50, 2N = 4 light) conditioned
    #    on the total energy (the conditional law is the microcanonical one); <M V^2> = 2E/f with f = 2N + 1 = 5, not 2E/(f - 1)
    ms = np.array([50.0, 1, 1, 1, 1]); E0 = 2.5; acc = []
    for _ in range(20):
        pp = rng.standard_normal((200000, 5)) * np.sqrt(ms); E = (pp ** 2 / (2 * ms)).sum(axis=1)
        acc.append(pp[np.abs(E - E0) < 0.025, 0] ** 2 / ms[0])
    acc = np.concatenate(acc); r_ = acc.mean() / (2 * E0 / 5)
    ok[f"<M V^2> = 2E/(2N+1) at fixed E (conditioned canonical sample, {len(acc)} points, 4 %; 2E/(2N) is 25 % off)"] = abs(r_ - 1) < 0.04
    for k, v in ok.items(): print(f"| {k} | {'PASS' if v else '**FAIL**'} |")
    print(f"\nSELF-TEST: {'PASS' if all(ok.values()) else 'FAIL'} ({sum(ok.values())} of {len(ok)})")
    return all(ok.values())


# ---------------------------------------------------------------- the tables
def load(D):
    out = {}
    for e in ENG:
        for M in MASSES:
            d = os.path.join(D, "testG", e, CELL, f"m_{M}")
            r = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
            if r["dt"].nunique() != 1 or r["n"].nunique() != 1: sys.exit(f"STOP: more than one dt or n in {d}")
            z = np.load(os.path.join(d, "acf_runs.npz")); nl = min(len(z[f"run{k}"]) for k in r["run"])
            acfs = np.array([z[f"run{k}"][:nl] for k in r["run"]], dtype=float)
            out[(e, M)] = dict(r=r, acf=acfs, dt=float(r["dt"].iloc[0]), n=int(r["n"].iloc[0]),
                               log=run_sections(open(os.path.join(d, "run.log"), errors="replace").read()))
    return out


def main():
    if "--selftest" in sys.argv: sys.exit(0 if selftest() else 1)
    D = next((a for a in sys.argv[1:] if not a.startswith("--")), os.path.join(HS, "experiments_gen3_gate_koa_261009"))
    print("# POST HOC, INFORMATION: Test G on KOA, the tables of decision 27 task B (261012 sec. 4.7.56). They change no verdict:\n"
          "# TEST G (KOA): FAIL (DIFFERENCE) at M = 50 stands.\n")
    print("## Self-test first\n\n| check | result |\n|---|---|")
    if not selftest(): sys.exit("STOP: the self-test failed")
    X = load(D)
    nu_ref = {M: 0.5 * (X[("gen2", M)]["r"]["nu"].mean() + X[("gen3", M)]["r"]["nu"].mean()) for M in MASSES}
    S, ND = {}, {}
    for (e, M), v in X.items():
        ND[(e, M)] = nu_d_each(v["acf"], v["dt"], nu_ref[M])
        S[(e, M)] = block_stats(v["r"]["nu"].to_numpy(float), ND[(e, M)], v["acf"], v["dt"], nu_ref[M])
    x = {M: T.x_of(M, L0) for M in MASSES}
    # ---- B1
    print("\n## B1 (POST HOC, INFORMATION): per engine and mass, the frequency estimators with jackknife SEs\n")
    print("| engine | M | trajectories | nu_d fit failures | jackknife replicates failed | mean nu_arg | mean nu_d | tau_r | Q | nu_0 | "
          "-1/(4Q^2) [%] |\n|---|---|---|---|---|---|---|---|---|---|---|")
    for e in ENG:
        for M in MASSES:
            f, s, nf = S[(e, M)]
            print(f"| {e} | {M} | {len(X[(e, M)]['r'])} | {int((~np.isfinite(ND[(e, M)])).sum())} | {nf} | {f[0]:.7f} +- {s[0]:.7f} | "
                  f"{f[1]:.7f} +- {s[1]:.7f} | {f[2]:.1f} +- {s[2]:.1f} | {f[3]:.2f} +- {s[3]:.2f} | {f[4]:.7f} +- {s[4]:.7f} | {-100 / (4 * f[3] ** 2):+.3f} |")
    print("\nsec. 16 (Test T legacy, 4a model, averaged ACF) at M = 50: Q 9.3, -1/(4Q^2) = -0.288 %.")
    print("x_M = K(alpha)/(2 pi L_eff) (tests_20260913.x_of, L0 = 10): "
          + "; ".join(f"M {M}: alpha {M / (2 * NS):g}, K {T.k_root(M / (2 * NS)):.5f}, x {x[M]:.7f}" for M in MASSES))
    print("\n### c_k = nu_k / x_M; delta_M = c_M / c_ref - 1, c_ref the inverse-variance weighted mean of c_M over M >= 300 of the same "
          "engine and estimator\n")
    print("| estimator | engine | c_ref | delta_50 [%] | delta_300 [%] | delta_1500 [%] | delta_2000 [%] |\n|---|---|---|---|---|---|---|")
    D50 = {}
    for k, lab in ((0, "nu_arg"), (1, "nu_d"), (4, "nu_0")):
        for e in ENG:
            c = {M: S[(e, M)][0][k] / x[M] for M in MASSES}; s = {M: S[(e, M)][1][k] / x[M] for M in MASSES}
            d, cref = deltas(c, s, (300, 1500, 2000)); D50[(lab, e)] = d[50]
            print(f"| {lab} | {e} | {cref:.5f} | " + " | ".join(f"{100 * d[M][0]:+.3f} +- {100 * d[M][1]:.3f}" for M in MASSES) + " |")
    print(f"\n### delta_50 per engine and gen3 - gen2; methods sec. 16: {SEC16[0]:+.3f} +- {SEC16[1]:.3f} % with nu_0 (Test T legacy, gen2, "
          "against the nine-mass weighted slope)\n")
    print("| estimator | gen2 delta_50 [%] | z vs 0 | z vs sec. 16 | gen3 delta_50 [%] | z vs 0 | z vs sec. 16 | gen3 - gen2 [%] | z |\n"
          "|---|---|---|---|---|---|---|---|---|")
    for lab in ("nu_arg", "nu_d", "nu_0"):
        row = []
        for e in ENG:
            v, sv = 100 * D50[(lab, e)][0], 100 * D50[(lab, e)][1]
            row += [f"{v:+.3f} +- {sv:.3f}", f"{v / sv:+.2f}", f"{(v - SEC16[0]) / math.hypot(sv, SEC16[1]):+.2f}"]
        a, b = D50[(lab, "gen2")], D50[(lab, "gen3")]; dd, sd = b[0] - a[0], math.hypot(a[1], b[1])
        print(f"| {lab} | " + " | ".join(row) + f" | {100 * dd:+.3f} +- {100 * sd:.3f} | {dd / sd:+.2f} |")
    # ---- B2
    print("\n## B2 (POST HOC, INFORMATION): the argmax bins, nu_arg = k df with df = 1/(n dt), per mass in both arms\n")
    for M in MASSES:
        v2, v3 = X[("gen2", M)], X[("gen3", M)]
        if (v2["n"], v2["dt"]) != (v3["n"], v3["dt"]): print(f"M = {M}: n or dt differ between the arms")
        df = 1 / (v2["n"] * v2["dt"]); h2, h3 = hist_bins(v2["r"]["nu"], df), hist_bins(v3["r"]["nu"], df)
        k2, k3 = np.rint(v2["r"]["nu"] / df), np.rint(v3["r"]["nu"] / df)
        off = max(abs(v2["r"]["nu"] / df - k2).max(), abs(v3["r"]["nu"] / df - k3).max())
        print(f"M = {M}: n {v2['n']}, dt {v2['dt']:.7f}, df = {df:.5e} = {100 * df / nu_ref[M]:.3f} % of nu; largest distance of a nu_arg "
              f"from its bin {off:.1e} bins; SD in bins: gen2 {k2.std(ddof=1):.2f}, gen3 {k3.std(ddof=1):.2f}; mean gen3 - gen2 "
              f"{k3.mean() - k2.mean():+.3f} bins = {100 * (k3.mean() - k2.mean()) * df / nu_ref[M]:+.3f} % of nu\n")
        print("| bin k | nu = k df | gen2 | gen3 |\n|---|---|---|---|")
        for b in sorted(set(h2) | set(h3)): print(f"| {b} | {b * df:.7f} | {h2.get(b, 0)} | {h3.get(b, 0)} |")
        print(f"\n(the two arms' counts agree bin by bin: {'yes' if h2 == h3 else 'no'})\n")
    # ---- B3
    print("## B3 (POST HOC, INFORMATION): pairing by start (the same exact seed, the same initial state), gen3 - gen2 per seed\n")
    print("| M | common seeds | same run index | seeding lines byte-identical | nu_arg: mean diff [% of nu] | SE | same-seed corr. | "
          "nu_d: mean diff [% of nu] | SE | same-seed corr. | unpaired SE, nu_arg [%] | unpaired SE, nu_d [%] |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
    for M in MASSES:
        v2, v3 = X[("gen2", M)], X[("gen3", M)]
        r2 = v2["r"].assign(nd=ND[("gen2", M)]); r3 = v3["r"].assign(nd=ND[("gen3", M)])
        m = r2.merge(r3, on="seed", suffixes=("_2", "_3"))
        same_run = int((m["run_2"] == m["run_3"]).sum())
        ident = sum(1 for rr2, rr3 in zip(m["run_2"], m["run_3"]) if v2["log"].get(rr2, {}).get("nstart", 0) > 0
                    and v2["log"][rr2]["start"] == v3["log"].get(rr3, {}).get("start"))
        row = []
        for col in ("nu", "nd"):
            mm = m[np.isfinite(m[f"{col}_2"]) & np.isfinite(m[f"{col}_3"])]; d = 100 * (mm[f"{col}_3"] - mm[f"{col}_2"]) / nu_ref[M]
            row += [f"{d.mean():+.3f}", f"{d.std(ddof=1) / math.sqrt(len(d)):.3f}", f"{np.corrcoef(mm[f'{col}_2'], mm[f'{col}_3'])[0, 1]:+.3f}"]
        se_u = [100 * math.hypot(np.nanstd(a, ddof=1) / math.sqrt(np.isfinite(a).sum()), np.nanstd(b, ddof=1) / math.sqrt(np.isfinite(b).sum()))
                / nu_ref[M] for a, b in ((r2["nu"], r3["nu"]), (r2["nd"], r3["nd"]))]
        print(f"| {M} | {len(m)} | {same_run} | {ident} | " + " | ".join(row) + f" | {se_u[0]:.3f} | {se_u[1]:.3f} |")
    # ---- B4
    print(f"\n## B4 (POST HOC, INFORMATION): the divider's kinetic temperature; equipartition at fixed E over f = 2N + 1 = {F_DOF}: "
          f"<M V^2> = 2 E_tot / {F_DOF}\n")
    print("### (i) exact, one instant per trajectory: [EDMD-ENERGY] 2 end of record (HD_KE_TRACE=1)\n")
    print("| engine | M | trajectories | mean M V^2 = 2 E_div | SE | expected 2 E_tot / 201 | ratio | (mean - expected) / SE |\n|---|---|---|---|---|---|---|---|")
    pool = {e: ([], []) for e in ENG}
    for e in ENG:
        for M in MASSES:
            L = [v for v in X[(e, M)]["log"].values() if np.isfinite(v["ediv"])]
            mv = np.array([2 * v["ediv"] for v in L]); ex = float(np.mean([2 * v["etot"] / F_DOF for v in L]))
            pool[e][0].extend(mv); pool[e][1].append(ex); se = mv.std(ddof=1) / math.sqrt(len(mv))
            print(f"| {e} | {M} | {len(mv)} | {mv.mean():.4f} | {se:.4f} | {ex:.6f} | {mv.mean() / ex:.4f} | {(mv.mean() - ex) / se:+.2f} |")
        mv = np.array(pool[e][0]); ex = float(np.mean(pool[e][1])); se = mv.std(ddof=1) / math.sqrt(len(mv))
        print(f"| {e} | all four | {len(mv)} | {mv.mean():.4f} | {se:.4f} | {ex:.6f} | {mv.mean() / ex:.4f} | {(mv.mean() - ex) / se:+.2f} |")
    print("\n### (ii) over the record, a proxy: M <(dx/dt)^2> over one and two sample intervals, from the fetched ACF (same estimator in both arms)\n")
    print("| M | dt | omega dt | gen2, k = 1 | gen3, k = 1 | gen3/gen2 - 1, k = 1 [%] | gen2, k = 2 | gen3, k = 2 | gen3/gen2 - 1, k = 2 [%] |\n"
          "|---|---|---|---|---|---|---|---|---|")
    for M in MASSES:
        row = []; dtM = X[("gen2", M)]["dt"]
        for k in (1, 2):
            a = fd_proxy(X[("gen2", M)]["acf"], dtM, M, k); b = fd_proxy(X[("gen3", M)]["acf"], dtM, M, k)
            sa, sb = a.std(ddof=1) / math.sqrt(len(a)), b.std(ddof=1) / math.sqrt(len(b))
            q = b.mean() / a.mean(); sq = q * math.hypot(sa / a.mean(), sb / b.mean())
            row += [f"{a.mean():.4f} +- {sa:.4f}", f"{b.mean():.4f} +- {sb:.4f}", f"{100 * (q - 1):+.2f} +- {100 * sq:.2f}"]
        print(f"| {M} | {dtM:.4f} | {2 * math.pi * nu_ref[M] * dtM:.3f} | " + " | ".join(row) + " |")
    print(f"\n(expected <M V^2> = 2 x 100 / {F_DOF} = {200 / F_DOF:.5f}; the proxy is biased low by the decorrelation within k dt)")
    # ---- B5
    print("\n## B5 (POST HOC, INFORMATION): [EDMD-HEALTH] lines in the fetched logs (printed only when a gen2 counter is nonzero)\n")
    print("| engine | M | ##RUN sections | [EDMD-HEALTH] lines |\n|---|---|---|---|")
    for e in ENG:
        for M in MASSES:
            L = X[(e, M)]["log"]; print(f"| {e} | {M} | {len(L)} | {sum(v['health'] for v in L.values())} |")


if __name__ == "__main__":
    main()
