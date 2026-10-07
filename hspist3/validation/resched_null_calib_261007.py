#!/usr/bin/env python3
"""##CHRIS 2026-10-07 (261012 sec. 4.4.10, plan-author decision of 2026-10-07, parts A3 and B): analysis only, no simulation,
on the fetched replay of build 73fc07f (hspist3/experiments_resched_gate_261005/) and the 279282b campaign data of the same cells.

A3  numbers for the corrections to sec. 4.4.9: the size of one quantization step in a 25-seed mean; the light-mass point
    estimates; whether the smoke pilot's seeds are among the campaign's (and the pilot reported next to the replay); the power
    of a 25-seed confirmatory test; the arguments for noise (look-elsewhere numbers, opposite signs, no trend with mass) and
    against (the per-mass chi2).
B1  permutation test, both replayed method-B cells: within each mass the 25 old and 25 new per-seed nu are relabelled at random
    (1e5 relabelings, fixed RNG seed); per relabeling the per-mass z, their chi2 over the nine masses, max |z|, and c_s and
    k_S^dyn with the registered estimators (paper1_populate_cs_err_20261002.slope_with_errors: unweighted through-origin slope;
    paper1_confinement_results_261004.identity: inverse-variance mean of (M + 2N_s/3) omega^2 / 2 over alpha >= 5), checked
    against the registered values first. Permutation p next to the nominal p.
B2  the chi2 of the registered single-c_s through-origin fit over the nine masses, old and new, both cells, with the per-mass
    residuals in SE units.
B3  the sorted per-seed nu, old and new, at pi/8 for alpha = 0.5, 5 and 15.
B4  the contact audit resolved by class over the 565 replay trajectories: per cell and class the maximum and the trajectory
    that holds it.
usage (from hspist3/):  python3 validation/resched_null_calib_261007.py
"""
import contextlib, glob, io, math, os, re, sys
import numpy as np, pandas as pd
from scipy.stats import chi2 as CHI2, norm, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import paper1_confinement_results_261004 as R
import tests_20260913 as T

NEW = os.path.join(HS, "experiments_resched_gate_261005")
CELLS = ("e0p10_H_H10_L39.25", "epi8_H_H10_L10")
NPERM, RNG_SEED = 100000, 20261007
SMOKE_NU = {50: 0.07040729, 100: 0.05624957, 200: 0.04478494, 300: 0.03664223, 500: 0.02967332, 750: 0.02414966,
            1000: 0.02112673, 1500: 0.01718229, 2000: 0.01506197}   # [DATA] KOA smoke log, build 73fc07f (261012 sec. 4.4.7)
CONTACT = re.compile(r"\[EDMD-CONTACT\] executed events (\d+); max abs\(contact distance\) \[px\]: disk-disk (\S+), "
                     r"outer walls (\S+), divider (\S+), pistons (\S+)\n")
CLASSES = ("disk-disk", "outer walls", "divider", "pistons")


@contextlib.contextmanager
def new_tree():
    old = (R.HS, R.ALLOW_NEW_BUILD); R.HS = NEW; R.ALLOW_NEW_BUILD = True
    try: yield
    finally: R.HS, R.ALLOW_NEW_BUILD = old


def load():
    """Registered estimators on old and new data, plus the per-seed nu arrays."""
    with contextlib.redirect_stdout(io.StringIO()):
        CS = {c["cid"]: c for c in R.cells()}
    out = {}
    for cid in CELLS:
        c_old = dict(CS[cid]); c_new = dict(CS[cid])
        with contextlib.redirect_stdout(io.StringIO()):
            R.method_B(c_old); R.method_A(c_old); R.identity(c_old)
            for k in ("static", "s_kT", "kT", "F2term"): c_new[k] = c_old[k]
            with new_tree(): R.method_B(c_new); R.identity(c_new)
        out[cid] = dict(c=CS[cid], old=c_old, new=c_new,
                        nu_old={r["M"]: np.asarray(r["nus"], float) for r in c_old["B"]},
                        nu_new={r["M"]: np.asarray(r["nus"], float) for r in c_new["B"]})
    return out


def estimators(means, ses, Ms, x, Ns):
    """Vectorised registered estimators. means, ses: (..., 9). Returns c_s (unweighted through-origin slope) and k_S^dyn."""
    cs = (means * x).sum(-1) / (x * x).sum()
    heavy = np.array([any(abs(M / (2.0 * Ns) - a) < 1e-9 for a in R.HEAVY) for M in Ms])   # as R.identity: alpha in HEAVY
    Mh = np.array([M + 2.0 * Ns / 3.0 for M in Ms])
    om = 2 * math.pi * means; k = Mh * om * om / 2.0; s = k * 2 * ses / means; w = 1.0 / s ** 2
    w = np.where(heavy, w, 0.0)
    return cs, (w * k).sum(-1) / w.sum(-1)


def section_a3(D):
    print("## A3 -- numbers for the corrections to sec. 4.4.9\n")
    print("### Size of one quantization step in a 25-seed mean\n")
    print("| cell | M | alpha | bin df/nu [%] | one seed moving one bin: change of the 25-seed mean [%] | per-mass SE of the old mean [%] | ratio |")
    print("|---|---|---|---|---|---|---|")
    for cid in CELLS:
        for M in sorted(D[cid]["nu_old"]):
            o, n = D[cid]["nu_old"][M], D[cid]["nu_new"][M]; u = np.unique(np.r_[o, n])
            df = float(np.median(np.diff(u))); step = 100 * df / o.mean() / len(o); se = 100 * o.std(ddof=1) / math.sqrt(len(o)) / o.mean()
            print(f"| {cid} | {M} | {M / (2 * D[cid]['c']['Ns']):g} | {100 * df / o.mean():.3f} | {step:.4f} | {se:.3f} | {step / se:.2f} |")
    print("\n### Light masses: point estimates (new - old)/old with SE, and the smoke pilot's shift\n")
    print("| cell | alpha | (new - old)/old [%] | SE [%] | 95 % interval [%] | smoke-pilot shift vs campaign mean [%] | its distance from the replay estimate [SE] |")
    print("|---|---|---|---|---|---|---|")
    for cid in CELLS:
        for M in (50, 100):
            o, n = D[cid]["nu_old"][M], D[cid]["nu_new"][M]
            d = (n.mean() - o.mean()) / o.mean(); s = math.hypot(o.std(ddof=1) / 5, n.std(ddof=1) / 5) / o.mean()
            sm = (SMOKE_NU[M] - o.mean()) / o.mean() if cid == "epi8_H_H10_L10" else float("nan")
            print(f"| {cid} | {M / (2 * D[cid]['c']['Ns']):g} | {100 * d:+.2f} | {100 * s:.2f} | [{100 * (d - 1.96 * s):+.2f}, {100 * (d + 1.96 * s):+.2f}] | "
                  f"{100 * sm:+.2f} | {(sm - d) / s if math.isfinite(sm) else float('nan'):+.2f} |")
    print("\n### The smoke pilot's seeds against the campaign's (epi8_H_H10_L10)\n")
    camp = {}
    for l in open(os.path.join(HS, "cluster", "confinement_20261013", "tasks_B_epi8_H_H10_L10.txt")):
        f = l.split(); camp.setdefault(int(f[2]), set()).add(int(f[4]))
    pilot = {M: T.run_seed(20261013, 0, mi, 0) for mi, M in enumerate(T.A1_MASSES)}
    print("| M | smoke pilot seed (run_seed(20261013, 0, m, 0)) | among the campaign's 25 seeds of this mass | among all 225 |\n|---|---|---|---|")
    allc = set().union(*camp.values())
    for M in T.A1_MASSES:
        print(f"| {M} | {pilot[M]} | {'yes' if pilot[M] in camp[M] else 'no'} | {'yes' if pilot[M] in allc else 'no'} |")
    print("\n### The smoke pilot reported next to the replay (one more new-engine trajectory per mass)\n")
    print("| M | replay new mean nu (25) | SE | smoke pilot nu | pilot - replay mean [SD of one seed] | new mean with the pilot (26) | SE |")
    print("|---|---|---|---|---|---|---|")
    for M in T.A1_MASSES:
        n = D["epi8_H_H10_L10"]["nu_new"][M]; n26 = np.r_[n, SMOKE_NU[M]]
        print(f"| {M} | {n.mean():.6f} | {n.std(ddof=1) / 5:.6f} | {SMOKE_NU[M]:.6f} | {(SMOKE_NU[M] - n.mean()) / n.std(ddof=1):+.2f} | "
              f"{n26.mean():.6f} | {n26.std(ddof=1) / math.sqrt(26):.6f} |")
    print("\n### Power of a 25-seed confirmatory test (sec. 4.4.9 option 2) if the true k_S^dyn shift equals the observed one\n")
    dz = 0.0356 / 0.0155   # observed difference and sigma_diff of G-E4 number 5 (sec. 4.4.9)
    p3 = norm.sf(3 - dz) + norm.cdf(-3 - dz); p2 = norm.cdf(2 - dz) - norm.cdf(-2 - dz)
    print(f"expected z = 0.0356 / 0.0155 = {dz:.3f}; P(|z| >= 3) = {p3:.3f}; P(|z| < 2) = {p2:.3f}; P(2 <= |z| < 3) = {1 - p3 - p2:.3f}")
    print("\n### Arguments for noise and against, each computed here\n")
    p1 = 2 * norm.sf(2)
    print(f"P(at least one of 9 independent |z| >= 2) = 1 - (1 - {p1:.4f})^9 = {1 - (1 - p1) ** 9:.3f}")
    q = 2 * norm.sf(2.81)
    for n in (18, 27): print(f"P(max |z| >= 2.81 among {n} independent numbers) = 1 - (1 - {q:.5f})^{n} = {1 - (1 - q) ** n:.3f}")
    for cid in CELLS:
        o, n = D[cid]["old"], D[cid]["new"]
        print(f"{cid}: c_s {100 * (n['cs'] / o['cs'] - 1):+.2f} %, k_S^dyn {100 * (n['kS'] / o['kS'] - 1):+.2f} % "
              f"({'opposite' if (n['cs'] - o['cs']) * (n['kS'] - o['kS']) < 0 else 'same'} directions)")
        zs = []; signs = ""
        for M in sorted(D[cid]["nu_old"]):
            a, b = D[cid]["nu_old"][M], D[cid]["nu_new"][M]
            z = (b.mean() - a.mean()) / math.hypot(a.std(ddof=1) / 5, b.std(ddof=1) / 5); zs.append(z); signs += "+" if z > 0 else "-"
        rho, pr = spearmanr(sorted(D[cid]["nu_old"]), zs)
        print(f"  signs of the nine per-mass differences, M ascending: {signs}; Spearman rho(M, z) = {rho:+.2f} (p = {pr:.2f}); "
              f"against: chi2 = {np.sum(np.square(zs)):.2f} / 9, nominal p = {CHI2.sf(np.sum(np.square(zs)), 9):.4f}")


def section_b1(D):
    print(f"\n## B1 -- permutation test ({NPERM} relabelings per cell, numpy default_rng({RNG_SEED}))\n")
    rng = np.random.default_rng(RNG_SEED)
    print("| cell | statistic | observed | nominal p | permutation p |\n|---|---|---|---|---|")
    for cid in CELLS:
        d = D[cid]; Ms = sorted(d["nu_old"]); Ns = d["c"]["Ns"]; x = np.asarray(d["old"]["x"], float)
        mo = np.array([d["nu_old"][M].mean() for M in Ms]); so = np.array([d["nu_old"][M].std(ddof=1) / math.sqrt(len(d["nu_old"][M])) for M in Ms])
        mn = np.array([d["nu_new"][M].mean() for M in Ms]); sn = np.array([d["nu_new"][M].std(ddof=1) / math.sqrt(len(d["nu_new"][M])) for M in Ms])
        cso, kso = estimators(mo, so, Ms, x, Ns); csn, ksn = estimators(mn, sn, Ms, x, Ns)
        chk = [abs(cso - d["old"]["cs"]), abs(csn - d["new"]["cs"]), abs(kso - d["old"]["kS"]) / kso, abs(ksn - d["new"]["kS"]) / ksn]
        if max(chk) > 1e-12: sys.exit(f"STOP: vectorised estimators do not reproduce the registered ones for {cid}: {chk}")
        z = (mn - mo) / np.hypot(so, sn); chi_obs = float((z * z).sum()); mz_obs = float(np.abs(z).max())
        dcs, dks = csn - cso, ksn - kso
        zcs = dcs / math.hypot(d["old"]["cs_err"], d["new"]["cs_err"]); zks = dks / math.hypot(d["old"]["s_kS"], d["new"]["s_kS"])
        PMo = np.empty((NPERM, 9)); PSo = np.empty((NPERM, 9)); PMn = np.empty((NPERM, 9)); PSn = np.empty((NPERM, 9))
        for j, M in enumerate(Ms):
            pool = np.r_[d["nu_old"][M], d["nu_new"][M]]; na = len(d["nu_old"][M])
            idx = np.argsort(rng.random((NPERM, len(pool))), axis=1); v = pool[idx]
            a, b = v[:, :na], v[:, na:]
            PMo[:, j], PSo[:, j] = a.mean(1), a.std(1, ddof=1) / math.sqrt(na)
            PMn[:, j], PSn[:, j] = b.mean(1), b.std(1, ddof=1) / math.sqrt(len(pool) - na)
        pz = (PMn - PMo) / np.hypot(PSo, PSn); pchi = (pz * pz).sum(1); pmz = np.abs(pz).max(1)
        pcso, pkso = estimators(PMo, PSo, Ms, x, Ns); pcsn, pksn = estimators(PMn, PSn, Ms, x, Ns)
        rows = [("chi2 of the nine per-mass z", chi_obs, CHI2.sf(chi_obs, 9), (pchi >= chi_obs).mean()),
                ("max |z| over nine masses", mz_obs, 1 - (1 - 2 * norm.sf(mz_obs)) ** 9, (pmz >= mz_obs).mean()),
                ("k_S^dyn new - old (z)", f"{dks:+.6g} ({zks:+.2f})", 2 * norm.sf(abs(zks)), (np.abs(pksn - pkso) >= abs(dks)).mean()),
                ("c_s new - old (z)", f"{dcs:+.6g} ({zcs:+.2f})", 2 * norm.sf(abs(zcs)), (np.abs(pcsn - pcso) >= abs(dcs)).mean())]
        for name, obs, pn, pp in rows:
            print(f"| {cid} | {name} | {obs if isinstance(obs, str) else f'{obs:.3f}'} | {pn:.4f} | {pp:.4f} |")
    print("\n(estimator check: the vectorised c_s and k_S^dyn reproduce the registered values of both cells, old and new, to 1e-12)")


def section_b2(D):
    from paper1_populate_cs_err_20261002 import slope_with_errors
    print("\n## B2 -- the registered single-c_s through-origin fit across the nine masses (slope_with_errors; dof = 8)\n")
    print("| cell | data | c_s | chi2 (= chi2_red x 8) | chi2_red | per-mass residuals (y - c_s x)/SE, M ascending |\n|---|---|---|---|---|---|")
    for cid in CELLS:
        for lab in ("old", "new"):
            c = D[cid][lab]; x = np.asarray(c["x"], float); y = np.array([r["nu"] for r in c["B"]]); sy = np.array([r["se"] for r in c["B"]])
            s, err, errs, chi2r = slope_with_errors(x, y, sy); res = (y - s * x) / sy
            print(f"| {cid} | {lab} | {s:.5f} | {chi2r * 8:.1f} | {chi2r:.2f} | {' '.join(f'{v:+.1f}' for v in res)} |")


def section_b3(D):
    print("\n## B3 -- sorted per-seed nu at pi/8 (epi8_H_H10_L10), old and new\n")
    for M in (50, 500, 1500):
        for lab, key in (("old", "nu_old"), ("new", "nu_new")):
            v = np.sort(D["epi8_H_H10_L10"][key][M])
            print(f"alpha = {M / 100:g}, {lab} (n = {len(v)}, mean {v.mean():.6f}): " + " ".join(f"{u:.6f}" for u in v))
        print()


def section_b4():
    print("## B4 -- contact audit by class over the replay trajectories (max abs(contact distance) at executed events, px)\n")
    print("| cell | class | maximum | held by | trajectories with an audit line |\n|---|---|---|---|---|")
    RB = os.path.join(NEW, R.REL_B)
    for cid in CELLS:
        best = {k: (-1.0, "") for k in CLASSES}; n = 0
        for d in sorted(glob.glob(os.path.join(RB, cid, "m_*"))):
            for sec in open(os.path.join(d, "run.log"), errors="ignore").read().split("##RUN")[1:]:
                hdr = sec.split("\n", 1)[0].strip(); m = CONTACT.search(sec)
                if not m: continue
                n += 1
                for k, v in zip(CLASSES, m.groups()[1:]):
                    if float(v) > best[k][0]: best[k] = (float(v), f"{os.path.basename(d)}, {hdr}")
        for k in CLASSES: print(f"| {cid} | {k} | {best[k][0]:.3e} | {best[k][1]} | {n} |")
    A = os.path.join(NEW, "experiments_energy_transfer", "paper1_confinement_Afix_261004", "epi8_H_H10_L10")
    best = {k: (-1.0, "") for k in CLASSES}; n = 0
    for f in sorted(glob.glob(os.path.join(A, "x_*", "run_*.log"))):
        m = CONTACT.search(open(f, errors="ignore").read())
        if not m: continue
        n += 1
        for k, v in zip(CLASSES, m.groups()[1:]):
            if float(v) > best[k][0]: best[k] = (float(v), os.path.relpath(f, A))
    for k in CLASSES: print(f"| AF epi8_H_H10_L10 | {k} | {best[k][0]:.3e} | {best[k][1]} | {n} |")


def main():
    D = load()
    print("# Null calibration of the engine replay (261012 sec. 4.4.10) -- analysis only\n")
    section_a3(D); section_b1(D); section_b2(D); section_b3(D); section_b4()


if __name__ == "__main__":
    main()
