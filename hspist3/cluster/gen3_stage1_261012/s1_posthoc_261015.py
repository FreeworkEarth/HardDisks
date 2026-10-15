#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 55 item 3): POST-HOC INFORMATION on stage 1 Part B, dated 2026-10-15. POST HOC: it never
stands in place of the verdicts of record (decision 55 item 2, the registered script's output). Read only, from the fetched folder, with
the full report's numbers: s1_verdict_261012.analyse(<root>) over all six groups (not --interim, its jackknife) and its own functions
(wls, fits_for, verdict_k); nothing is fitted or estimated in any other way, except where a line says so (the bimodality coefficient).
  (a) Verdict K's statistics on Phi without eta 0.60 (11 densities), and without 0.60 and 0.50 (10): X^2, dof, p, max |z|.
  (b) eta 0.60: per size psi6 (mean +- SE), E1 (c_s +- sigma, masses, chi2/dof) and Delta; the N 100 cell's per-seed psi6: the
      record means (psi6_rec, per compartment psi6_rec_0/1) and the instantaneous value at the hold's end (psi6_hold_end_0/1) of every
      used trajectory: n, mean, SD, min, quartiles, max, a histogram, per mass, and the bimodality coefficient
      BC = (g^2 + 1) / (k + 3 (n-1)^2 / ((n-2)(n-3))) (g, k: sample skewness and excess kurtosis, bias-corrected; BC > 5/9 hints at two
      modes or strong skew; a heuristic, stated as such); the same for N 100 at 0.55 and 0.63 and N 400 at 0.60, for comparison;
      the F1 fit without N 100 (as printed by the verdict script) and the F2 fit on the three sizes 400, 900, 1600 (3 points,
      3 parameters: no degree of freedom; its sigma_a is the propagated point errors only).
  (c) Each density's z^2, its share of Verdict K's X^2.
  (d) eta 0.50: the four Delta_N +- sigma, the F2 curve at each size, the residuals and residual / sigma, their chi2 and p.
  (e) The table for the paper (one CSV): per density, a (F2) +- sigma_a, z, F2 p, a_Z +- its error, and the per-size Delta_N +- sigma.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_posthoc_261015.py --data <the fetched folder> --csv <new CSV path>
"""
import argparse, csv, math, os, sys
import numpy as np, pandas as pd
from scipy.stats import chi2 as CHI2, skew, kurtosis
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s1_verdict_261012 as V

NS = (100, 400, 900, 1600)


def kstats(fits, etas):
    K = V.verdict_k({e: fits[e] for e in etas})
    zmax = max(K["z"].items(), key=lambda kv: abs(kv[1]))
    return (f"{len(etas)} densities: X^2 = {K['X2']:.2f} ({K['dof']} dof), p = {K['p']:.4f}, max |z| = {abs(zmax[1]):.2f} (eta {zmax[0]:.4g}); "
            f"|z| >= 3.5 at: {', '.join(f'{e:.4g}' for e in K['bad']) or 'none'}")


def traj(root, res, C, cid, cols):
    """per used trajectory (the selection of s1_verdict_261012.per_cell: folders not flagged, clean runs of red_nu.csv) the s1_runs columns"""
    c = C[cid]; out = []
    for M in sorted(c["masses"]):
        r = res["mass"].get((cid, M)); d = os.path.join(root, c["g"], "merged", "stage1", "gen3", cid, f"m_{M}")
        p1, pr, ph = (os.path.join(d, f) for f in ("s1_runs.csv", "red_nu.csv", "runs_health.csv"))
        if r is None or r["flagged"] or not os.path.exists(p1) or not os.path.exists(pr): continue
        s1 = pd.read_csv(p1)
        if not len(s1): continue
        ok = set(pd.read_csv(pr)["run"].astype(int))
        if os.path.exists(ph):
            h = pd.read_csv(ph)
            if len(h): ok &= {int(x.run) for x in h.itertuples() if int(x.gen3_clean) == 1}
        s = s1[s1["run"].isin(sorted(ok))]
        for x in s.itertuples(): out.append(dict(M=M, alpha=c["masses"][M], **{k: float(getattr(x, k)) for k in cols}))
    return pd.DataFrame(out)


def describe(v, nb=12):
    v = np.asarray(v, float); n = len(v)
    g = skew(v, bias=False); k = kurtosis(v, fisher=True, bias=False); bc = (g * g + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
    q = np.quantile(v, [0, 0.25, 0.5, 0.75, 1]); h, e = np.histogram(v, bins=nb)
    return (f"n {n}, mean {v.mean():.4f}, SD {v.std(ddof=1):.4f}, min {q[0]:.4f}, quartiles {q[1]:.4f} / {q[2]:.4f} / {q[3]:.4f}, max {q[4]:.4f}; "
            f"skewness {g:+.2f}, excess kurtosis {k:+.2f}, BC {bc:.3f} ({'> 5/9' if bc > 5 / 9 else '<= 5/9'}); "
            f"histogram ({nb} equal bins {e[0]:.3f} .. {e[-1]:.3f}): {' '.join(str(x) for x in h)}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--csv", required=True); a = ap.parse_args()
    if os.path.exists(a.csv): sys.exit(f"REFUSED: {a.csv} exists")
    res = V.analyse(a.data, interim=False, jackknife=True, quiet=True); C = V.registered()
    fc, fz = res["fitsC"], res["fitsZ"]; K = res["K"]
    print("# POST HOC (2026-10-15; decision 55 item 3), stage 1 Part B: information, never in place of the verdicts of record (decision 55 item 2)\n")
    print(f"groups present: {', '.join(res['groups_present'])}; stopped by the stop rule: {res['stopped'] or 'none'}\n")
    # (a)
    phi = list(V.PHI)
    print("## (a) Verdict K's statistics on Phi, with and without eta 0.60 and 0.50 (the registered F2 fits, s1_verdict_261012.verdict_k)\n")
    print(f"- all of Phi (the verdict of record): {kstats(fc, phi)}")
    print(f"- without 0.60: {kstats(fc, [e for e in phi if abs(e - 0.60) > 1e-9])}")
    print(f"- without 0.60 and 0.50: {kstats(fc, [e for e in phi if abs(e - 0.60) > 1e-9 and abs(e - 0.50) > 1e-9])}")
    # (b)
    print("\n## (b) eta 0.60\n\n| N | psi6 +- SE (n) | E1 c_s +- sigma | masses | chi2/dof | Delta [%] | low-Q | hold t | record z |\n|---|---|---|---|---|---|---|---|---|")
    kr = V.kr_cs(0.60)
    for N in NS:
        c = res["cell"][V.GEN.cell_id("B", N, 0.60, "default")]; E = c["E1"]
        print(f"| {N} | {c['psi6'][0]:.4f} +- {c['psi6'][1]:.4f} ({c['psi6'][2]}) | {E['cs']:.5f} +- {E['sig']:.5f} | {E['n']} | {E['chi2']:.1f}/{E['dof']} | "
              f"{100 * (E['cs'] / kr - 1):+.3f} +- {100 * E['sig'] / kr:.3f} | {'FLAG' if E['lowQ'] else '-'} | {c['hold_t'][0]:+.2f} | {c['rec_z'][0]:+.2f} |")
    cols = ["psi6_rec", "psi6_rec_0", "psi6_rec_1", "psi6_hold_end_0", "psi6_hold_end_1"]
    print("\nthe per-seed psi6 of the N 100 cell at eta 0.60 and, for comparison, of N 100 at 0.55 and 0.63 and N 400 at 0.60 (every used trajectory, all masses):\n")
    for N, e in ((100, 0.60), (100, 0.55), (100, 0.63), (400, 0.60)):
        t = traj(a.data, res, C, V.GEN.cell_id("B", N, e, "default"), cols)
        print(f"- N {N}, eta {e:.2f}:")
        print(f"  - the record means psi6_rec (one per trajectory): {describe(t.psi6_rec)}")
        print(f"  - the compartments' record means psi6_rec_0 and psi6_rec_1 (two per trajectory): {describe(np.r_[t.psi6_rec_0, t.psi6_rec_1])}; "
              f"correlation of the two compartments {np.corrcoef(t.psi6_rec_0, t.psi6_rec_1)[0, 1]:+.2f}")
        print(f"  - the instantaneous value at the hold's end, psi6_hold_end_0 and _1 (two per trajectory): {describe(np.r_[t.psi6_hold_end_0, t.psi6_hold_end_1])}")
        gm = t.groupby("alpha").psi6_rec.agg(["mean", "std", "count"])
        print("  - the record means per mass, mean (SE): " + ", ".join(f"alpha {al:g} {m['mean']:.4f} ({m['std'] / math.sqrt(m['count']):.4f})" for al, m in gm.iterrows()))
        if (N, e) == (100, 0.60):
            g = t.groupby("M")
            print("  - the record means per mass (record length 200 predicted periods of each mass):\n\n    | M | alpha | n | mean | SD | min | max |\n    |---|---|---|---|---|---|---|")
            for M, x in g: print(f"    | {M} | {x.alpha.iloc[0]:g} | {len(x)} | {x.psi6_rec.mean():.4f} | {x.psi6_rec.std(ddof=1):.4f} | {x.psi6_rec.min():.4f} | {x.psi6_rec.max():.4f} |")
            print()
    f = fc[0.60]; F1x = f["F1x"]
    print(f"\nthe F1 fit without N 100 at eta 0.60 (as the verdict script prints it): a = {100 * F1x['beta'][0]:+.3f} +- {100 * F1x['se'][0]:.3f} %, "
          f"b = {100 * F1x['beta'][1]:+.2f} %, chi2 = {F1x['chi2']:.2f} ({F1x['dof']} dof), p = {F1x['p']:.3f}")
    pts = []
    for N in (400, 900, 1600):
        E = res["cell"][V.GEN.cell_id("B", N, 0.60, "default")]["E1"]; pts.append((N, E["cs"] / kr - 1, E["sig"] / kr))
    u = np.array([p[0] ** -0.5 for p in pts]); y = np.array([p[1] for p in pts]); s = np.array([p[2] for p in pts])
    W = V.wls(u, y, s, 2)
    print(f"the F2 fit on the three sizes 400, 900, 1600 at eta 0.60 (3 points, 3 parameters, NO degree of freedom: the curve passes through the "
          f"three points; sigma_a is the propagated point errors only, no goodness of fit exists): a = {100 * W['beta'][0]:+.3f} +- {100 * W['se'][0]:.3f} %, "
          f"b = {100 * W['beta'][1]:+.2f} %, d = {100 * W['beta'][2]:+.2f} %, chi2 = {W['chi2']:.2e} ({W['dof']} dof)")
    # (c)
    print("\n## (c) Each density's share of Verdict K's X^2 (z = a / sigma_a of F2)\n\n| eta | z | z^2 | share of X^2 |\n|---|---|---|---|")
    for e in phi:
        z = K["z"][e]; print(f"| {e:.4g} | {z:+.2f} | {z * z:.2f} | {100 * z * z / K['X2']:.1f} % |")
    print(f"| all | | {K['X2']:.2f} | 100 % |")
    # (d)
    f = fc[0.50]; F2 = f["F2"]; b = F2["beta"]; kr = V.kr_cs(0.50)
    print("\n## (d) eta 0.50: the four Delta_N and the F2 residuals\n\n| N | u = N^(-1/2) | Delta [%] | sigma [%] | F2 at u [%] | residual [%] | residual / sigma |\n|---|---|---|---|---|---|---|")
    rr = []
    for N in NS:
        E = res["cell"][V.GEN.cell_id("B", N, 0.50, "default")]["E1"]; uu = N ** -0.5; d = E["cs"] / kr - 1; sg = E["sig"] / kr
        fv = b[0] + b[1] * uu + b[2] * uu * uu; rr.append((d - fv) / sg)
        print(f"| {N} | {uu:.5f} | {100 * d:+.3f} | {100 * sg:.3f} | {100 * fv:+.3f} | {100 * (d - fv):+.3f} | {(d - fv) / sg:+.2f} |")
    print(f"\nsum of (residual / sigma)^2 = {sum(x * x for x in rr):.2f} = F2's chi2 {F2['chi2']:.2f} with {F2['dof']} dof, p = {F2['p']:.4f} (form not adequate: p < 0.01)")
    # (e)
    rows = []
    for e in V.GEN.ETA_B:
        F2 = fc[e].get("F2"); Z2 = fz[e].get("F2"); r = dict(eta=e, set="Phi" if e in V.PHI else "approach")
        r.update(F2_a_pct=100 * F2["beta"][0] if F2 else math.nan, F2_sigma_a_pct=100 * F2["se"][0] if F2 else math.nan,
                 z=F2["beta"][0] / F2["se"][0] if F2 else math.nan, F2_p=F2["p"] if F2 else math.nan,
                 form=("not adequate" if F2["p"] < 0.01 else "ok") if F2 else "-",
                 aZ_pct=100 * Z2["beta"][0] if Z2 else math.nan, aZ_sigma_pct=100 * Z2["se"][0] if Z2 else math.nan)
        for N in NS:
            E = res["cell"][V.GEN.cell_id("B", N, e, "default")]["E1"]; k = V.kr_cs(e)
            r[f"Delta_N{N}_pct"] = 100 * (E["cs"] / k - 1) if E else math.nan; r[f"Delta_N{N}_sigma_pct"] = 100 * E["sig"] / k if E else math.nan
        rows.append(r)
    os.makedirs(os.path.dirname(os.path.abspath(a.csv)), exist_ok=True)
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\n## (e) The table for the paper: {a.csv} ({len(rows)} rows: every density of block B; the registered numbers, arranged POST HOC)\n")
    print("| eta | set | a (F2) [%] | z | F2 p | a_Z [%] | Delta N 100 [%] | N 400 | N 900 | N 1600 |\n|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['eta']:.4g} | {r['set']} | {r['F2_a_pct']:+.3f} +- {r['F2_sigma_a_pct']:.3f} | {r['z']:+.2f} | {r['F2_p']:.3f} | {r['aZ_pct']:+.4f} +- {r['aZ_sigma_pct']:.4f} | "
              + " | ".join(f"{r[f'Delta_N{N}_pct']:+.3f} +- {r[f'Delta_N{N}_sigma_pct']:.3f}" for N in NS) + " |")


if __name__ == "__main__":
    main()
