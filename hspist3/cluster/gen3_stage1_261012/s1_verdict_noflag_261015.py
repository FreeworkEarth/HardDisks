#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 54 item 2; registration 261012 sec. 4.8.6.2, 4.8.7.1, 4.8.7.2, 4.8.7.6): WHAT THE VERDICT
SCRIPT'S report() LEAVES OUT of the registered printout, printed from s1_verdict_261012's own analyse() and functions, unchanged: nothing
is fitted, estimated or selected in any other way. Its docstring promises "the weighted mean a-bar over Phi (F2 and F1), also without
flagged cells and densities"; report() prints it with all cells only.
  (1) Verdict K's z per density (z = a / sigma_a of F2): the 12 numbers of Phi that X^2 and max |z| are made of; the approach cells' z as
      information (no verdict).
  (2) The flags: stationarity (|hold t| >= 3 or |record z| >= 3; every cell of B, D and S), low-Q (a used mass with Q_k < 5), "form not
      adequate" (F2 p < 0.01), START-DEPENDENT, and the (cell, mass) folders the group rule leaves out.
  (3) Every pooled result also without the flagged cells and densities ("Flagged cells stay in; every pooled result is also printed
      without them", sec. 4.8.6.2 and 4.8.7.6; sec. 4.8.7.1 for the densities): the fits per eta (F2, F1, F1 without N 100) and the
      pressure fits with the flagged cells (stationarity or low-Q) left out, the rows that change; Verdict K's X^2, p and max |z| and the
      weighted mean a over Phi (F2 and F1) without the flagged densities and cells. INFORMATION: the registered Verdict K is the one the
      verdict script prints.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_verdict_noflag_261015.py --data <the fetched folder>
"""
import argparse, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s1_verdict_261012 as V


def stat_flag(c):
    t, z = c["hold_t"][0], c["rec_z"][0]
    return bool((math.isfinite(t) and abs(t) >= V.T_STAT) or (math.isfinite(z) and abs(z) >= V.T_STAT))


def fmt_fit(f):
    F2, F1, F1x = f.get("F2"), f.get("F1"), f.get("F1x")
    return (f"{f['n']} ({', '.join(str(n) for n in f['Ns'])}) | "
            + (f"{100 * F2['beta'][0]:+.3f} +- {100 * F2['se'][0]:.3f} | {F2['p']:.3f}" if F2 else "no F2 (fewer than four usable sizes) | -") + " | "
            + (f"{100 * F1['beta'][0]:+.3f} +- {100 * F1['se'][0]:.3f} | {F1['p']:.3f}" if F1 else "- | -") + " | "
            + (f"{100 * F1x['beta'][0]:+.3f} +- {100 * F1x['se'][0]:.3f}" if F1x else "-"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); a = ap.parse_args()
    res = V.analyse(a.data, interim=False, jackknife=True, quiet=True)
    print("# What the verdict script's report() leaves out, from its own analyse() and functions (s1_verdict_noflag_261015.py)\n")
    print(f"groups present: {', '.join(res['groups_present'])}; stopped by the stop rule: {res['stopped'] or 'none'}\n")
    K = res["K"]
    print("## (1) Verdict K's z per density (z = a / sigma_a of F2, sec. 4.8.7.2)\n\n| eta | set | F2 a [%] | sigma_a [%] | z | F2 p |\n|---|---|---|---|---|---|")
    for e in V.GEN.ETA_B:
        F2 = res["fitsC"][e].get("F2")
        if not F2: print(f"| {e:.4g} | {'Phi' if e in V.PHI else 'approach'} | no F2 | | | |"); continue
        print(f"| {e:.4g} | {'Phi' if e in V.PHI else 'approach (information)'} | {100 * F2['beta'][0]:+.3f} | {100 * F2['se'][0]:.3f} | "
              f"{F2['beta'][0] / F2['se'][0]:+.2f}{' (>= 3.5)' if e in V.PHI and abs(F2['beta'][0] / F2['se'][0]) >= V.Z_K else ''} | {F2['p']:.3f} |")
    if K:
        print(f"\nOn Phi: X^2 = {K['X2']:.2f} ({K['dof']} dof), p = {K['p']:.4f}, max |z| = {K['maxz']:.2f}; |z| < 2 at {sum(abs(v) < 2 for v in K['z'].values())} "
              f"of {len(K['z'])}; |z| >= 3.5 at: {', '.join(f'{e:.4g}' for e in K['bad']) or 'none'}")
    # (2) the flags
    print("\n## (2) The flags\n")
    cells = res["cell"]
    st = [cid for cid, c in cells.items() if stat_flag(c)]
    ht = max(((abs(c["hold_t"][0]), cid, c["hold_t"][0]) for cid, c in cells.items() if math.isfinite(c["hold_t"][0])))
    rz = max(((abs(c["rec_z"][0]), cid, c["rec_z"][0]) for cid, c in cells.items() if math.isfinite(c["rec_z"][0])))
    print(f"- stationarity (sec. 4.8.7.6, |t| >= 3 or |z| >= 3), over all {len(cells)} cells of B, D and S: {', '.join(st) or 'none'}; "
          f"the largest |hold t| {ht[2]:+.2f} ({ht[1]}), the largest |record z| {rz[2]:+.2f} ({rz[1]}); tests: "
          f"{sum(math.isfinite(c['hold_t'][0]) for c in cells.values())} hold, {sum(math.isfinite(c['rec_z'][0]) for c in cells.values())} record")
    lq = [cid for cid, c in cells.items() if c["E1"] and c["E1"]["lowQ"]]
    for cid in lq:
        ms = [(M, r["Q"]) for (c2, M), r in res["mass"].items() if c2 == cid and r["alpha"] >= 1 and not r["flagged"] and math.isfinite(r.get("Q", math.nan)) and r["Q"] < 5]
        print(f"- low-Q (sec. 4.8.6.2, a used mass with Q_k < 5): {cid}: " + ", ".join(f"M {M} Q_k {q:.2f}" for M, q in sorted(ms)))
    if not lq: print("- low-Q: none")
    fna = [e for e in V.GEN.ETA_B if res["fitsC"][e].get("F2") and res["fitsC"][e]["F2"]["p"] < 0.01]
    print(f"- form not adequate (sec. 4.8.7.1, F2 p < 0.01): " + (", ".join(f"{e:.4g} (p {res['fitsC'][e]['F2']['p']:.4f}, {'Phi' if e in V.PHI else 'approach'})" for e in fna) or "none"))
    dep = [k for k, v in res["S"].items() if v["independent"] is False]
    print(f"- START-DEPENDENT (sec. 4.8.7.5): {', '.join(f'N {N} eta {e:.2f}' for N, e in dep) or 'none'} (of {len(res['S'])} cells)")
    gr = sorted(((cid, M, r["lost"], r["fitfail"], r["excluded"]) for (cid, M), r in res["mass"].items() if r["flagged"]))
    print(f"- the group rule (sec. 4.8.5; left out of their cells' fits by the script itself): {len(gr)} (cell, mass) folders: "
          + "; ".join(f"{cid} M {M} ({lost} lost: {ff} fit failures, {ex} excluded)" for cid, M, lost, ff, ex in gr))
    # (3) without the flagged cells and densities
    print("\n## (3) Every pooled result also without the flagged cells and densities (information; the registered verdict is the script's)\n")
    flagged_cells = set(st) | set(lq)
    print(f"flagged cells left out: {', '.join(sorted(flagged_cells)) or 'none'}; flagged densities left out of the pooled statements: "
          f"{', '.join(f'{e:.4g}' for e in fna) or 'none'}\n")
    depset = set(dep); fitsC, fitsZ = {}, {}
    for e in V.GEN.ETA_B:                                            # the loop of analyse(), with the flagged cells left out
        pc, pz = [], []
        for N in V.NS:
            if (N, e) in depset: continue
            cell = cells.get(V.GEN.cell_id("B", N, e, "default"))
            if not cell or cell["cid"] in flagged_cells: continue
            if cell["E1"]:
                kr = V.kr_cs(e); pc.append((N, cell["E1"]["cs"] / kr - 1, cell["E1"]["sig"] / kr))
            zk = V.z_kr(e)
            if math.isfinite(cell["Z"][0]): pz.append((N, cell["Z"][0] / zk - 1, cell["Z"][1] / zk))
        fitsC[e] = V.fits_for(pc); fitsZ[e] = V.fits_for(pz)
    ch = [e for e in V.GEN.ETA_B if fmt_fit(fitsC[e]) != fmt_fit(res["fitsC"][e])]
    print("the finite-size fits per eta without the flagged cells, the rows that change (E1, Delta [%]):\n\n"
          "| eta | sizes | F2 a | F2 p | F1 a | F1 p | F1 without N 100: a |\n|---|---|---|---|---|---|---|")
    for e in ch: print(f"| {e:.4g} | {fmt_fit(fitsC[e])} |")
    print(f"\n{len(V.GEN.ETA_B) - len(ch)} of {len(V.GEN.ETA_B)} rows are unchanged (no flagged cell at those densities)")
    chz = [e for e in V.GEN.ETA_B if fmt_fit(fitsZ[e]) != fmt_fit(res["fitsZ"][e])]
    print("\nthe pressure fits without the flagged cells, the rows that change (Delta^Z [%]):\n\n| eta | sizes | a_Z | b_Z | F2 chi2 |\n|---|---|---|---|---|")
    for e in chz:
        F2 = fitsZ[e].get("F2")
        print(f"| {e:.4g} | {fitsZ[e]['n']} ({', '.join(str(n) for n in fitsZ[e]['Ns'])}) | "
              + (f"{100 * F2['beta'][0]:+.4f} +- {100 * F2['se'][0]:.4f} | {100 * F2['beta'][1]:+.3f} | {F2['chi2']:.2f} |" if F2 else "no F2 (fewer than four usable sizes) | | |"))
    print(f"\n{len(V.GEN.ETA_B) - len(chz)} of {len(V.GEN.ETA_B)} rows are unchanged")
    keep = [e for e in V.PHI if e not in fna]
    K2 = V.verdict_k({e: fitsC[e] for e in keep})
    if K2:
        print(f"\nVerdict K's statistics on Phi without the flagged densities and cells ({K2['dof']} densities: {', '.join(f'{e:.4g}' for e in keep)}): "
              f"X^2 = {K2['X2']:.2f} ({K2['dof']} dof), p = {K2['p']:.4f}, max |z| = {K2['maxz']:.2f}; |z| >= 3.5 at: {', '.join(f'{e:.4g}' for e in K2['bad']) or 'none'} "
              f"(information; the rule's thresholds would read {'KR-CONSISTENT' if K2['consistent'] else 'KR-INCONSISTENT'})")
    for lab, ks in (("Phi, all cells and densities (as the script prints it)", (V.PHI, res["fitsC"])), ("Phi without the flagged densities and cells", (keep, fitsC))):
        es, fc = ks
        a2 = V.wmean([(fc[e]["F2"]["beta"][0], fc[e]["F2"]["se"][0]) for e in es if "F2" in fc[e]])
        a1 = V.wmean([(fc[e]["F1"]["beta"][0], fc[e]["F1"]["se"][0]) for e in es if "F1" in fc[e]])
        print(f"weighted mean a over {lab}: F2 {100 * a2[0]:+.3f} +- {100 * a2[1]:.3f} % ({a2[2]}); F1 {100 * a1[0]:+.3f} +- {100 * a1[1]:.3f} % ({a1[2]})")


if __name__ == "__main__":
    main()
