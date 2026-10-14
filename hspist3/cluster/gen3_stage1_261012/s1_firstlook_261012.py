#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 46 item 5): THE FIRST LOOK AT THE PHYSICS of the groups present, labelled INTERIM, on the Mac from
the fetched folder (read only), from s1_verdict_261012.analyse (interim, with its jackknife: the INTERIM tables' numbers) and nothing else:
  per eta, for N = 100 and N = 400 side by side: c_s of E1 and of E0 against KR, Delta = c_s / c_s^KR - 1 [%] with its error (E1: sigma; E0:
  err_scaled) and z = Delta / error;
  per eta: Z per cell against Z_KR, Delta^Z = Z / Z_KR - 1 [%] with its SE;
  the low-Q flags: per cell, the used masses (alpha >= 1, not flagged) with Q_k < 5, and the smallest Q_k of every cell.
No fit over N, no verdict, no prediction (sec. 4.8.9: the fits and Verdict K wait for all four sizes). Phi = the 12 densities <= 0.60.
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_firstlook_261012.py --data <the fetched folder> [--sizes 100,400]
"""
import argparse, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s1_verdict_261012 as V


def cell(res, N, e): return res["cell"].get(V.GEN.cell_id("B", N, e, "default"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--sizes", default="100,400"); a = ap.parse_args()
    Ns = [int(x) for x in a.sizes.split(",")]
    res = V.analyse(a.data, interim=True, jackknife=True, quiet=True)
    print(f"# INTERIM -- a first look at the physics (s1_firstlook_261012.py; groups present: {', '.join(res['groups_present']) or 'none'}; "
          "no fit over N, no verdict)\n")
    print("## c_s against KR per cell: Delta = c_s / c_s^KR - 1 [%] +- its error (E1: sigma; E0: err_scaled), z = Delta / error\n")
    print("| eta | set | c_s^KR | " + " | ".join(f"N {N}: E1 Delta [%] (z; masses) | N {N}: E0 Delta [%] (z)" for N in Ns) + " |\n|---|---|---|"
          + "---|" * (2 * len(Ns)))
    for e in V.GEN.ETA_B:
        kr = V.kr_cs(e); out = []
        for N in Ns:
            c = cell(res, N, e)
            if not c: out += ["-", "-"]; continue
            e1, e0 = c["E1"], c["E0"]
            out.append(f"{100 * (e1['cs'] / kr - 1):+.3f} +- {100 * e1['sig'] / kr:.3f} ({(e1['cs'] / kr - 1) / (e1['sig'] / kr):+.1f}; {e1['n']})" if e1 else "none (< 5 masses)")
            out.append(f"{100 * (e0['cs'] / kr - 1):+.3f} +- {100 * e0['err_scaled'] / kr:.3f} ({(e0['cs'] / kr - 1) / (e0['err_scaled'] / kr):+.1f})" if e0 else "-")
        print(f"| {e:.4g} | {'Phi' if e in V.PHI else 'approach'} | {kr:.5f} | " + " | ".join(out) + " |")
    print("\n## Z per cell against Z_KR: Delta^Z = Z / Z_KR - 1 [%] +- SE (n trajectories)\n")
    print("| eta | set | Z_KR | " + " | ".join(f"N {N}: Z +- SE | N {N}: Delta^Z [%] (z)" for N in Ns) + " |\n|---|---|---|" + "---|" * (2 * len(Ns)))
    for e in V.GEN.ETA_B:
        zk = V.z_kr(e); out = []
        for N in Ns:
            c = cell(res, N, e)
            if not c or not math.isfinite(c["Z"][0]): out += ["-", "-"]; continue
            z, se, n = c["Z"]
            out += [f"{z:.5f} +- {se:.5f} ({n})", f"{100 * (z / zk - 1):+.3f} +- {100 * se / zk:.3f} ({(z / zk - 1) / (se / zk):+.1f})"]
        print(f"| {e:.4g} | {'Phi' if e in V.PHI else 'approach'} | {zk:.5f} | " + " | ".join(out) + " |")
    print("\n## The low-Q flags (sec. 4.8.6.2: a cell in which any used mass has Q_k < 5)\n")
    print("| N | eta | low-Q flag | used masses with Q_k < 5 (M: Q_k) | smallest Q_k of the used masses (M) |\n|---|---|---|---|---|")
    for N in Ns:
        for e in V.GEN.ETA_B:
            c = cell(res, N, e)
            if not c: continue
            rows = [r for (cid, M), r in res["mass"].items() if cid == c["cid"] and r["alpha"] >= 1 and not r["flagged"] and math.isfinite(r.get("Q", math.nan))]
            low = sorted((r["M"], r["Q"]) for r in rows if r["Q"] < 5)
            mn = min(rows, key=lambda r: r["Q"]) if rows else None
            print(f"| {N} | {e:.4g} | {'FLAG' if c['E1'] and c['E1']['lowQ'] else '-'} | {', '.join(f'{M}: {q:.2f}' for M, q in low) or '-'} | "
                  + (f"{mn['Q']:.2f} ({mn['M']})" if mn else "-") + " |")


if __name__ == "__main__":
    main()
