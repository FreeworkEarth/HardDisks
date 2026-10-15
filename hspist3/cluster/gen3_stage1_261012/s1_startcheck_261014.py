#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (plan-author decision 49 item 3): THE START CHECK OF BLOCK S AT INTERIM, WITHOUT THE VERDICT.
s1_verdict_261012.py --interim prints the per-cell tables only (sec. 4.8.9); its analyse() computes Verdict S's start check (sec. 4.8.7.5)
all the same. This prints, per block S cell (N, eta), the four registered quantities per start -- (i) Z and (ii) psi6 over the 20
trajectories of both masses, (iii) nu_d-bar at alpha = 1 and (iv) at alpha = 10 over 10 each -- with the registered
z = (mean_0 - mean_turned) / sqrt(SE_0^2 + SE_turned^2) as analyse() computes it (res["S"]), and, as information, the main block's default
start (Z, psi6) where its group is present. No verdict: the START-INDEPENDENT / START-DEPENDENT column is not printed, and the threshold
is not applied. From the fetched folder only (read only).
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_startcheck_261014.py --data <a folder holding S, and the B groups present>
"""
import argparse, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s1_verdict_261012 as V


def fm(v, fmt):
    return format(v, fmt) if isinstance(v, (int, float)) and math.isfinite(v) else "-"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); a = ap.parse_args()
    res = V.analyse(a.data, interim=True, quiet=True)
    print(f"# INTERIM -- block S: the start check without the verdict (s1_startcheck_261014.py; groups present: "
          f"{', '.join(res['groups_present']) or 'none'})\n")
    print("per start: the mean +- SE (n); z = (start0 - turned) / sqrt(SE_0^2 + SE_turned^2), from s1_verdict_261012.analyse; "
          "nu_d-bar per internal time unit (24 units = 1 sigma-time); the main block's default start as information\n")
    print("| N | eta | Z start0 | Z turned | z Z | psi6 start0 | psi6 turned | z psi6 | nu_d-bar alpha 1: start0 | turned | z | "
          "nu_d-bar alpha 10: start0 | turned | z | main block, default start: Z, psi6 |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (N, e), v in sorted(res["S"].items()):
        c0 = res["cell"].get(V.GEN.cell_id("S", N, e, "0")); ct = res["cell"].get(V.GEN.cell_id("S", N, e, "turned"))
        z = v["z"]; cols = []
        for k in ("Z", "psi6"):
            fmt = ".5f" if k == "Z" else ".4f"
            for c in (c0, ct):
                m, se, n = c[k][0], c[k][1], (c[k][2] if len(c[k]) > 2 else "")
                cols.append(f"{fm(m, fmt)} +- {fm(se, fmt)} ({n})")
            cols.append(fm(z[k], "+.2f"))
        for al, lab in ((1, "nud_a1"), (10, "nud_a10")):
            M = int(round(al * N))
            for c in (c0, ct):
                r = res["mass"].get((c["cid"], M))
                cols.append(f"{fm(r['nud'], '.6g')} +- {fm(r['nud_se'], '.2g')} ({r['present'] - r['excluded'] - r['fitfail']})" if r else "-")
            cols.append(fm(z[lab], "+.2f"))
        b = res["cell"].get(V.GEN.cell_id("B", N, e, "default"))
        cols.append(f"{fm(b['Z'][0], '.5f')}, {fm(b['psi6'][0], '.4f')}" if b else "- (its group not present)")
        print(f"| {N} | {e:.2f} | " + " | ".join(cols) + " |")
    zs = [abs(x) for v in res["S"].values() for x in v["z"].values() if math.isfinite(x)]
    print(f"\n{len(res['S'])} cells; {len(zs)} finite z of {4 * len(res['S'])}; the largest |z|: {max(zs):.2f}" if zs else "\nno finite z")


if __name__ == "__main__":
    main()
