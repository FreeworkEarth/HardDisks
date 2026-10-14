#!/usr/bin/env python3
"""##CHRIS 2026-10-13 (plan-author decision 46 item 5, decision 45 item 5 d; registration 261012 sec. 4.8.7.7): THE CONTINUITY CHECKS AS
INFORMATION LINES of the groups present, labelled INTERIM, on the Mac from the fetched folder (read only). The INTERIM tables do not print
them (s1_verdict_261012.report stops after the per-cell table in --interim); here the same analysis (s1_verdict_261012.analyse, interim,
its jackknife off: C-a, C-b and C-c use E0, nu_arg and Z, none of which depends on it) and the full report's own three lines, the same
f-strings, unchanged. Information, no verdict (sec. 4.8.7.7: "a miss is reported as a difference, it fails nothing").
usage (from hspist3/): python3 cluster/gen3_stage1_261012/s1_continuity_261012.py --data <the fetched folder>
"""
import argparse, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s1_verdict_261012 as V


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); a = ap.parse_args()
    res = V.analyse(a.data, interim=True, jackknife=False, quiet=True)
    print(f"# INTERIM -- the continuity checks (sec. 4.8.7.7) as information lines (s1_continuity_261012.py; groups present: "
          f"{', '.join(res['groups_present']) or 'none'}; expected |z| < 3; a miss is reported as a difference, it fails nothing)\n")
    c = res["cell"].get(V.GEN.cell_id("B", 100, math.pi / 8, "default"))
    if c and c["E0"]: print(f"- C-a: E0 at (100, pi/8) {c['E0']['cs']:.5f} +- {c['E0']['err_scaled']:.5f} against A1 v2 {V.REF_CA[0]} +- {V.REF_CA[1]}: "
                            f"z = {V.zdiff(c['E0']['cs'], c['E0']['err_scaled'], *V.REF_CA):+.2f}")
    for M, (m, s) in V.REF_CB.items():
        r = res["mass"].get((V.GEN.cell_id("B", 100, math.pi / 8, "default"), M))
        if r: print(f"- C-b: nu_arg at (100, pi/8), M = {M}: {r['nu_arg']:.7f} +- {r['nu_arg_se']:.2e} against Test G on KOA {m} +- {s}: "
                    f"z = {V.zdiff(r['nu_arg'], r['nu_arg_se'], m, s):+.2f}")
    c = res["cell"].get(V.GEN.cell_id("B", 400, math.pi / 8, "default"))
    if c: print(f"- C-c: Z at (400, pi/8) {c['Z'][0]:.5f} +- {c['Z'][1]:.5f} against {V.REF_CC[0]} +- {V.REF_CC[1]}: z = {V.zdiff(c['Z'][0], c['Z'][1], *V.REF_CC):+.2f}")


if __name__ == "__main__":
    main()
