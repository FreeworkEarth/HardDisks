#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 35, item 4c; 261012 sec. 4.7.68): THE TEST R TABLES for the update of the meeting figures 1
and 2 (task M1). Existing data only, no runs. Every number comes from the committed Test R evidence (main 7ed5e102), read with git show:
the per-arm tables merged/testR/<gen2|gen3>/epi8_H_H10_L10/m_<M>/nud_runs.csv (KOA's reduce_testR_koa.py; nu_arg = red_nu.csv's registered
argmax, nu_d the registered fit, fit_ok) and the printed verdict testR_verdict_output.txt.
  fig1_testR_koa_means.csv       per engine, mass and estimator: n, mean, SE = SD/sqrt(n), mean -+ 1.96 SE
  fig1_testR_koa_difference.csv  per mass and estimator: d = mean gen3 / mean gen2 - 1, SE(d), the 95 % interval d -+ 1.959964 SE(d), z,
                                 by the verdict's own dstat() (imported unchanged), and the registered per-mass outcome
  fig2_deltaM_testR.csv          the verdict's printed delta_M table (INFORMATION; against c_300): estimator, engine, M, delta [%], SE [%]
CHECKS (printed; the script stops on a failure): every cell of the verdict's registered table (means, d, SE(d), interval, z) reproduced
at the printed precision; the per-mass outcomes and the verdict line as printed; the delta table: 3 estimators x 2 engines x 2 masses.
usage (from hspist3/): python3 experiments_m1_figures_261011/m1_tables_testR_261011.py [--out DIR]
"""
import argparse, io, math, os, re, subprocess, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE); ROOT = os.path.dirname(HS)
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_testR_261011")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import testR_koa_verdict_261011 as V                              # the registered verdict's dstat, mass_verdict and constants
COMMIT = "7ed5e102"                                               # Report R1's evidence commit
LOCR = "hspist3/experiments_gen3_testR_koa_261011"
FAILS = []


def git(*a): return subprocess.run(["git", "-C", ROOT, *a], check=True, capture_output=True, text=True).stdout


def check(label, ok, detail=""):
    print(f"- {label}: {'PASS' if ok else 'FAIL'}{(' -- ' + detail) if detail else ''}")
    if not ok: FAILS.append(label)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.join(HERE, "tables")); a = ap.parse_args()
    full = git("rev-parse", COMMIT).strip()
    print(f"# M1 update (decision 35, item 4c): the Test R tables, printed by experiments_m1_figures_261011/m1_tables_testR_261011.py; "
          f"existing data only\n\nthe source: main {full[:8]} ({LOCR}), read with git show\n")
    nud = {(e, M): pd.read_csv(io.StringIO(git("show", f"{full}:{LOCR}/merged/testR/{e}/{V.CID}/m_{M}/nud_runs.csv"))).sort_values("run")
           for e in V.ENG for M in V.MASSES}
    means, diffs = [], []
    for M in V.MASSES:
        out = {}
        for lab, col in (("nu_arg", "nu_arg"), ("nu_d", "nu_d")):
            arr = {e: (nud[(e, M)][col][nud[(e, M)]["fit_ok"] == 1] if col == "nu_d" else nud[(e, M)][col]).to_numpy(float) for e in V.ENG}
            for e in V.ENG:
                v = arr[e]; se = v.std(ddof=1) / math.sqrt(len(v))
                means.append(dict(engine=e, M=M, estimator=lab, n=len(v), mean=v.mean(), se=se, lo95=v.mean() - 1.96 * se, hi95=v.mean() + 1.96 * se))
            s = V.dstat(arr["gen3"], arr["gen2"]); out[lab] = s
            diffs.append(dict(M=M, estimator=lab, n_gen3=s["n3"], n_gen2=s["n2"], mean_gen3=s["m3"], mean_gen2=s["m2"], d_pct=100 * s["d"],
                              se_pct=100 * s["se"], lo95_pct=100 * s["lo"], hi95_pct=100 * s["hi"], z=s["z"]))
        o = V.mass_verdict(out["nu_arg"]["z"], out["nu_d"]["z"], out["nu_d"]["lo"], out["nu_d"]["hi"])
        for r in diffs[-2:]: r["outcome"] = o
    Mn, D = pd.DataFrame(means), pd.DataFrame(diffs)
    os.makedirs(a.out, exist_ok=True)
    Mn.to_csv(os.path.join(a.out, "fig1_testR_koa_means.csv"), index=False); D.to_csv(os.path.join(a.out, "fig1_testR_koa_difference.csv"), index=False)
    vt = git("show", f"{full}:{LOCR}/testR_verdict_output.txt")
    reg = vt.split("## The registered numbers", 1)[1].split("\n## ", 1)[0]          # the registered table's section only
    rows = [l for l in reg.splitlines() if re.match(r"^\| (50|100|300) \| nu_(arg|d) \|", l)]
    check("the verdict's registered table: 6 rows", len(rows) == 6)
    for l in rows:
        c = [x.strip() for x in l.strip().strip("|").split("|")]; d = D[(D.M == int(c[0])) & (D.estimator == c[1])].iloc[0]
        mine = [str(d.n_gen3), str(d.n_gen2), f"{d.mean_gen3:.7f}", f"{d.mean_gen2:.7f}", f"{d.d_pct:+.4f}", f"{d.se_pct:.4f}",
                f"[{d.lo95_pct:+.4f}, {d.hi95_pct:+.4f}]", f"{d.z:+.2f}"]
        check(f"M = {c[0]} {c[1]} reproduced", c[2:] == mine, "" if c[2:] == mine else f"printed {c[2:]} / here {mine}")
    line = "; ".join(f"M = {M}: {D[D.M == M].outcome.iloc[0]}" for M in V.MASSES)
    check("the per-mass outcomes as printed", line in vt, line)
    check("the verdict line as printed", "TEST R: ENGINES AGREE" in vt)
    drows = re.findall(r"^\| (nu_arg|nu_d|nu_0) \| (gen2|gen3) \| ([-+][\d.]+) \+- ([\d.]+) \| ([-+][\d.]+) \+- ([\d.]+) \|$", vt, re.M)
    dl = []
    for est, e, v50, s50, v100, s100 in drows:
        for M, v, s in ((50, v50, s50), (100, v100, s100)):
            dl.append(dict(estimator=est, engine=e, M=M, delta_pct=float(v), se_pct=float(s), lo95_pct=float(v) - 1.96 * float(s),
                           hi95_pct=float(v) + 1.96 * float(s), reference="c_300 (the one mass >= 300 of Test R)", source=f"Report R1 ({full[:8]})"))
    check("the verdict's delta table: 3 estimators x 2 engines x 2 masses", len(dl) == 12, f"{len(dl)} rows")
    pd.DataFrame(dl).to_csv(os.path.join(a.out, "fig2_deltaM_testR.csv"), index=False)
    print("\n| M | estimator | d [%] | 95 % interval [%] | z | outcome |\n|---|---|---|---|---|---|")
    for _, d in D.iterrows(): print(f"| {d.M} | {d.estimator} | {d.d_pct:+.4f} | [{d.lo95_pct:+.4f}, {d.hi95_pct:+.4f}] | {d.z:+.2f} | {d.outcome} |")
    print(f"\nCHECKS: {'all PASS' if not FAILS else 'FAIL: ' + '; '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
