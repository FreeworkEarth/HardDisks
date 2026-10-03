#!/usr/bin/env python3
"""##CHRIS 2026-10-10: does -ffp-contract=off move Paper 1's physics? One cell, eta = 0.112200
(the canonical cell nearest 0.10), 9 masses x 3 seeds on the no-fuse v1 binary, against the
canonical 25-seed run with contraction on. The estimator is IMPORTED from the canonical script
(cell(): TD = 200, X_EDGE = 2.5, same _prefix/_spectrum/argmax), so the binary is the only variable.
Pass rule, as pre-registered in the RED branch: c_s within 2 sigma (combined) of the 260919 value."""
import os, sys, math, csv, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
from paper1_populate_cs_err_20261002 import cell, slope_with_errors
ETA, L0 = 0.112200, 34.9999
# The canonical A1v2 data live in T.DROOT/<leaf>/m_<M>/ -- the populate script's own convention --
# NOT in T.A1_DIR (that is the older r25 campaign, 37.5 periods/trace, which TD = 200 cannot read;
# pointing at it produced c_s = 5.28 +- 2.45 with chi2_red 560 and a void "pass"). Retracted.
CANON  = os.path.join(T.DROOT, "eta_0p112200")
NOFUSE = os.path.join(T.ROOT, "nofuse_check_20261010", "eta_0p112200")
REF = {}
with open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv")) as fh:
    for r in csv.DictReader(fh):
        # ##CHRIS 2026-10-02 (Task H): select by the RECORDED eta. Since methods sec. 14.2 the table's `eta` column holds
        # eta_true (0.112267 for this cell) and `eta_rec` the recorded value; the old line was
        #     if abs(float(r["eta"]) - ETA) < 1e-6: REF = r
        if abs(float(r.get("eta_rec") or r["eta"]) - ETA) < 1e-6: REF = r
cs_ref, e_ref = float(REF["c_s"]), float(REF["c_s_err_scaled"])
# the table's c_s is at L_eff,true (sec. 14.2), so the recomputation uses the same acoustic length; a table without the
# column (pre-correction) gives the recorded l_eff. The binary comparison is a ratio, so its sigma is geometry-independent.
LE = float(REF.get("L_eff_true") or T.l_eff(L0))

def per_mass(d, label, max_runs=None):
    out = {}
    for M in T.A1_MASSES:
        runs = T.cell_runs(os.path.join(d, f"m_{M}"), M)
        if max_runs: runs = [q for q in runs if q[0] < max_runs]
        out[M] = cell((ETA, L0, M, runs))
    return out
A = per_mass(CANON, "canonical (25 seeds, contraction on)")
A3 = per_mass(CANON, "canonical, runs 0-2 only", max_runs=3)
B = per_mass(NOFUSE, "no-fuse (3 seeds)")

print(f"reference (260919 table): c_s = {cs_ref} +- {e_ref} (scaled)\n")
print("| M | nu canonical (25) | nu canonical runs 0-2 | nu no-fuse (3) | no-fuse vs canonical-25 |")
print("|---|---|---|---|---|")
for M in T.A1_MASSES:
    a, a3, b = A[M], A3[M], B[M]
    ea = a["sd"]/math.sqrt(a["n"]); eb = (b["sd"]/math.sqrt(b["n"])) if b["n"] > 1 else float("nan")
    z = abs(b["nu"]-a["nu"])/math.hypot(ea, eb) if np.isfinite(eb) else float("nan")
    print(f"| {M} | {a['nu']:.6f} ± {ea:.6f} (n={a['n']}) | {a3['nu']:.6f} (n={a3['n']}) | {b['nu']:.6f} ± {eb:.6f} (n={b['n']}) | {z:.1f}σ |")

def cs_of(D):
    ms = [M for M in T.A1_MASSES if D[M]["n"] > 1]
    x = np.array([T.x_of(M, L0) * T.l_eff(L0) / LE for M in ms]); y = np.array([D[M]["nu"] for M in ms])   # old: T.x_of(M, L0)
    sy = np.array([D[M]["sd"]/math.sqrt(D[M]["n"]) for M in ms])
    return slope_with_errors(x, y, sy)
sA, eA, eAs, chA = cs_of(A); sB, eB, eBs, chB = cs_of(B)
print(f"\nc_s canonical, recomputed here : {sA:.5f} ± {eAs:.5f} (chi2_red {chA:.2f})   [table: {cs_ref} ± {e_ref}]")
# SELF-CHECK: the estimator must reproduce the published value before the binary comparison means anything.
rep = abs(sA - cs_ref) <= 5e-6 * max(1.0, abs(cs_ref))
print(f"estimator reproduces the table to 5e-6? {'YES' if rep else '**NO -- comparison is VOID**'}")
print(f"c_s no-fuse, 3 seeds           : {sB:.5f} ± {eBs:.5f} (chi2_red {chB:.2f})")
zt = abs(sB - cs_ref)/math.hypot(eBs, e_ref)
print(f"\n### no-fuse vs 260919 table: |Δc_s| = {abs(sB-cs_ref):.5f}  ->  **{zt:.2f} sigma combined**  ->  {('PASS' if zt < 2 else 'FAIL') if rep else 'VOID (estimator self-check failed)'} (rule: < 2)")
