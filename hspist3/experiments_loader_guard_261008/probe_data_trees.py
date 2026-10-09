# ##CHRIS 2026-10-08 (261012 sec. 4.7.4, decision 2): probe of the loader guard on the real data trees (read-only): every
# directory holding an accelerated run record must be refused; every directory under the data roots of the paper scripts
# must pass. usage (from hspist3/): python3 experiments_loader_guard_261008/probe_data_trees.py
import os, sys, re
HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
sys.path.insert(0, os.path.join(HS, "validation"))
import edmd_acc_guard as G
import provenance_edmd_acc_261009 as P
rows = P.scan_tree()
acc_dirs = sorted({os.path.dirname(os.path.join(P.REPO, rel)) for rel, cls, h in rows
                   if cls in ("run record", "run log", "summary csv") and any(k == "ACCELERATED" for k, _ in h)})
def refused(d):
    try: G.guard(d); return False
    except G.AcceleratedRunError: return True
# 1. every directory that holds an accelerated record (outside the experiments root itself) must be refused
inroot = [d for d in acc_dirs if not os.path.basename(d).startswith("experiments")]
miss = [d for d in inroot if not refused(d)]
print(f"directories holding an accelerated record: {len(acc_dirs)} ({len(acc_dirs) - len(inroot)} of them an experiments root, never scanned by design)")
print(f"  of the others refused by the guard: {len(inroot) - len(miss)} / {len(inroot)}")
for d in miss[:10]: print("  NOT REFUSED:", os.path.relpath(d, P.REPO))
for d in acc_dirs:
    if os.path.basename(d).startswith("experiments"): print("  experiments root with an accelerated record:", os.path.relpath(d, P.REPO))
# 2. every directory under the data roots the paper scripts read must pass
ROOTS = ["experiments_speed_of_sound/EDMD/mode1_normalized_units/00_eta_sweep_ROMAN/" + x for x in
         ("A1v2_20260914", "famB_20260911", "A2_dilute50_20260917", "A2_dilute_20260916", "A2_long200_20260915", "A2_topup_20260912",
          "A2_alpha2_20260912", "campaign_r25_psi6_20260823", "routeA_lowdensity_20260912", "confinement_B_20261013",
          "confinement_pilot_20261013", "tests_20260913")]
ROOTS += ["experiments_energy_transfer/" + x for x in sorted(os.listdir(os.path.join(HS, "experiments_energy_transfer")))
          if x.startswith(("paper1_confinement_", "level"))]
ROOTS += ["experiments_resched_gate_261005"] if os.path.isdir(os.path.join(HS, "experiments_resched_gate_261005")) else []
nd = nref = 0; bad = []
for r in ROOTS:
    base = os.path.join(HS, r)
    if not os.path.isdir(base): print("  (missing root)", r); continue
    for dp, dn, fn in os.walk(base):
        nd += 1
        if refused(dp): nref += 1; bad.append(os.path.relpath(dp, HS))
print(f"\npaper data roots: {len(ROOTS)}; directories under them: {nd}; refused: {nref}")
for b in bad[:10]: print("  REFUSED:", b)
