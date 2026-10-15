#!/usr/bin/env python3
"""##CHRIS 2026-10-14 (decision 49 addendum 1 item 4): for option (i), the legacy held path with L0 rounded to the 1/48 grid, the KR
stiffness N_s c_KR^2 / L_eff^2 at the rounded L0 against the registered L0, per N = 400 cell (the size of the rescaling the identity's
two sides would need, by the KR model). Read-only; the functions of validation/paper1_confinement_prereg_20261012.py."""
import importlib.util, math, os, sys
HS = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
spec = importlib.util.spec_from_file_location("prereg", os.path.join(HS, "validation", "paper1_confinement_prereg_20261012.py"))
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage1_261012")); import gen_stage1_koa_261012 as GEN
N, NS, H = 400, 200, 20.0
print("| eta | L0 | L0 on the 1/48 grid | eta there | k_S^KR there / at L0 - 1 |\n|---|---|---|---|---|")
for e in GEN.ETA_B:
    L0 = N * math.pi / (8 * H * e); g = round(L0 * 48) / 48; eg = P.eta_of(NS, H, g)
    k = lambda L, et: NS * P.cs2(et) / P.leff(L) ** 2
    print(f"| {e:.4g} | {L0:.6f} | {g:.6f} | {eg:.6f} | {100 * (k(g, eg) / k(L0, e) - 1):+.3f} % |")
