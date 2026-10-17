#!/usr/bin/env python3
"""Independent verification of D57 addendum 3 task 3, part 1: the geometry selection. Read-only on the repo."""
import math, os, sys
import numpy as np, pandas as pd

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
OUT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HS, "cluster/gen3_stage2_261016/s2_vacancy_table_261015_heights.csv")
CC = os.path.join(HS, "experiments_gen3_d57_261016/geometry_cost/geometry.csv")

d = pd.read_csv(SRC, dtype={"pin": str, "H_text": str})
KREF = {900: 720, 1600: 960}
TG = [-0.006, 0.0, 0.006]

# sanity: pin column vs arrangement suffix
suf = d["arrangement"].str.split().str[-1]
print("pin column vs arrangement suffix mismatches:", int((suf != d["pin"]).sum()), "of", len(d))
# sanity: H_text == k/24, L0 from formula
Hk = d["k_px"] / 24.0
print("max |float(H_text) - k/24|:", float((d["H_text"].astype(float) - Hk).abs().max()))
L0f = d["N"] * math.pi / (8 * Hk * d["eta"])
print("max |L0_table - N pi/(8 H eta)|:", float((d["L0"] - L0f).abs().max()))
print("max |aspect_table - (L0-0.025)/H|:", float((d["aspect_W_over_H"] - (L0f - 0.025) / Hk).abs().max()))

sel = []
notes = []
def choose(N, eta, pin, vac_ok):
    q = d[(d.N == N) & (np.abs(d.eta - eta) < 1e-9) & (d.pin == pin) & vac_ok(d.vacancies) & (d.height_px_ok == 1)]
    kr = KREF[N]
    dist = {a: int((g.k_px - kr).abs().min()) for a, g in q.groupby("arrangement")}
    order = sorted(dist.items(), key=lambda t: t[1])
    best = order[0][1]
    winners = [a for a, v in order if v == best]
    runner = order[1] if len(order) > 1 else (None, None)
    notes.append(dict(N=N, eta=eta, pin=pin, n_arr=len(order), best=winners, best_dist=best,
                      second=runner[0], second_dist=runner[1], all=order))
    arr = winners[0]
    g = q[q.arrangement == arr].copy()
    for t in TG:
        g["m"] = (g.e - t).abs()
        mn = g["m"].min()
        cand = g[np.isclose(g["m"], mn, rtol=0, atol=1e-12)].copy()
        cand["kd"] = (cand.k_px - kr).abs()
        cand = cand.sort_values("kd")
        r = cand.iloc[0]
        k = int(r.k_px); H = k / 24.0; L0 = N * math.pi / (8 * H * eta)
        sel.append(dict(N=N, eta=eta, pin=pin, arrangement=arr, n_tied_arr=len(winners), arr_kmin=int(g.k_px.min()), arr_kmax=int(g.k_px.max()),
                        arr_dist=best, target_e=t, k=k, H=H, H_text=r.H_text, L0=L0, aspect=float(r.aspect_W_over_H),
                        aspect_recalc=(L0 - 0.025) / H, e=float(r.e), gap=float(r.gap_sigma), vacancies=int(r.vacancies),
                        miss=abs(float(r.e) - t), n_e_ties=len(cand),
                        second_best_miss=float(np.sort(g["m"].unique())[1]) if g["m"].nunique() > 1 else np.nan))

for N in (900, 1600):
    for eta in (0.70, 0.72, 0.74, 0.76, 0.78, 0.80):
        choose(N, eta, "turned", lambda v: (v >= 0) & (v <= 1))
for eta in (0.76, 0.78):
    choose(900, eta, "0", lambda v: v == 0)

S = pd.DataFrame(sel)
S.to_csv(os.path.join(OUT, "my_geometry.csv"), index=False, float_format="%.12g")
print("\nselected rows:", len(S))
print("cells with a tie in the arrangement distance:", [(n["N"], n["eta"], n["pin"], n["best"]) for n in notes if len(n["best"]) > 1])
print("\nper cell: winner, dist, runner-up, dist")
for n in notes:
    print(f"  N {n['N']} eta {n['eta']:.2f} pin {n['pin']}: {n['best']} d={n['best_dist']} | 2nd {n['second']} d={n['second_dist']} | n_arr={n['n_arr']}")
print("rows with an e tie:", int((S.n_e_ties > 1).sum()))
print("duplicate k within a cell:", int((S.groupby(['N', 'eta', 'pin']).k.nunique() != 3).sum()))
print("\npoorly met targets (|e - target| > 0.002):")
print(S[S.miss > 0.002][["N", "eta", "pin", "target_e", "k", "e", "miss"]].to_string() if (S.miss > 0.002).any() else "  none")
print("max miss:", S.miss.max(), " at", S.loc[S.miss.idxmax(), ["N", "eta", "pin", "target_e", "k"]].to_dict())

# compare with CC
C = pd.read_csv(CC, dtype={"pin": str, "H_text": str})
key = ["N", "eta", "pin", "target_e"]
S["eta_r"] = S.eta.round(6); C["eta_r"] = C.eta.round(6)
m = S.merge(C, on=["N", "eta_r", "pin", "target_e"], suffixes=("_me", "_cc"), how="outer", indicator=True)
print("\nmerge:", m["_merge"].value_counts().to_dict())
def rel(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))
print("arrangement equal:", bool((m.arrangement_me == m.arrangement_cc).all()))
print("k equal:", bool((m.k_me == m.k_cc).all()))
print("vacancies equal:", bool((m.vacancies_me == m.vacancies_cc).all()))
print("H_text equal:", bool((m.H_text_me == m.H_text_cc).all()))
print("arr range equal:", bool(((m.arr_kmin.astype(str) + "-" + m.arr_kmax.astype(str)) == m.arr_k_range).all()))
print("arr dist equal:", bool((m.arr_dist == m.arr_dist_to_kref).all()))
print("max |e diff|:", float((m.e_me - m.e_cc).abs().max()), " max rel:", rel(m.e_me, m.e_cc))
print("max |gap diff|:", float((m.gap_me - m.gap_cc).abs().max()), " max rel:", rel(m.gap_me, m.gap_cc))
print("max rel L0 diff:", rel(m.L0_me, m.L0_cc), " max |L0 diff|:", float((m.L0_me - m.L0_cc).abs().max()))
print("max rel H diff:", rel(m.H_me, m.H_cc))
print("max |aspect diff| (table value):", float((m.aspect_me - m.aspect_cc).abs().max()))
print("max |aspect_recalc - CC aspect|:", float((m.aspect_recalc - m.aspect_cc).abs().max()))
print("max |miss diff|:", float((m.miss_me - m.miss_cc).abs().max()))
S.drop(columns="eta_r").to_csv(os.path.join(OUT, "my_geometry.csv"), index=False, float_format="%.12g")
pd.DataFrame([dict(N=n["N"], eta=n["eta"], pin=n["pin"], winner="|".join(n["best"]), dist=n["best_dist"],
                   second=n["second"], second_dist=n["second_dist"], n_arr=n["n_arr"],
                   all=";".join(f"{a}:{v}" for a, v in n["all"])) for n in notes]).to_csv(os.path.join(OUT, "my_arrangement_ranking.csv"), index=False)
