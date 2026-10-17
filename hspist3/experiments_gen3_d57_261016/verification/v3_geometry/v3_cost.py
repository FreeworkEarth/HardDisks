#!/usr/bin/env python3
"""Independent verification of D57 addendum 3 task 3, part 2: the cost. Uses my own geometry (my_geometry.csv). Read-only on the repo."""
import math, os, sys
import numpy as np, pandas as pd

HS = "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3"
OUT = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(HS, "cluster/gen3_stage1_261012"), os.path.join(HS, "validation"), HS):
    sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import plot_speed_of_sound_edmd as SOS

G = pd.read_csv(os.path.join(OUT, "my_geometry.csv"), dtype={"pin": str, "H_text": str})
CCC = pd.read_csv(os.path.join(HS, "experiments_gen3_d57_261016/geometry_cost/cost_runs.csv"), dtype={"pin": str})

US900 = {"EPYC 9654": 2.139, "Gold 6342": 2.209, "EPYC 7543P": 2.284, "Gold 6248": 2.764, "E5-2670 v3": 3.555, "E5-4640 v2": 3.935}
W = {"EPYC 7543P": 32, "EPYC 9654": 59, "E5-2670 v3": 125, "E5-4640 v2": 10, "Gold 6248": 133, "Gold 6342": 41}
assert sum(W.values()) == 400
def us(model, N):
    if N == 900: return US900[model]
    r = 4.580 / 4.242 if model == "E5-4640 v2" else 3.836 / 3.533
    return US900[model] * r

def Zlow(eta):
    if eta < 0.72 - 1e-12:
        return float(SOS.Z_liu_global(np.array([eta]))[0])
    return float(SOS.Z_liu_solid(np.array([eta]))[0])
def Zup(eta):
    return Zlow(eta) * (1.268 if eta >= 0.72 - 1e-12 else 1.0)

print("Z values: eta, lower, upper; and Z_liu_global vs Z_liu_solid at each eta")
for eta in (0.70, 0.72, 0.74, 0.76, 0.78, 0.80):
    print(f"  {eta:.2f}  {Zlow(eta):.6f}  {Zup(eta):.6f}   global {float(SOS.Z_liu_global(np.array([eta]))[0]):.6f} solid {float(SOS.Z_liu_solid(np.array([eta]))[0]):.6f}")

rows = []
for g in G.to_dict("records"):
    N, eta, pin, k = int(g["N"]), float(g["eta"]), g["pin"], int(g["k"])
    H = k / 24.0
    L0 = N * math.pi / (8.0 * H * eta)
    fpin = 1.002 if pin == "turned" else 1.171
    for a in (1, 2, 5, 10, 20):
        M = a * N
        nup = GEN.nu_pred_binary(N, H, L0, M)
        Trec = math.ceil((200.0 / nup) / GEN.DTF) * GEN.DTF
        for lab, Z in (("lower", Zlow(eta)), ("upper", Zup(eta))):
            ev = (1.2004 * (Z - 1.0) + 0.0112) * N * (5000.0 + Trec) * fpin
            r = dict(N=N, eta=eta, pin=pin, k=k, target_e=g["target_e"], alpha=a, M=M, nu_pred=nup, T_rec=Trec, Z=Z, model=lab, events=ev)
            # core-seconds over the mix: sum over models of weight x events x us
            r["coresec_mix_per_run"] = sum(W[m] / 400.0 * ev * us(m, N) * 1e-6 for m in W)
            for m in W:
                r["min_" + m] = ev * us(m, N) * 1e-6 / 60.0
            rows.append(r)
C = pd.DataFrame(rows)
C.to_csv(os.path.join(OUT, "my_cost_runs.csv"), index=False, float_format="%.12g")
NRUN = 10
print("\nruns:", len(C[C.model == "lower"]) * NRUN)

# mix means
for N in (900, 1600):
    print(f"mix mean us/event N {N}: {sum(W[m] / 400 * us(m, N) for m in W):.6f}")
print("us/event at N 1600 per model:", {m: round(us(m, 1600), 6) for m in US900})

# totals
tot = {lab: float((C[C.model == lab].coresec_mix_per_run * NRUN).sum() / 3600) for lab in ("lower", "upper")}
print(f"\nTOTAL core-hours lower {tot['lower']:.4f} upper {tot['upper']:.4f}; wall at 64 cores {tot['lower']/64:.4f} .. {tot['upper']/64:.4f} h")
U = C[C.model == "upper"]
print(f"longest single run (upper) E5-2670 v3: {U['min_E5-2670 v3'].max():.4f} min; E5-4640 v2: {U['min_E5-4640 v2'].max():.4f} min")
i = U['min_E5-4640 v2'].idxmax(); print("  at", U.loc[i, ["N", "eta", "pin", "k", "alpha", "events", "T_rec"]].to_dict())
print(f"P2 time-limit rule: ceil(1.5 x {U['min_E5-4640 v2'].max():.3f} / 5) x 5 = {int(math.ceil(1.5 * U['min_E5-4640 v2'].max() / 5) * 5)} min")

# per cell
print("\nper cell: N eta pin | Zl Zu | runs | events upper min..max [M] | mix minutes upper median (max) | v3 / v2 longest | core-h lower / upper")
percell = []
for (N, eta, pin), g in C.groupby(["N", "eta", "pin"]):
    l, u = g[g.model == "lower"], g[g.model == "upper"]
    mixmin = (u.coresec_mix_per_run / 60).to_numpy()
    chl = float(l.coresec_mix_per_run.sum() * NRUN / 3600); chu = float(u.coresec_mix_per_run.sum() * NRUN / 3600)
    rec = dict(N=N, eta=eta, pin=pin, Zl=l.Z.iloc[0], Zu=u.Z.iloc[0], runs=len(l) * NRUN, ev_min=u.events.min() / 1e6, ev_max=u.events.max() / 1e6,
               mix_med=float(np.median(np.repeat(mixmin, NRUN))), mix_max=float(mixmin.max()), v3max=float(u['min_E5-2670 v3'].max()),
               v2max=float(u['min_E5-4640 v2'].max()), ch_lower=chl, ch_upper=chu)
    percell.append(rec)
    print(f"  {N} {eta:.2f} {pin:6s} | {rec['Zl']:.4f} {rec['Zu']:.4f} | {rec['runs']} | {rec['ev_min']:.3f} .. {rec['ev_max']:.3f} | "
          f"{rec['mix_med']:.4f} ({rec['mix_max']:.4f}) | {rec['v3max']:.4f} / {rec['v2max']:.4f} | {chl:.4f} / {chu:.4f}")
pd.DataFrame(percell).to_csv(os.path.join(OUT, "my_cost_percell.csv"), index=False, float_format="%.10g")

print("\nper model (upper): N 900 median (max), N 1600 median (max) minutes per run")
for m in sorted(US900, key=US900.get):
    s = []
    for N in (900, 1600):
        v = U[U.N == N]["min_" + m].to_numpy()
        s.append(f"{np.median(np.repeat(v, NRUN)):.4f} ({v.max():.4f})")
    print(f"  {m:11s} us {us(m,900):.4f} / {us(m,1600):.4f}  N900 {s[0]}  N1600 {s[1]}")

# compare row by row with CC's cost_runs.csv
CCC["eta_r"] = CCC.eta.round(6); C["eta_r"] = C.eta.round(6)
m = C.merge(CCC, on=["N", "eta_r", "pin", "k", "target_e", "alpha", "model"], suffixes=("_me", "_cc"), how="outer", indicator=True)
print("\ncost_runs merge:", m["_merge"].value_counts().to_dict())
def rel(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); return float(np.max(np.abs(a - b) / np.abs(b)))
for col in ("M", "nu_pred", "T_rec", "Z", "events"):
    print(f"  max rel diff {col}: {rel(m[col + '_me'], m[col + '_cc']):.3e}")
print(f"  max rel diff min_mix: {rel(m.coresec_mix_per_run / 60, m.min_mix):.3e}")
print(f"  max rel diff min_v2: {rel(m['min_E5-4640 v2'], m.min_v2):.3e}")
print(f"  max rel diff min_v3: {rel(m['min_E5-2670 v3'], m.min_v3):.3e}")
print("  runs all 10:", bool((m.runs == 10).all()))
