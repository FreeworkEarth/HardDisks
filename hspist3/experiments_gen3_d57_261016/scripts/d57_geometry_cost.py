#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 57, addendum 3, task 3; notes 261012 sec. 4.7.140): THE REGISTRATION'S GEOMETRY TABLE AND THE
DESIGN'S COST. No run; nothing registered (the registration is the plan author's).
  GEOMETRY, from s2_vacancy_table's long CSV (cluster/gen3_stage2_261016/s2_vacancy_table_261015_heights.csv: every height H = k/24 sigma of
  the range, the binary's lattice fit replicated, the strtof check):
    - the 12 cells eta 0.70, 0.72, 0.74, 0.76, 0.78, 0.80 x N 900, 1600, pin turned: among the arrangements the fit chooses with at most 1
      vacancy per compartment, the one whose heights lie nearest the stage 1 height (k 720 for N 900, 960 for N 1600; the distance = the
      smallest |k - k_ref| over its heights); within it, for each target start strain -0.6 %, 0 and +0.6 %, the height whose start strain e
      is nearest the target (a tie: the k nearer k_ref);
    - the isotropy check, (900, 0.76) and (900, 0.78), pin 0: the same with 0 vacancies.
    Each row: H (the --height text), k, L0 = N pi / (8 H eta), the aspect W / H (W = L0 - 0.025), e, the gap, the vacancies, |e - target|.
  COST of the design: 3 heights x 5 masses (alpha 1, 2, 5, 10, 20) x 10 seeds per cell; 12 cells turned and 2 pin 0 (2100 runs).
    - events per run: P2's package model (gen_p2_koa_261016.plan): r = 1.2004 (Z - 1) + 0.0112 per disk and sigma-time times N (T_eq + T_rec)
      times the rehearsal's factor per pin (1.002 turned, 1.171 pin 0), T_eq = 5000, T_rec = 200 periods of the binary's nu_pred
      (gen_stage1_koa_261012.nu_pred_binary) rounded up to whole steps. Z: D56-1's planning values (d56_a_timing.z_eos: Liu's global branch
      at 0.70, Liu's solid branch from 0.72) as the LOWER value, and those times 1.268 from 0.72 on (P2's package value) as the UPPER;
      P2's turned cells executed 0.84-1.01 of the upper (Report D57-2, table B).
    - microseconds per event per CPU model: P2's measured medians at N 900 (Report D57-2, table A); at N 1600 times D56-1's measured
      N 1600 / N 900 ratio (E5-2670 v3: 3.836 / 3.533; E5-4640 v2: 4.580 / 4.242; the v3 ratio for the four models D56-1 did not see).
    - the KOA mix: P2's runs per model (32 EPYC 7543P, 59 EPYC 9654, 125 E5-2670 v3, 10 E5-4640 v2, 133 Gold 6248, 41 Gold 6342 of 400).
    - core-hours = the sum over the runs of events x the mix's mean microseconds per event; the wall time at 64 cores = core-hours / 64;
      the longest run on the slowest model (E5-4640 v2) and on E5-2670 v3; the merge scaled from P2's (47 min for 400 runs, of which
      p2_reduce 2542 s).
Writes geometry_cost/: geometry.csv, cost_runs.csv, geometry_cost.md.
usage (from hspist3/): python3 experiments_gen3_d57_261016/scripts/d57_geometry_cost.py [--out <a new folder>]
"""
import math, os, sys
import numpy as np, pandas as pd
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d57_261016")
for p in (os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import plot_speed_of_sound_edmd as SOS
CSV = os.path.join(HS, "cluster", "gen3_stage2_261016", "s2_vacancy_table_261015_heights.csv")
OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else os.path.join(E, "geometry_cost")
ETAS = (0.70, 0.72, 0.74, 0.76, 0.78, 0.80); KREF = {900: 720, 1600: 960}; TARGETS = (-0.006, 0.0, 0.006)
ALPHAS = (1, 2, 5, 10, 20); NSEED = 10; TEQ = 5000.0; RATE_A, RATE_B, ZX = 1.2004, 0.0112, 1.268; MEAS = {"0": 1.171, "turned": 1.002}
US900 = {"EPYC 9654": 2.139, "Gold 6342": 2.209, "EPYC 7543P": 2.284, "Gold 6248": 2.764, "E5-2670 v3": 3.555, "E5-4640 v2": 3.935}   # Report D57-2 A
MIX = {"EPYC 7543P": 32, "EPYC 9654": 59, "E5-2670 v3": 125, "E5-4640 v2": 10, "Gold 6248": 133, "Gold 6342": 41}                  # of 400
N1600 = {"E5-2670 v3": 3.836 / 3.533, "E5-4640 v2": 4.580 / 4.242}                                                                  # D56-1 (a)
P2_MERGE_MIN, P2_RUNS = (23 * 60 + 9 + 59 / 60) - (22 * 60 + 22 + 45 / 60), 400


def z_eos(eta):
    a = np.array([eta])
    return float(SOS.Z_liu_global(a)[0]) if eta < 0.72 - 1e-12 else float(SOS.Z_liu_solid(a)[0])


def pick(df, N, eta, pin, maxvac):
    q = df[(df["N"] == N) & np.isclose(df["eta"], eta) & (df["pin"] == pin) & (df["vacancies"] >= 0) & (df["vacancies"] <= maxvac) & (df["height_px_ok"] == 1)]
    best = None
    for arr, g in q.groupby("arrangement"):
        dist = int((g["k_px"] - KREF[N]).abs().min())
        if best is None or dist < best[0]: best = (dist, arr, g)
    dist, arr, g = best; out = []
    for t in TARGETS:
        i = (g["e"] - t).abs().add((g["k_px"] - KREF[N]).abs() * 1e-12).idxmin(); r = g.loc[i]
        out.append(dict(N=N, eta=eta, pin=pin, arrangement=arr, vacancies=int(r["vacancies"]), arr_k_range=f"{g['k_px'].min()}-{g['k_px'].max()}",
                        arr_dist_to_kref=dist, target_e=t, k=int(r["k_px"]), H_text=str(r["H_text"]), H=int(r["k_px"]) / 24.0,
                        L0=N * math.pi / (8.0 * (int(r["k_px"]) / 24.0) * eta), aspect=float(r["aspect_W_over_H"]), e=float(r["e"]), gap=float(r["gap_sigma"]),
                        miss=abs(float(r["e"]) - t)))
    return out


def main():
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    df = pd.read_csv(CSV, dtype={"pin": str, "H_text": str}); rows = []
    for N in (900, 1600):
        for eta in ETAS: rows += pick(df, N, eta, "turned", 1)
    for eta in (0.76, 0.78): rows += pick(df, 900, eta, "0", 0)
    G = pd.DataFrame(rows)
    dup = G.groupby(["N", "eta", "pin"])["k"].nunique(); assert (dup == 3).all(), f"two targets share a height: {dup[dup != 3]}"
    mix_n = sum(MIX.values()); us = {N: sum(MIX[m] / mix_n * US900[m] * (N1600.get(m, N1600['E5-2670 v3']) if N == 1600 else 1.0) for m in MIX) for N in (900, 1600)}
    slow = {N: US900["E5-4640 v2"] * (N1600["E5-4640 v2"] if N == 1600 else 1.0) for N in (900, 1600)}
    v3 = {N: US900["E5-2670 v3"] * (N1600["E5-2670 v3"] if N == 1600 else 1.0) for N in (900, 1600)}
    C = []
    for g in G.to_dict("records"):
        zl = z_eos(g["eta"]); zu = zl * (ZX if g["eta"] >= 0.72 - 1e-12 else 1.0)
        for a in ALPHAS:
            M = int(round(a * g["N"])); nup = GEN.nu_pred_binary(g["N"], g["H"], g["L0"], M); Trec = math.ceil((200.0 / nup) / GEN.DTF) * GEN.DTF
            for lab, z in (("lower", zl), ("upper", zu)):
                ev = (RATE_A * (z - 1.0) + RATE_B) * g["N"] * (TEQ + Trec) * MEAS[g["pin"]]
                C.append(dict(N=g["N"], eta=g["eta"], pin=g["pin"], k=g["k"], target_e=g["target_e"], alpha=a, M=M, nu_pred=nup, T_rec=Trec, Z=z, model=lab, events=ev,
                              min_mix=ev * us[g["N"]] / 6e7, min_v2=ev * slow[g["N"]] / 6e7, min_v3=ev * v3[g["N"]] / 6e7, runs=NSEED))
    C = pd.DataFrame(C); os.makedirs(OUT)
    G.to_csv(os.path.join(OUT, "geometry.csv"), index=False, float_format="%.10g"); C.to_csv(os.path.join(OUT, "cost_runs.csv"), index=False, float_format="%.8g")
    L = ["# The registration's geometry table and the design's cost (decision 57 addendum 3, task 3; d57_geometry_cost.py)", "",
         "## The geometry", "",
         "Per cell: the arrangement with at most 1 vacancy (pin 0: 0 vacancies) whose heights lie nearest the stage 1 height (k 720 / 960), and in it the "
         "heights whose START strain is nearest -0.6 %, 0, +0.6 % (s2_vacancy_table_261015_heights.csv; every height passes the strtof check).", "",
         "| N | eta | pin | arrangement | vacancies | its heights (k) and their distance to the stage 1 height | target e | H [sigma] (k) | L0 | aspect W/H | e (start) | |e - target| | gap [sigma] |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for g in G.to_dict("records"):
        L.append(f"| {g['N']} | {g['eta']:.2f} | {g['pin']} | {g['arrangement']} | {g['vacancies']} | {g['arr_k_range']} ({g['arr_dist_to_kref']}) | {100 * g['target_e']:+.1f} % | "
                 f"{g['H_text']} ({g['k']}) | {g['L0']:.6f} | {g['aspect']:.4f} | {100 * g['e']:+.3f} % | {100 * g['miss']:.3f} % | {g['gap']:.4f} |")
    L += ["", "## The cost", "",
          f"3 heights x 5 masses x {NSEED} seeds per cell; 12 cells turned and 2 pin 0: {int(C[C['model'] == 'lower']['runs'].sum())} runs. Events per run: P2's "
          f"package model with D56-1's planning Z (lower) and that Z x {ZX} from 0.72 on (upper; P2's turned cells executed 0.84-1.01 of the upper). "
          f"Microseconds per event: P2's measured medians per CPU model (N 1600: x D56-1's N 1600 / N 900 ratio); the KOA mix of P2: "
          f"mean {us[900]:.3f} us at N 900, {us[1600]:.3f} at N 1600.", "",
          "| N | eta | pin | Z lower / upper | runs | events per run [M], alpha 1 .. 20 (upper) | minutes per run, the mix (upper): median (largest) | "
          "the longest run (upper): E5-2670 v3 / E5-4640 v2 [min] | core-hours lower / upper |", "|---|---|---|---|---|---|---|---|---|"]
    tot = {"lower": 0.0, "upper": 0.0}
    for (N, eta, pin), g in C.groupby(["N", "eta", "pin"], sort=True):
        u, l = g[g["model"] == "upper"], g[g["model"] == "lower"]
        chl, chu = float((l["min_mix"] * l["runs"]).sum() / 60), float((u["min_mix"] * u["runs"]).sum() / 60); tot["lower"] += chl; tot["upper"] += chu
        mm = np.repeat(u["min_mix"].to_numpy(), u["runs"].to_numpy())
        L.append(f"| {N} | {eta:.2f} | {pin} | {l['Z'].iloc[0]:.2f} / {u['Z'].iloc[0]:.2f} | {int(l['runs'].sum())} | "
                 f"{u['events'].min() / 1e6:.0f} .. {u['events'].max() / 1e6:.0f} | {np.median(mm):.2f} ({mm.max():.2f}) | {u['min_v3'].max():.2f} / {u['min_v2'].max():.2f} | "
                 f"{chl:.1f} / {chu:.1f} |")
    nr = int(C[C["model"] == "lower"]["runs"].sum()); up = C[C["model"] == "upper"]
    L += ["", "The minutes per run per CPU model (upper events; each run's events x the model's microseconds per event; over the design's runs):", "",
          "| CPU model | P2 runs on it (of 400) | us per event N 900 / N 1600 | N 900: minutes per run, median (largest) | N 1600: minutes per run, median (largest) |",
          "|---|---|---|---|---|"]
    for m in sorted(US900, key=US900.get):
        f16 = N1600.get(m, N1600["E5-2670 v3"]); cells = []
        for N, k in ((900, 1.0), (1600, f16)):
            g = up[up["N"] == N]; mm = np.repeat((g["events"] * US900[m] * k / 6e7).to_numpy(), g["runs"].to_numpy()); cells.append(f"{np.median(mm):.2f} ({mm.max():.2f})")
        L.append(f"| {m} | {MIX[m]} | {US900[m]:.3f} / {US900[m] * f16:.3f} | {cells[0]} | {cells[1]} |")
    L += ["", f"- total: {nr} runs; {tot['lower']:.0f} core-hours (lower) to {tot['upper']:.0f} (upper) at P2's KOA mix; at 64 cores {tot['lower'] / 64:.1f} to "
              f"{tot['upper'] / 64:.1f} h of wall time (the array's packing, no queue);",
          f"- the longest single run (upper): {up['min_v3'].max():.1f} min on E5-2670 v3, {up['min_v2'].max():.1f} min on E5-4640 v2 "
          f"(P2's rule, 1.5 x the longest on the slower model rounded up to 5 min: {int(math.ceil(1.5 * up['min_v2'].max() / 5) * 5)} min);",
          f"- for comparison, P2: 400 runs, 36.9 core-hours used against its package's 52.9 (upper);",
          f"- the merge: P2's took {P2_MERGE_MIN:.0f} min for {P2_RUNS} runs (p2_reduce 2542 s of it); scaled by the run count to {nr} runs: about "
          f"{P2_MERGE_MIN * nr / P2_RUNS / 60:.1f} h (N 1600's trajectories are longer; the merge's 24 h limit)."]
    open(os.path.join(OUT, "geometry_cost.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
