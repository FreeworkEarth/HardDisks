#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.24; plan-author decision 12, part 2): P2's frequencies and c_s next to A1 v2's. PILOT,
EXPLORATORY, no paper use, no new runs: the existing P2 data (sec. 4.7.19) and A1 v2's canonical trajectories.
Per eta_true and mass (M = 50, 300, 2000): P2's mean nu on the full record with its SE (12 seeds), A1 v2's for the same eta and mass
(its 25 seeds less the loader's discards), and their ratio. The estimator is the registered, canonical one
(paper1_populate_cs_err_20261002.cell: per trajectory the argmax of the periodogram over the first TD = 200 predicted periods, bins
from round(TD / X_EDGE) = 80 on; per cell the mean and SD / sqrt(n)). A1 v2 goes through cell() itself; P2's traces are compressed
(.csv.gz, sec. 4.7.19 decision 10), so this script runs the same steps with a loader that also reads gzip, after the same
provenance guard. Per eta: c_s by the canonical through-origin slope of the per-mass nu against x_M = k(M / (2 N_s)) / (2 pi
L_eff), with the canonical error (propagated; scaled by max(1, sqrt(chi2_red))), for P2 and for A1 v2 on the same three masses,
and for A1 v2 on all nine (the canonical cell set). The box: P2 ran in A1 v2's true box, 2 L0_true with L0_true = L0 - delta/2
(the binary's 1/48-sigma grid), so L_eff = L0_true - 2r - t/2 for both (the box-truncation correction of methods sec. 14, the
canonical table's L_eff_true).
GATES, printed: (1) this script's loader and per-trajectory estimator, run on A1 v2's traces, give cell()'s per-cell nu exactly;
(2) A1 v2's nine-mass c_s and its two errors equal the canonical table (260919_A1v2_final_cs_vs_eta.csv) to 5e-6 relative.
Then the question of decision 12: c_s(0.7167) - c_s(0.7060) for P2 and for A1 v2, with its SE (the two cells independent).
Also psi6 at release and at the end, and nu, for every seed at eta 0.7167.
usage (from hspist3/ of the engine-gen3 worktree):
  python3 experiments_gen3_p2_261009/p2_nu_cs_table.py --out experiments_gen3_p2_261009/data > p2_nu_cs_table_output.txt
"""
import argparse, csv, glob, gzip, math, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p2_run as P2
sys.path.insert(0, os.path.join(P2.MAIN_HS, "validation")); sys.path.insert(0, P2.MAIN_HS)
import tests_20260913 as T
import edmd_acc_guard
import paper1_populate_cs_err_20261002 as C        # the canonical estimator: cell(), slope_with_errors(), box_delta(), TD, X_EDGE

CANON = T.plot_path(C.SRC19)                      # 260919_A1v2_final_cs_vs_eta.csv, the corrected canonical table
LEAF = {0.6905: "eta_0p690", 0.7060: "eta_0p705", 0.7113: "eta_0p710", 0.7167: "eta_0p715"}   # A1 v2's leaf per eta_true


def nu_of(p):
    """the canonical per-trajectory frequency (cell()'s three lines), from a .csv or .csv.gz trace, after the provenance guard"""
    import pandas as pd
    edmd_acc_guard.guard(p)
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", newline="") as fh:
        nup = float(next(csv.DictReader(fh))["Predicted_Frequency"])
    d = pd.read_csv(p, usecols=["Time", "Displacement(σ)"])
    t, x = d["Time"].to_numpy(float), d["Displacement(σ)"].to_numpy(float)
    dt = (t[-1] - t[0]) / (len(t) - 1)
    n = T._prefix(t, nup, C.TD)
    if n is None: return float("nan")
    P, df = T._spectrum(x[:n], dt)
    k = int(round(C.TD / C.X_EDGE))
    return (k + int(np.argmax(P[k:]))) * df


def mse(v):
    v = np.asarray(v, float)
    return float(v.mean()), (float(v.std(ddof=1)) / math.sqrt(len(v)) if len(v) > 1 else float("nan")), len(v)


def slope(Le, rows):
    """the canonical slope_with_errors on (x_M, mean nu, SE) with x_M at L_eff = Le; rows: [(M, mean, se)]"""
    x = np.array([T.k_root(M / (2.0 * T.N_SIDE)) / (2 * math.pi * Le) for M, _, _ in rows])
    return C.slope_with_errors(x, [m for _, m, _ in rows], [s for _, _, s in rows])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(); O = os.path.abspath(a.out)
    print("# P2's frequencies and c_s next to A1 v2's (261012 sec. 4.7.24; decision 12, part 2), printed by "
          "experiments_gen3_p2_261009/p2_nu_cs_table.py -- PILOT, EXPLORATORY, no paper use\n")
    table = {l["leaf"]: l for l in T.a1_leaf_table()}
    # ---- A1 v2: the canonical cell() on every (eta, mass) of the four leaves (all nine masses)
    tasks = []
    for eta_t, L0a, _ in P2.A1V2:
        leaf = table[LEAF[eta_t]]; L0 = float(leaf["L0"])
        assert abs(L0 - L0a) < 1e-12, (L0, L0a)
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            tasks.append((eta_t, L0, M, runs))
    # ---- P2: every seed's trace
    p2tr = {}
    for eta_t, L0a, _ in P2.A1V2:
        for M in P2.MASSES:
            for r in range(P2.NSEED):
                g = glob.glob(os.path.join(O, f"eta{eta_t:.4f}_M{M}", f"seed{r}", "wall_x_positions_*_run0.csv*"))
                assert len(g) == 1, (eta_t, M, r, g)
                p2tr[(eta_t, M, r)] = g[0]
    a1paths = [(eta_t, M, r, p) for eta_t, L0, M, runs in tasks if M in P2.MASSES for r, p, disc in runs if not disc]
    with Pool(a.workers) as pool:
        cells = pool.map(C.cell, tasks, chunksize=1)
        nus_p2 = dict(zip(p2tr.keys(), pool.map(nu_of, list(p2tr.values()), chunksize=4)))
        nus_a1 = pool.map(nu_of, [p for _, _, _, p in a1paths], chunksize=8)
    cellk = {(q["eta"], q["M"]): q for q in cells}
    # ---- gate 1: this script's loader and estimator on A1 v2's traces = cell()
    g1 = []
    for eta_t, _, _ in P2.A1V2:
        for M in P2.MASSES:
            mine = [v for (e, m, r, p), v in zip(a1paths, nus_a1) if e == eta_t and m == M]
            q = cellk[(eta_t, M)]
            g1.append(len(mine) == q["n"] and float(np.mean(mine)) == q["nu"])
    print(f"gate 1: this script's loader and estimator on A1 v2's traces give cell()'s mean nu exactly in "
          f"{sum(g1)} of {len(g1)} (eta, mass) cells{'' if all(g1) else '  **GATE FAILED**'}")
    # ---- gate 2: A1 v2 on all nine masses = the canonical table
    canon = {r["eta"]: r for r in csv.DictReader(open(CANON))}
    g2, c9 = [], {}
    for eta_t, L0a, _ in P2.A1V2:
        L0 = L0a; delta = C.box_delta(L0); LeT = T.l_eff(L0) - delta / 2
        rows = [(M, cellk[(eta_t, M)]["nu"], cellk[(eta_t, M)]["sd"] / math.sqrt(cellk[(eta_t, M)]["n"])) for M in T.A1_MASSES]
        s, e, es, x2 = slope(LeT, rows); c9[eta_t] = (s, e, es, x2, LeT, delta)
        etaT = table[LEAF[eta_t]]["eta"] * L0 / (L0 - delta / 2)
        r = canon.get(f"{etaT:.6f}")
        ok = r is not None and all(abs(float(r[k]) - v) <= 5e-6 * max(1.0, abs(v)) for k, v in (("c_s", s), ("c_s_err", e), ("c_s_err_scaled", es)))
        g2.append(ok)
    print(f"gate 2: A1 v2's nine-mass c_s, c_s_err and c_s_err_scaled equal the canonical table ({os.path.basename(CANON)}) to 5e-6 in "
          f"{sum(g2)} of {len(g2)} eta{'' if all(g2) else '  **GATE FAILED**'}\n")
    # ---- table 1: per eta_true and mass
    print("## 1. Per eta_true and mass: the mean frequency on the full record (200 predicted periods), canonical estimator\n")
    print("| eta_true | M | P2: n | P2: mean nu (SE) | A1 v2: n (discarded) | A1 v2: mean nu (SE) | ratio P2 / A1 v2 (SE) |\n|---|---|---|---|---|---|---|")
    p2rows, a1rows = {}, {}
    for eta_t, _, _ in P2.A1V2:
        for M in P2.MASSES:
            m2, s2, n2 = mse([nus_p2[(eta_t, M, r)] for r in range(P2.NSEED)])
            q = cellk[(eta_t, M)]; m1, s1 = q["nu"], q["sd"] / math.sqrt(q["n"])
            p2rows.setdefault(eta_t, []).append((M, m2, s2)); a1rows.setdefault(eta_t, []).append((M, m1, s1))
            rr = m2 / m1; sr = rr * math.sqrt((s2 / m2) ** 2 + (s1 / m1) ** 2)
            print(f"| {eta_t:.4f} | {M} | {n2} | {m2:.6f} ({s2:.6f}) | {q['n']} ({q['nd']}) | {m1:.6f} ({s1:.6f}) | {rr:.4f} ({sr:.4f}) |")
    # ---- table 2: c_s per eta
    print("\n## 2. Per eta_true: c_s by the canonical through-origin slope (error propagated; scaled = error x max(1, sqrt(chi2_red)))\n")
    print("| eta_true | box: P2 L0 (exact) / A1 v2 L0 as launched, delta, L0_true | L_eff P2 / A1 v2 [sigma] | c_s P2, M = 50, 300, 2000: "
          "c_s (err; scaled; chi2_red) | c_s A1 v2, same three masses | ratio P2 / A1 v2 (scaled SE) | c_s A1 v2, all nine masses (= the canonical table) |\n"
          "|---|---|---|---|---|---|---|")
    cs = {}
    for eta_t, L0a, _ in P2.A1V2:
        L0p = P2.grid_L0(L0a); LeP = T.l_eff(L0p)
        s9, e9, es9, x29, LeT, delta = c9[eta_t]
        sp = slope(LeP, p2rows[eta_t]); sa = slope(LeT, a1rows[eta_t]); cs[eta_t] = (sp, sa, (s9, e9, es9, x29))
        rr = sp[0] / sa[0]; sr = rr * math.sqrt((sp[2] / sp[0]) ** 2 + (sa[2] / sa[0]) ** 2)
        print(f"| {eta_t:.4f} | {L0p:.6f} / {L0a:.4f}, {delta:.6f}, {L0a - delta / 2:.6f} | {LeP:.6f} / {LeT:.6f} | "
              f"{sp[0]:.4f} ({sp[1]:.4f}; {sp[2]:.4f}; {sp[3]:.2f}) | {sa[0]:.4f} ({sa[1]:.4f}; {sa[2]:.4f}; {sa[3]:.2f}) | {rr:.4f} ({sr:.4f}) | "
              f"{s9:.4f} ({e9:.4f}; {es9:.4f}; {x29:.2f}) |")
    # ---- table 3: the question
    print("\n## 3. The question of decision 12: the side of the dip, c_s(0.7167) - c_s(0.7060) (SE from the two scaled errors)\n")
    print("| data | c_s(0.7060) | c_s(0.7167) | difference (SE) | relative to c_s(0.7060) [%] (SE) | difference / SE |\n|---|---|---|---|---|---|")
    for name, idx in (("P2, M = 50, 300, 2000", 0), ("A1 v2, the same three masses", 1), ("A1 v2, all nine masses (canonical)", 2)):
        a6, a7 = cs[0.7060][idx], cs[0.7167][idx]
        d = a7[0] - a6[0]; sd = math.sqrt(a7[2] ** 2 + a6[2] ** 2)
        rel = 100 * d / a6[0]; srel = 100 * math.sqrt((a7[2] / a6[0]) ** 2 + (a7[0] * a6[2] / a6[0] ** 2) ** 2)
        print(f"| {name} | {a6[0]:.4f} ({a6[2]:.4f}) | {a7[0]:.4f} ({a7[2]:.4f}) | {d:+.4f} ({sd:.4f}) | {rel:+.2f} ({srel:.2f}) | {d / sd:+.2f} |")
    # ---- table 4: psi6 at release, every seed at 0.7167
    print("\n## 4. eta_true 0.7167: psi6 (global) at the release and at the end, and nu, for every seed (P2's speed_of_sound_psi6.csv)\n")
    print("| M | seed index | run seed | psi6 at release | psi6 at the end | nu (full record) |\n|---|---|---|---|---|---|")
    for M in P2.MASSES:
        for r in range(P2.NSEED):
            row = list(csv.DictReader(open(os.path.join(O, f"eta0.7167_M{M}", f"seed{r}", "speed_of_sound_psi6.csv"))))[-1]
            print(f"| {M} | {r} | {row['seed']} | {float(row['psi6_global_hold']):.4f} | {float(row['psi6_global_end']):.4f} | {nus_p2[(0.7167, M, r)]:.6f} |")
    print("\n(PILOT, EXPLORATORY: no verdict rule. P2: gen3, M4 lattice start, divider held 1e4 sigma-time, 12 seeds per cell; A1 v2: gen2, "
          "the driver's grid seeding, a 2000-step hold, 25 seeds per cell. Both: N = 100, H = 10, the same true box per eta.)")


if __name__ == "__main__":
    main()
