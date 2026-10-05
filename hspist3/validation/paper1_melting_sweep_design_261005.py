#!/usr/bin/env python3
"""##CHRIS 2026-10-05 (Task EE): DRAFT design numbers for the melting size sweep (261012 sec. 4, DRAFT; not a pre-registration).

Cells: eta in {0.695, 0.700, 0.704, 0.708, 0.712, 0.716, 0.720} at N = 100, 400, 900 (N_s = 50, 200, 450), H = 10 sqrt(N/100),
L_0 on the 1/48 grid nearest N_s pi r^2/(eta H) (so 2 L_0 x 24 px is an integer: no box truncation); eta_true printed.
The N = 100 dip: from the canonical table (260919_A1v2_final_cs_vs_eta.csv, rows 0.695 <= eta_rec <= 0.720): depth
D = (c_max - c_min)/c_max between the window's maximum and minimum cells. Under M1 (Mayer-Wood loop) D(N) = D(100) (N/100)^-1/2.
Seeds: the dip at N = 900 is to be resolved at >= 5 sigma: sigma(D) <= D(900)/5, sigma(D) ~ sqrt(2) x the per-cell relative
error; the per-cell STATISTICAL error (c_s_err, seed-propagated, unscaled) is assumed independent of N and to fall as
seeds^-1/2 [INFERENCE: the confinement campaign's c_s_err did not trend with N_s at fixed seeds]. The chi2-scaled error
(c_s_err_scaled) of the N = 100 window cells is dominated by mass disagreement (chi2_red 1.8-14), which seeds do not reduce
[DERIVATION: err_scaled = err sqrt(chi2_red) -> the residual scatter when chi2_red >> 1]; it is printed as the floor.
Cost: per trajectory at N = 100 from the A1v2 run logs of the window cells (##RUN ... (s) / trajectories; Mac seconds) x the
measured KOA factor 1.804 (round_plan_261002.koa_speed); scaled to N by (N/100)^p per sigma-time and (N/100)^1/2 for the
200-period duration (period ~ L ~ sqrt(N) at fixed alpha). p = 2.87 is the measured local exponent at N_s 100 -> 200 of the
A-fixed pi/8 cells (fixed L_0, divider rate ~ N) [DATA, sec. 3.9 sacct]; p = 2.37 subtracts the 0.5 that the sweep geometry
saves (divider rate ~ H ~ sqrt(N)) [INFERENCE]. (a) the current engine, (b) a hypothetical 10x faster engine.
usage (from hspist3/): python3 validation/paper1_melting_sweep_design_261005.py
"""
import csv, glob, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS); sys.path.insert(0, os.path.join(HS, "cluster"))
import tests_20260913 as T
import round_plan_261002 as RP

ETAS = (0.695, 0.700, 0.704, 0.708, 0.712, 0.716, 0.720)
NS = (50, 200, 450)
MASS_ALPHA = (0.5, 1, 2, 3, 5, 7.5, 10, 15, 20)


def main():
    print("## Melting size sweep -- DRAFT design numbers (261012 sec. 4)\n")
    print("### Cells (H = 10 sqrt(N/100); L_0 on the 1/48 grid nearest the target eta)\n")
    print("| N | N_s | H | eta target | L_0 | 2 L_0 x 24 (px) | eta_true | masses M = alpha 2 N_s |\n|---|---|---|---|---|---|---|---|")
    for ns in NS:
        N = 2 * ns; H = 10 * math.sqrt(N / 100)
        for e in ETAS:
            L0 = round(ns * math.pi * 0.25 / (e * H) * 48) / 48; et = ns * math.pi * 0.25 / (H * L0)
            print(f"| {N} | {ns} | {H:g} | {e:.3f} | {L0:.6f} | {2 * L0 * 24:.1f} | {et:.6f} | {int(MASS_ALPHA[0] * 2 * ns)} ... {int(MASS_ALPHA[-1] * 2 * ns)} |")
    rows = [r for r in csv.DictReader(open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv"))) if 0.695 <= float(r["eta_rec"]) <= 0.7201]
    print("\n### The N = 100 window (canonical table, 25 seeds x 9 masses)\n")
    print("| eta_true | c_s | c_s_err (statistical) | rel. [%] | c_s_err_scaled | rel. [%] | chi2_red |\n|---|---|---|---|---|---|---|")
    for r in rows:
        c = float(r["c_s"])
        print(f"| {float(r['eta']):.4f} | {c:.4f} | {float(r['c_s_err']):.4f} | {100 * float(r['c_s_err']) / c:.2f} | {float(r['c_s_err_scaled']):.4f} | "
              f"{100 * float(r['c_s_err_scaled']) / c:.2f} | {float(r['chi2_red']):.2f} |")
    mx = max(rows, key=lambda r: float(r["c_s"])); mn = min([r for r in rows if float(r["eta"]) > float(mx["eta"])], key=lambda r: float(r["c_s"]))
    cmax, cmin = float(mx["c_s"]), float(mn["c_s"]); D = (cmax - cmin) / cmax
    rel_stat = math.sqrt(sum((float(r["c_s_err"]) / float(r["c_s"])) ** 2 for r in (mx, mn)) / 2)
    rel_scal = math.sqrt(sum((float(r["c_s_err_scaled"]) / float(r["c_s"])) ** 2 for r in (mx, mn)) / 2)
    print(f"\nN = 100: maximum at eta_true {float(mx['eta']):.4f} (c_s {cmax:.3f}), minimum at {float(mn['eta']):.4f} (c_s {cmin:.3f}); "
          f"depth D = (c_max - c_min)/c_max = {100 * D:.1f} %")
    print("\n### Expected depth and the seeds for a >= 5 sigma dip\n")
    print("| N | D under M1 [%] | sigma(D) needed [%] | per-cell rel. error needed [%] | seeds per mass (statistical error) | "
          "floor: N = 100 chi2-scaled rel. error [%] |\n|---|---|---|---|---|---|")
    seeds = {}
    for ns in NS:
        N = 2 * ns; DN = D * (N / 100) ** -0.5; sD = DN / 5; rel = sD / math.sqrt(2)
        nseed = max(25, math.ceil(25 * (rel_stat / rel) ** 2)); seeds[N] = nseed
        print(f"| {N} | {100 * DN:.2f} | {100 * sD:.2f} | {100 * rel:.2f} | {nseed} | {100 * rel_scal:.2f} |")
    # cost per trajectory at N = 100 in the window: A1v2 run logs (Mac) x KOA factor
    secs = n = 0
    for leaf in T.a1_leaf_table():
        if 0.695 <= leaf["eta"] <= 0.7201:
            for lg in glob.glob(os.path.join(T.DROOT, leaf["leaf"], "m_*", "run.log")):
                t = [int(x) for x in re.findall(r"##RUN .*?\((\d+) s\)", open(lg, errors="ignore").read())]; secs += sum(t); n += len(t)
    k = RP.koa_speed(); c100 = secs / n * k
    print(f"\ncost per trajectory at N = 100 in the window: {secs / n:.1f} s (Mac, {n} A1v2 trajectories) x KOA factor {k:.3f} = {c100:.1f} core-s")
    print("\n### Cost: 7 densities x 9 masses x seeds, per N (core-hours)\n")
    print("| N | seeds per mass | trajectories | core-s per trajectory, p = 2.87 | p = 2.37 | (a) current engine, p = 2.87 | (a) p = 2.37 | "
          "(b) 10x faster, p = 2.87 | (b) p = 2.37 |\n|---|---|---|---|---|---|---|---|---|")
    tot = {}
    for ns in NS:
        N = 2 * ns; ntr = 7 * 9 * seeds[N]; row = []
        for p in (2.87, 2.37):
            cs = c100 * (N / 100) ** p * (N / 100) ** 0.5; row.append(cs)
        ch = [ntr * cs / 3600 for cs in row]
        for p, v in zip((2.87, 2.37), ch): tot[p] = tot.get(p, 0) + v
        print(f"| {N} | {seeds[N]} | {ntr} | {row[0]:.0f} | {row[1]:.0f} | {ch[0]:.0f} | {ch[1]:.0f} | {ch[0] / 10:.0f} | {ch[1] / 10:.0f} |")
    print(f"\ntotal (a) current engine: {tot[2.87]:.0f} core-h (p = 2.87), {tot[2.37]:.0f} core-h (p = 2.37); "
          f"(b) 10x faster: {tot[2.87] / 10:.0f} / {tot[2.37] / 10:.0f} core-h")


if __name__ == "__main__":
    main()
