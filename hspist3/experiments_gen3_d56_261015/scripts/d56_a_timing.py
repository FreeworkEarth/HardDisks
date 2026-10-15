#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56, task 2 a; notes 261012 sec. 4.7.133): THE WALL TIME PER RUN AGAINST ETA, for the Stage 2
cells (eta 0.69, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80 at N 900 and 1600). Planning numbers, nothing registered.
  1. Stage 1's runs of B900 and B1600 (all nine masses, 25 seeds, 16 cells each: 7200 runs), from the fetched tables: the seconds and
     the node per run (runs_health.csv), the executed events (run.log), the measured Z (s1_runs.csv), the binary's nu_pred (red_nu.csv),
     the hold (run.log); the CPU per node from the timing tables (experiments_gen3_s1_261012/interim_261013/timing/timing_B*.tsv).
     The run's length T = T_eq + T_rec, T_eq = hold steps x 0.4/24, T_rec = ceil((200 / nu_pred) / DTF) DTF (the binary's own rule).
  2. The event rate per disk and sigma-time, r = events / (N T), fitted against the Enskog collision rate: r = a (Z - 1) + b with
     Z - 1 = 2 eta chi(eta) (chi the contact value): (i) Z the run's own measured Z; (ii) Z from KR (eta chi of the reference EOS, the
     form the plan author asked for). Least squares over the runs. The time per event, microseconds, per (N, CPU model), median over
     the dense cells (eta >= 0.63).
  3. Stage 2, per (eta, N): Z from the EOS (KR to 0.69; Liu's global branch at 0.70, inside the coexistence; Liu's solid branch from
     0.72: planning only, the solid's Z in these boxes is not the bulk one), and an upper planning value: the solid cells' Z times the
     pilot's measured excess at eta 0.76 (Z_pilot / Z_Liu, d56_c_analysis.py); the run length per mass by stage 1's rules (T_eq = 5000,
     200 periods of the binary's Henderson nu_pred, which runs high in the solid: 200 predicted periods are about 192 measured ones at
     0.76); events per run (alpha 5, the mean over the nine masses); minutes per run on the two CPU models of B900 and B1600;
     core-hours for 9 masses x 25 seeds and for 4 masses x 25 seeds (the middle four, alpha 3, 5, 7.5, 10); the merge time from B1600's
     record (4:44:18 by its files, 4:44:19 sacct): scaled by the run count, and by its dense cells' own per-run cost.
Writes into experiments_gen3_d56_261015/a_timing/ (a new folder): a_stage1_runs.csv, a_fits.txt, a_stage2_table.csv, a_stage2_table.md,
fig_a1_rate_and_time_per_event.png.
usage (from hspist3/): python3 experiments_gen3_d56_261015/scripts/d56_a_timing.py [--out <a new folder instead>]
"""
import glob, math, os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d56_261015")
for p in (os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_stage1_koa_261012 as GEN
import plot_speed_of_sound_edmd as SOS
KOA = os.path.join(HS, "experiments_gen3_stage1_koa_261012"); TIMING = os.path.join(HS, "experiments_gen3_s1_261012", "interim_261013", "timing")
MERGE = os.path.join(HS, "experiments_gen3_s1_261012", "interim_261013", "B1600", "e_merge_time.txt")
PILOT = os.path.join(E, "c_analysis", "c_runs.csv")
OUT = os.path.join(E, "a_timing")
ETAS, NS, MID4 = (0.69, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80), (900, 1600), (3, 5, 7.5, 10)
V3, V2 = "Intel(R) Xeon(R) CPU E5-2670 v3 @ 2.30GHz", "Intel(R) Xeon(R) CPU E5-4640 v2 @ 2.20GHz"


def t_rec(nup): return math.ceil((200.0 / nup) / GEN.DTF) * GEN.DTF


def z_eos(eta):
    a = np.array([eta])
    if eta <= 0.69 + 1e-12: return float(SOS.Z_kolafa_rottner_2006(a)[0]), "KR"
    if eta < 0.72 - 1e-12: return float(SOS.Z_liu_global(a)[0]), "Liu global (coexistence)"
    return float(SOS.Z_liu_solid(a)[0]), "Liu solid"


def stage1_runs():
    tm = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(TIMING, "timing_B*.tsv")))])
    cpu = dict(zip(tm["node"], tm["cpu"])); rows = []
    for g in ("B900", "B1600"):
        for m in sorted(glob.glob(os.path.join(KOA, g, "merged", "stage1", "gen3", "*", "m_*"))):
            cell = m.split(os.sep)[-2]; M = int(m.split("_")[-1]); N = int(cell[3:7]); eta = float(cell.split("eta")[1])
            h = pd.read_csv(os.path.join(m, "runs_health.csv")); s1 = pd.read_csv(os.path.join(m, "s1_runs.csv")); red = pd.read_csv(os.path.join(m, "red_nu.csv"))
            log = open(os.path.join(m, "run.log"), errors="replace").read()
            ev = {int(a): int(b) for a, b in re.findall(r"^##RUN \S+ \S+ run (\d+) seed \d+ .*?executed events (\d+)", log, re.M | re.S)}
            hold = {int(a): int(b) for a, b in re.findall(r"^##RUN \S+ \S+ run (\d+) seed \d+ .*?^Wall hold steps: (\d+)$", log, re.M | re.S)}
            for r in h.itertuples():
                rr = red[red["run"] == r.run]; ss = s1[s1["run"] == r.run]
                if not len(rr) or r.run not in ev: continue
                nup = float(rr["nu_pred"].iloc[0]); Teq = hold[r.run] * 0.4 / 24.0; Tr = t_rec(nup)
                rows.append(dict(group=g, cell=cell, N=N, eta=eta, M=M, run=int(r.run), seconds=int(r.seconds), node=r.node, cpu=cpu.get(r.node, "?"),
                                 events=ev[r.run], nu_pred=nup, T_eq=Teq, T_rec=Tr, T=Teq + Tr, Z=float(ss["Z"].iloc[0]) if len(ss) else math.nan,
                                 Z_kr=z_eos(eta)[0] if eta <= 0.69 + 1e-12 else math.nan))
    D = pd.DataFrame(rows); D["rate"] = D["events"] / (D["N"] * D["T"]); D["us_per_event"] = 1e6 * D["seconds"] / D["events"]
    return D


def lsq(x, y):
    A = np.vstack([x, np.ones_like(x)]).T; (a, b), *_ = np.linalg.lstsq(A, y, rcond=None); res = y - (a * x + b)
    return a, b, float(np.sqrt(np.mean((res / y) ** 2))), float(np.max(np.abs(res / y)))


def merge_record():
    txt = open(MERGE).read(); per = {}
    for m in re.finditer(r"^\| (B_N1600_eta[0-9.]+) \| 9 \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \|$", txt, re.M):
        per[float(m.group(1).split("eta")[1])] = int(m.group(6))
    scan = int(re.search(r"^\| scan .*? \| (\d+) \| \d+:\d+:\d+ \|", txt, re.M).group(1)); moves = int(re.search(r"^\| moves .*? \| (\d+) \| \d+:\d+:\d+ \|", txt, re.M).group(1))
    job = int(re.search(r"^\| \*\*the job\*\* .*? \| (\d+) \| ", txt, re.M).group(1))
    return per, scan, moves, job


def hms(s): s = int(round(s)); return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def main():
    global OUT
    if "--out" in sys.argv: OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    os.makedirs(OUT); D = stage1_runs(); D.to_csv(os.path.join(OUT, "a_stage1_runs.csv"), index=False, float_format="%.10g")
    L = [f"# d56_a_timing.py: the wall time per run against eta (decision 56 task 2 a)", "",
         f"stage 1 runs read: {len(D)} (B900 {sum(D.group == 'B900')}, B1600 {sum(D.group == 'B1600')}); CPU models: "
         + ", ".join(f"{k} {v}" for k, v in D["cpu"].value_counts().items())]
    a1, b1, rms1, mx1 = lsq(D["Z"].to_numpy() - 1, D["rate"].to_numpy())
    a2, b2, rms2, mx2 = lsq(D["Z_kr"].to_numpy() - 1, D["rate"].to_numpy())
    L += ["", "## The event rate per disk and sigma-time, r = events / (N T), against the collision rate (Z - 1 = 2 eta chi)", "",
          f"- (i) Z measured per run: r = {a1:.4f} (Z - 1) + {b1:.4f}; relative residual rms {100 * rms1:.2f} %, largest {100 * mx1:.2f} %",
          f"- (ii) Z from KR (eta chi of the reference EOS): r = {a2:.4f} (Z_KR - 1) + {b2:.4f}; relative residual rms {100 * rms2:.2f} %, largest {100 * mx2:.2f} %",
          f"- the stage 1 generator's planning rule (gen_stage1_koa_261012.events_per_sigma): r = 1.128 (Z_KR - 1) + 1.2", "",
          "per cell (the mean over its 225 runs): eta, N, measured Z, r, the residual of (i) and of (ii):", ""]
    pc = D.groupby(["N", "eta"]).agg(Z=("Z", "mean"), Z_kr=("Z_kr", "mean"), rate=("rate", "mean")).reset_index()
    for r in pc.itertuples():
        L.append(f"  {r.eta:.4f} {r.N:5d}  Z {r.Z:8.4f}  r {r.rate:8.4f}  (i) {100 * (a1 * (r.Z - 1) + b1 - r.rate) / r.rate:+6.2f} %  (ii) {100 * (a2 * (r.Z_kr - 1) + b2 - r.rate) / r.rate:+6.2f} %")
    dense = D[D["eta"] >= 0.63 - 1e-9]; tev = dense.groupby(["N", "cpu"])["us_per_event"].median()
    tq = dense.groupby(["N", "cpu"])["us_per_event"].quantile([0.1, 0.9]).unstack()
    L += ["", "## Microseconds per event (the dense cells, eta >= 0.63; median, 10 and 90 % quantiles)", ""]
    for (N, c), v in tev.items():
        L.append(f"- N {N}, {c}: {v:.3f} ({tq.loc[(N, c), 0.1]:.3f} .. {tq.loc[(N, c), 0.9]:.3f}), runs {len(dense[(dense.N == N) & (dense.cpu == c)])}")
    pil = pd.read_csv(PILOT); p76 = pil[pil["cell"].astype(str) == "76"]
    T76 = 5000.0 + t_rec(float(p76["nu_pred"].iloc[0])); r76 = float((p76["events"] / (900 * T76)).mean()); Z76 = float(p76["Z"].mean())
    zl76 = z_eos(0.76)[0]; exc = Z76 / zl76
    L += ["", "## The pilot (N 900, eta 0.76, alpha 5, the Mac; d56_c_analysis.py)", "",
          f"- measured: Z {Z76:.3f} (Liu solid {zl76:.3f}: x {exc:.3f}); r {r76:.4f} per disk and sigma-time; the fit (i) with the measured Z: "
          f"{a1 * (Z76 - 1) + b1:.4f} ({100 * (a1 * (Z76 - 1) + b1 - r76) / r76:+.2f} %); with Liu's Z: {a1 * (zl76 - 1) + b1:.4f} ({100 * (a1 * (zl76 - 1) + b1 - r76) / r76:+.2f} %)",
          f"- the Mac: {float((p76['mac_seconds'] / p76['events']).mean()) * 1e6:.3f} microseconds per event (five runs at once)"]
    per, scan, moves, job = merge_record(); n16 = 16 * 9 * 25
    dense_cells = {e: s for e, s in per.items() if e >= 0.60 - 1e-9}; per_run_dense = (np.mean(list(dense_cells.values())) / 225.0) + (scan + moves) / n16
    rows = []
    for N in NS:
        H = 10.0 * math.sqrt(N / 100.0)
        for eta in ETAS:
            L0 = N * math.pi / (8.0 * H * eta); Z, src = z_eos(eta); Zhi = Z * exc if eta >= 0.72 - 1e-12 else Z
            r_lo, r_hi = a1 * (Z - 1) + b1, a1 * (Zhi - 1) + b1
            per_a = []
            for al in GEN.ALPHA:
                nup = GEN.nu_pred_binary(N, H, L0, al * N); T = 5000.0 + t_rec(nup); per_a.append((al, nup, T))
            Ts = np.array([x[2] for x in per_a]); T5 = [x[2] for x in per_a if x[0] == 5][0]; Tm4 = np.array([x[2] for x in per_a if x[0] in MID4])
            row = dict(N=N, eta=eta, H=H, L0=L0, Z_eos=Z, Z_source=src, Z_plan_hi=Zhi, rate_eos=r_lo, rate_hi=r_hi, T_alpha5=T5, T_mean9=float(Ts.mean()),
                       T_min9=float(Ts.min()), T_max9=float(Ts.max()))
            for lab, r in (("eos", r_lo), ("hi", r_hi)):
                ev5 = r * N * T5; ev9 = r * N * Ts; ev4 = r * N * Tm4
                row[f"events_alpha5_{lab}"] = ev5; row[f"events_mean9_{lab}"] = float(ev9.mean())
                for cl, c in (("v3", V3), ("v2", V2)):
                    us = float(tev.loc[(N, c)])
                    row[f"min_alpha5_{lab}_{cl}"] = ev5 * us / 6e7; row[f"min_mean9_{lab}_{cl}"] = float(ev9.mean()) * us / 6e7
                    row[f"coreh_9x25_{lab}_{cl}"] = float(ev9.sum()) * 25 * us / 3.6e9; row[f"coreh_4x25_{lab}_{cl}"] = float(ev4.sum()) * 25 * us / 3.6e9
            rows.append(row)
    S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "a_stage2_table.csv"), index=False, float_format="%.6g")
    md = ["| eta | N | Z (EOS) | Z (upper) | events per disk and sigma-time | run length alpha 5 (mean of 9) [sigma-time] | events per run, alpha 5 | "
          "minutes per run, alpha 5 (v3 / v2) | minutes per run, mean of 9 (v3) | core-hours 9 x 25 (v3) | core-hours 4 x 25 (v3) |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in S.itertuples():
        md.append(f"| {r.eta:.2f} | {r.N} | {r.Z_eos:.2f} ({r.Z_source}) | {r.Z_plan_hi:.2f} | {r.rate_eos:.2f} .. {r.rate_hi:.2f} | {r.T_alpha5:.0f} ({r.T_mean9:.0f}) | "
                  f"{r.events_alpha5_eos / 1e6:.0f} .. {r.events_alpha5_hi / 1e6:.0f} M | {r.min_alpha5_eos_v3:.1f} .. {r.min_alpha5_hi_v3:.1f} / {r.min_alpha5_eos_v2:.1f} .. {r.min_alpha5_hi_v2:.1f} | "
                  f"{r.min_mean9_eos_v3:.1f} .. {r.min_mean9_hi_v3:.1f} | {r.coreh_9x25_eos_v3:.0f} .. {r.coreh_9x25_hi_v3:.0f} | {r.coreh_4x25_eos_v3:.0f} .. {r.coreh_4x25_hi_v3:.0f} |")
    for N in NS:
        s = S[S.N == N]
        md.append(f"| all 7 | {N} | | | | | | | | {s.coreh_9x25_eos_v3.sum():.0f} .. {s.coreh_9x25_hi_v3.sum():.0f} | {s.coreh_4x25_eos_v3.sum():.0f} .. {s.coreh_4x25_hi_v3.sum():.0f} |")
    md.append(f"| all 14 | | | | | | | | | {S.coreh_9x25_eos_v3.sum():.0f} .. {S.coreh_9x25_hi_v3.sum():.0f} (v2: {S.coreh_9x25_eos_v2.sum():.0f} .. {S.coreh_9x25_hi_v2.sum():.0f}) | "
              f"{S.coreh_4x25_eos_v3.sum():.0f} .. {S.coreh_4x25_hi_v3.sum():.0f} (v2: {S.coreh_4x25_eos_v2.sum():.0f} .. {S.coreh_4x25_hi_v2.sum():.0f}) |")
    mg = ["", "## The merge, from B1600's record (scan, moves and per cell: experiments_gen3_s1_261012/interim_261013/B1600/e_merge_time.txt)", "",
          f"- B1600: the job {hms(job)} for {n16} runs (16 cells x 9 masses x 25 seeds); the scan and the moves {scan + moves} s; the cells with eta <= 0.25 "
          f"took {sum(s for e, s in per.items() if e <= 0.25 + 1e-9)} s of the {sum(per.values())} s of per-cell work, the six cells eta 0.60-0.69 {sum(dense_cells.values())} s "
          f"({min(dense_cells.values())}-{max(dense_cells.values())} s each, 225 runs)",
          f"- per run, from the dense cells: {per_run_dense:.3f} s (their mean per-cell work / 225 + the scan and the moves per run)"]
    for N in NS:
        for nm, nmass in (("9 masses", 9), ("4 masses", 4)):
            n = 7 * nmass * 25
            mg.append(f"- Stage 2, N {N}, 7 cells x {nm} x 25 seeds = {n} runs: by the run count {hms(job * n / n16)}; by the dense cells' cost {hms(per_run_dense * n)}")
    L += ["", "## Stage 2 (planning; minutes and core-hours at the dense cells' median microseconds per event; v3 = E5-2670 v3, v2 = E5-4640 v2)", "",
          "The lower value of each range: Z from the EOS (KR, Liu); the upper: the solid cells' Z times the pilot's measured excess "
          f"(x {exc:.3f} at 0.76). The run length: T_eq 5000 + 200 periods of the binary's nu_pred per mass; 4 masses: alpha 3, 5, 7.5, 10.", ""] + md + mg
    open(os.path.join(OUT, "a_stage2_table.md"), "w").write("\n".join(md + mg) + "\n")
    open(os.path.join(OUT, "a_fits.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15, 5.6))
    for N, mk, col, ms in ((1600, "D", "#0b2545", 5), (900, "^", "#2a6fbb", 10)):
        p = pc[pc.N == N]; ax.plot(p.eta, p.rate, mk, color=col, ms=ms, mfc="none" if N == 900 else col, mew=1.3,
                                   label=f"stage 1, N {N} (the mean over 225 runs per cell)")
    e = np.linspace(0.02, 0.80, 300); zz = np.array([z_eos(x)[0] for x in e])
    ax.plot(e, a1 * (zz - 1) + b1, color="#c0392b", lw=1.5, label=f"r = {a1:.3f} (Z - 1) + {b1:.3f}, Z from KR to 0.69, Liu beyond")
    zh = np.array([z_eos(x)[0] * (exc if x >= 0.72 else 1.0) for x in e]); ax.plot(e, a1 * (zh - 1) + b1, color="#c0392b", lw=1.0, ls="--", label=f"the solid's Z x {exc:.3f} (the pilot's excess)")
    ax.plot([0.76], [r76], "*", color="#d9822b", ms=14, mec="k", label="the pilot (N 900, eta 0.76; Mac, 2 runs)")
    for x in ETAS: ax.axvline(x, color="0.85", lw=0.6, zorder=0)
    ax.set_xlabel("packing fraction eta"); ax.set_ylabel("events per disk and sigma-time"); ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="upper left")
    ax.set_title("the event rate against the collision rate (Z - 1 = 2 eta chi); Stage 2's etas as grey lines", fontsize=10)
    for (N, c), mk in (((900, V3), "^"), ((1600, V3), "D"), ((900, V2), "v"), ((1600, V2), "s")):
        q = D[(D.N == N) & (D.cpu == c)].groupby("eta")["us_per_event"].median()
        bx.plot(q.index, q.values, mk + "-", lw=0.8, ms=5, label=f"N {N}, {'E5-2670 v3' if c == V3 else 'E5-4640 v2'}")
    bx.set_xlabel("packing fraction eta"); bx.set_ylabel("microseconds per event (median per cell)"); bx.grid(alpha=0.25); bx.legend(fontsize=8)
    bx.set_title("KOA's time per event, stage 1 B900 and B1600", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig_a1_rate_and_time_per_event.png"), dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
