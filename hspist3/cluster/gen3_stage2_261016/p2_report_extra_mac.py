#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 57, addendum 2, task 1; notes 261012 sec. 4.7.138): THE TABLES ADDED TO REPORT D57-2, from the
fetched tree (hspist3/experiments_gen3_stage2_koa_261016/P2) and the folder of p2_report_mac.py (--report, for its p2_strain.csv). No
interpretation: the numbers and the rules stated here. Writes into a NEW folder --out:
  A. the time per CPU model: the CPU model of each run from its task log (job_logs/gen3-p2_<job>_<chunk>.out, "chunk <j> on <node> (<CPU
     model>)"; chunk j = task line j), the run's seconds and executed events from its cell's run.log (p2_runs.csv: the worker's clock from
     the binary's start to the run's record, whole seconds); per model: runs, minutes per run (median, range; and at alpha 20), and
     microseconds per executed event (median, range) against the D56-1 timing model (3.533 on E5-2670 v3, 4.242 on E5-4640 v2); the
     core-hours used against the package's prediction;
  B. the executed events against the package's prediction (gen_p2_koa_261016.plan(): r = 1.2004 (Z_upper - 1) + 0.0112 events per disk and
     sigma-time, Z_upper = Liu's solid branch x 1.268, times N (T_eq + T_rec) and the rehearsal's factor per pin, 1.171 for pin 0 and 1.002
     for turned), per cell (median, range; the medians per mass), with the cell's measured Z (the hold's last quarter, s1_runs.csv);
  C. the orientation per cell: the close-packed line direction per compartment (p2_lattice.config_strain: x or y) in the initial state, the
     hold's last snapshot and the last snapshot; the record's psi6 orientation (psi6_t, circular mean mod 60) within 7.5 deg of 0 or of 30;
     per compartment whose line direction changes, the first snapshot with the other direction (p2_frames.csv; the direction from theta0:
     y when theta0 mod 60 lies within 15 deg of 30, as config_strain decides) and whether it lies in the hold or the record;
  D. for every (cell, mass, compartment) of Report D57-2's test (iv) with R > 2: the run farthest from the others' median end strain, its
     end strain, its record mean and its strain at the hold's end, and the range of the other runs.
Writes p2_runs_cpu.csv (one row per run), p2_extra.md.
usage (from hspist3/): python3 cluster/gen3_stage2_261016/p2_report_extra_mac.py --data <fetched P2 folder> --report <p2_report_mac.py's folder>
       --out <new folder>
"""
import glob, math, os, re, sys, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
for p in (HERE, os.path.join(CL, "gen3_stage1_261012"), os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import gen_p2_koa_261016 as GP          # its TASKS follows this script's --out (read at its import): the task list is read from HERE
warnings.filterwarnings("ignore")
ETAS, GEOMS, PINS, ALPHAS = GP.ETAS, GP.GEOMS, GP.PINS, GP.ALPHAS
CELLS = [(eta, G, pin) for eta in ETAS for G in GEOMS for pin in PINS]


def cid(eta, G, pin): return f"P2_N{GP.N:04d}_eta{eta:.6f}_{G}_pin{pin}"
def short(m): x = re.search(r"(EPYC \w+|E5-\d+ v\d|Gold \d+)", m or ""); return x.group(1) if x else (m or "?")
def fx(v, d=3): return "-" if v is None or not np.isfinite(v) else f"{v:.{d}f}"


def mr(v, d=3):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return "-" if not len(v) else f"{np.median(v):.{d}f} ({v.min():.{d}f} to {v.max():.{d}f})"


def line_of(theta):
    t = np.mod(np.asarray(theta, float), 60.0); return np.where(np.isnan(t), np.nan, np.where(np.abs(t - 30.0) < 15.0, 90.0, 0.0))


def main():
    av = sys.argv; data = os.path.abspath(av[av.index("--data") + 1]); rep = os.path.abspath(av[av.index("--report") + 1])
    out = os.path.abspath(av[av.index("--out") + 1])
    if os.path.exists(out): sys.exit(f"REFUSED: {out} exists (a new folder only)")
    G3 = os.path.join(data, "merged", "stage2", "gen3")
    R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(G3, "*", "m_*", "p2_runs.csv")))], ignore_index=True)
    TL = [l.split() for l in open(os.path.join(HERE, "tasks_p2_koa_261016.txt")) if l.strip()]; key = {(t[1].split("/")[2], int(t[2]), int(t[3])): j for j, t in enumerate(TL, 1)}
    cpu, node, logs = {}, {}, sorted(glob.glob(os.path.join(data, "job_logs", "gen3-p2_*_*.out")))
    for f in logs:
        for l in open(f, errors="replace"):
            m = re.match(r"chunk (\d+) on (\S+) \((.*)\): group", l)
            if m: cpu[int(m.group(1))] = m.group(3); node[int(m.group(1))] = m.group(2)
    R["chunk"] = [key.get((c, int(M), int(r)), -1) for c, M, r in zip(R["cell"], R["M"], R["run"])]
    R["cpu_model"] = R["chunk"].map(cpu); R["cpu"] = R["cpu_model"].map(short); R["node_task_log"] = R["chunk"].map(node)
    R["seed_task_line"] = [int(TL[j - 1][4]) if j > 0 else -1 for j in R["chunk"]]
    cells, rows, _ = GP.plan(); pred = {(r["cid"], r["M"]): r for r in rows}
    R["events_pred"] = [pred[(c, int(M))]["ev"] for c, M in zip(R["cell"], R["M"])]; R["Z_upper_model"] = [pred[(c, int(M))]["Zu"] for c, M in zip(R["cell"], R["M"])]
    R["events_over_pred"] = R["events"] / R["events_pred"]; R["us_per_event"] = 1e6 * R["seconds"] / R["events"]; R["alpha"] = R["M"] / GP.N
    os.makedirs(out)
    R[["cell", "eta", "G", "pin", "M", "alpha", "run", "seed", "seed_task_line", "chunk", "node", "node_task_log", "cpu_model", "cpu", "seconds", "events",
       "events_pred", "events_over_pred", "us_per_event", "Z", "Z_upper_model"]].to_csv(os.path.join(out, "p2_runs_cpu.csv"), index=False, float_format="%.8g")
    L = ["# The tables added to Report D57-2 (p2_report_extra_mac.py; decision 57 addendum 2, task 1)", "",
         f"- runs: {len(R)}; task logs with a CPU line: {len(cpu)} of {len(logs)}; runs whose node in run.log equals their task log's: "
         f"{int((R['node'] == R['node_task_log']).sum())}; runs whose seed equals their task line's: {int((R['seed'] == R['seed_task_line']).sum())}", ""]
    L += ["## A. The time per CPU model", "",
          "Minutes per run: the worker's clock from the binary's start to the run's record (the binary and its file moves; run.log). Microseconds per "
          "executed event: those seconds over the run's executed events. The D56-1 timing model: 3.533 on E5-2670 v3, 4.242 on E5-4640 v2.", "",
          "| CPU model | runs | minutes per run: median (range) | minutes at alpha 20: median (largest) | microseconds per event: median (range) | "
          "median against E5-2670 v3's |", "|---|---|---|---|---|---|"]
    ref = float(R[R["cpu"] == "E5-2670 v3"]["us_per_event"].median()) if (R["cpu"] == "E5-2670 v3").any() else math.nan
    for c, g in sorted(R.groupby("cpu"), key=lambda it: it[1]["us_per_event"].median()):
        a20 = g[np.isclose(g["alpha"], 20)]["seconds"] / 60
        L.append(f"| {c} | {len(g)} | {mr(g['seconds'] / 60, 2)} | {fx(a20.median(), 2)} ({fx(a20.max(), 2)}) | {mr(g['us_per_event'], 3)} | "
                 f"{fx(g['us_per_event'].median() / ref, 3)} |")
    pv3 = sum(r["min_v3"] for r in rows) * GP.NSEED / 60; pv2 = sum(r["min_v2"] for r in rows) * GP.NSEED / 60
    L += ["", f"- core-hours used (the runs' seconds): {R['seconds'].sum() / 3600:.2f}; the package's prediction {pv3:.1f} on E5-2670 v3, {pv2:.1f} on E5-4640 v2; "
              f"the longest run {R['seconds'].max() / 60:.2f} min ({R.loc[R['seconds'].idxmax(), 'cpu']}, {R.loc[R['seconds'].idxmax(), 'cell'].replace('P2_N0900_', '')}, "
              f"alpha {R.loc[R['seconds'].idxmax(), 'alpha']:g}) against the array's limit of 25 min", ""]
    L += ["## B. The executed events against the package's prediction, per cell", "",
          "The prediction: gen_p2_koa_261016.plan() (r = 1.2004 (Z_upper - 1) + 0.0112 per disk and sigma-time with Z_upper = Liu's solid branch x 1.268, "
          "times N (T_eq + T_rec) and the rehearsal's factor per pin: 1.171 for pin 0, 1.002 for turned). Z measured: the hold's last quarter (mean over the runs).", "",
          "| eta | G | pin | Z_upper (model) | Z measured | executed / predicted: median (range) | median per mass (alpha 1, 2, 5, 10, 20) |", "|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        g = R[R["cell"] == cid(eta, G, pin)]
        L.append(f"| {eta} | {G} | {pin} | {fx(g['Z_upper_model'].mean(), 2)} | {fx(g['Z'].mean(), 2)} | {mr(g['events_over_pred'], 3)} | "
                 + ", ".join(fx(g[np.isclose(g['alpha'], a)]['events_over_pred'].median(), 3) for a in ALPHAS) + " |")
    F = []
    for f in sorted(glob.glob(os.path.join(G3, "*", "m_*", "p2_frames.csv"))):
        x = pd.read_csv(f); x["cell"] = f.split(os.sep)[-3]; x["M"] = int(f.split(os.sep)[-2][2:]); F.append(x)
    F = pd.concat(F, ignore_index=True) if F else pd.DataFrame()
    L += ["", "## C. The orientation per cell", "",
          "The close-packed line direction per compartment (x or y; p2_lattice.config_strain) in the initial state, the hold's last snapshot and the "
          "last snapshot, counted over the cell's 100 compartments (50 runs x 2); the record's psi6 orientation (psi6_t, the circular mean mod 60) "
          "within 7.5 deg of 0 (lines along x) or of 30 (along y); for the compartments whose direction changes, the first snapshot with the other "
          "direction (the direction from theta0: y when theta0 mod 60 lies within 15 deg of 30), its time from the start [sigma-time] and its phase "
          "(the hold lasts 5000).", "",
          "| eta | G | pin | start: x / y | the hold's end: x / y | the end: x / y | record orientation: near 0 / near 30 | compartments that change | "
          "the first snapshot with the other direction: t, median (range) | in the hold / in the record |", "|---|---|---|---|---|---|---|---|---|---|"]
    for eta, G, pin in CELLS:
        c = cid(eta, G, pin); g = R[R["cell"] == c]; cnt = []
        for ph in ("start", "hold", "end"):
            v = g[[f"{ph}_line_0", f"{ph}_line_1"]].to_numpy(float).ravel(); cnt.append(f"{int((v == 0).sum())} / {int((v == 90).sum())}")
        od = g[["orient_rec_0", "orient_rec_1"]].to_numpy(float).ravel()
        n0, n30 = int((np.abs(((od + 30.0) % 60.0) - 30.0) < 7.5).sum()), int((np.abs(od - 30.0) < 7.5).sum())
        ts, ph_ = [], []
        fr = F[F["cell"] == c] if len(F) else F
        for r in g.to_dict("records"):
            h = fr[(fr["M"] == r["M"]) & (fr["run"] == r["run"])].sort_values("frame") if len(fr) else fr
            for k in (0, 1):
                if not len(h) or not np.isfinite(r[f"start_line_{k}"]): continue
                ln = line_of(h[f"theta0_{k}"]); ch = np.where(np.isfinite(ln) & (ln != r[f"start_line_{k}"]))[0]
                if len(ch): ts.append(float(h["t"].iloc[ch[0]])); ph_.append(int(h["phase"].iloc[ch[0]]))
        L.append(f"| {eta} | {G} | {pin} | {cnt[0]} | {cnt[1]} | {cnt[2]} | {n0} / {n30} | {len(ts)} of {2 * len(g)} | {mr(ts, 0) if ts else '-'} | "
                 f"{sum(p == 0 for p in ph_)} / {sum(p == 1 for p in ph_)} |")
    ST = pd.read_csv(os.path.join(rep, "p2_strain.csv")); ST = ST[ST["alpha"].astype(str) != "all"]; D = []
    for r in ST.to_dict("records"):
        for k in (0, 1):
            if np.isfinite(r[f"R_{k}"]) and r[f"R_{k}"] > 2.0:
                g = R[(R["cell"] == r["cell"]) & np.isclose(R["alpha"], float(r["alpha"]))].copy(); med = g[f"end_e_{k}"].median()
                i = (g[f"end_e_{k}"] - med).abs().idxmax(); o = g.drop(index=i); x = g.loc[i]
                fr = F[(F["cell"] == r["cell"]) & (F["M"] == int(x["M"])) & (F["run"] == int(x["run"])) & (F["phase"] == 1)] if len(F) else F
                D.append(f"| {r['cell'].replace('P2_N0900_', '')} | {float(r['alpha']):g} | {'left' if k == 0 else 'right'} | {fx(r[f'R_{k}'], 2)} | {int(x['run'])} "
                         f"({int(x['seed'])}) | {fx(100 * x[f'end_e_{k}'], 3)} | {fx(100 * fr[f'e_{k}'].mean() if len(fr) else math.nan, 3)} | "
                         f"{fx(100 * x[f'hold_e_{k}'], 3)} | {fx(100 * o[f'end_e_{k}'].min(), 3)} to {fx(100 * o[f'end_e_{k}'].max(), 3)} |")
    L += ["", "## D. The entries of test (iv) with R > 2", "",
          "Per (cell, mass, compartment) with R > 2: the run farthest from the median end strain of its 10, its end strain, its mean over the record's "
          "snapshots, its strain at the hold's end, and the end strains of the other 9 (strains in %).", ""]
    L += (["| cell | alpha | compartment | R | run (seed) | end e | record mean e | e at the hold's end | the other 9: end e |", "|---|---|---|---|---|---|---|---|---|"] + D
          if D else ["none"])
    open(os.path.join(out, "p2_extra.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
