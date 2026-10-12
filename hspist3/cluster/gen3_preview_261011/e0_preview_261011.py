#!/usr/bin/env python3
"""##CHRIS 2026-10-11 (plan-author decision 29, task P; 261012 sec. 4.7.59): THE REDUCTION OF THE PREVIEW SWEEP -- PREVIEW, INFORMATION;
no verdict, never part of stage 1, never in a paper. The canonical estimator E0 of validation/paper1_populate_cs_err_20261002.py; its
functions are imported unchanged, and its arithmetic is main()'s:
  per trajectory  the frequency of E0's cell(): TD = 200 periods, X_EDGE = 2.5, T._load/_prefix/_spectrum, the argmax above
                  k = TD / X_EDGE. On KOA the merge runs reduce_B.py, which computes exactly this per trajectory (it imports TD and X_EDGE
                  from paper1_populate_cs_err; its docstring's gate). The Mac therefore reads red_nu.csv, and no trace is fetched.
  per (cell, M)   the mean frequency over the seeds, its SD (ddof 1) and n, as cell() returns them.
  per cell        the box-truncation correction of methods sec. 14, as main() writes the canonical table:
                    delta = box_delta(L0);
                    eta = eta_rec L0 / (L0 - delta/2);
                    L_eff,true = L0 - 2r - t/2 - delta/2 and x_M = K(alpha) / (2 pi L_eff,true);
                    c_s, c_s_err, chi2_red and c_s_err_scaled from slope_with_errors(x_M, mean nu, SD/sqrt(n)), with the slope and
                    c_s_scatter_mass from T.slope (the canonical central value);
                    KR and dev_KR_pct at the true eta where eta_rec <= 0.69.
  recorded eta and lengths per cell (methods sec. 14) come from the binary's own [EDMD3-RUN] lines in run.log: L0, eta, the box in
                  pixels and exact_box. They are checked against eta_rec = N pi r^2 / (2 L0 H) and against box_delta(L0), which reproduces
                  SIM_WIDTH = (int)(2 L0 24) in float32. eta_rec is the recorded eta rounded to 6 decimals, the precision of the trace's
                  eta column that E0 reads for A1 v2.
GATE (--gate-a1v2): this module's per-cell code is fed A1 v2's own trajectories. It must reproduce, string for string, E0's 13
printed columns of the canonical table 260919_A1v2_final_cs_vs_eta.csv in all its cells: eta, L0, c_s, c_s_err, c_s_err_scaled,
c_s_scatter_mass, chi2_red, n_masses, trajectories_used, trajectories_discarded, eta_rec, delta_sigma, L_eff_true. The per-trajectory
frequencies come from paper1_populate_cs_err.cell() itself, on every leaf and mass of DROOT, as main() builds its tasks.
The theory columns KR and dev_KR_pct are compared too, and printed as INFORMATION with their largest difference; they are not part of
the gate.
  WHY THEY ARE NOT: the first run of this gate (all 15 columns) reproduced 27 rows in full. The other 8 differed only in KR (one unit
  in the 5th decimal) and dev_KR_pct (one unit in the 3rd). kr_cs is unchanged since b9f2ef86 (2026-09-19). The KR function itself
  entered git only at 4f5fd314 (2026-10-02 19:08), after the table's last write (e823187a, 15:46 the same day), so its state at that
  write is not recorded.
usage (from hspist3/):
  python3 cluster/gen3_preview_261011/e0_preview_261011.py --gate-a1v2 [--jobs 6]
  python3 cluster/gen3_preview_261011/e0_preview_261011.py --preview <fetched merged/ folder> --tasks <task list> --out <csv>
"""
import argparse, csv, glob, math, os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import tests_20260913 as T
import paper1_populate_cs_err_20261002 as E0

FIELDS = ["eta", "L0", "c_s", "c_s_err", "c_s_err_scaled", "c_s_scatter_mass", "chi2_red", "n_masses", "trajectories_used",
          "trajectories_discarded", "KR", "dev_KR_pct", "eta_rec", "delta_sigma", "L_eff_true"]   # the canonical table's columns
THEORY = ("KR", "dev_KR_pct")                                                                  # information, not the gate
E0_FIELDS = [f for f in FIELDS if f not in THEORY]
RUN_RE = re.compile(r"^\[EDMD3-RUN\] L0=(\S+) M=(\d+) run=\d+ seed=(\d+): .*? box=(\d+) x (\d+) px exact_box=(\d) eta=(\S+) N=(\d+) N_L=(\d+)",
                    re.M)


def e0_row(L0_str, eta_rec, cells):
    """the canonical table's row for one density from its per-(eta, M) cells [dict(M, nu, sd, n, nd)], by main()'s arithmetic"""
    L0 = float(L0_str)
    cs = sorted([c for c in cells if c["n"] > 0], key=lambda c: c["M"])
    if len(cs) < 3: return None
    y = np.array([q["nu"] for q in cs])
    sy = np.array([(q["sd"] / math.sqrt(q["n"])) if q["n"] > 1 else np.nan for q in cs])
    delta = E0.box_delta(L0)
    LeT = T.l_eff(L0) - delta / 2
    xT = np.array([T.k_root(q["M"] / (2.0 * T.N_SIDE)) / (2 * math.pi * LeT) for q in cs])
    sT, scatT, nm = T.slope(xT, y)
    s2, errT, err_scaled, chi2T = E0.slope_with_errors(xT, y, sy)
    if not (abs(s2 - sT) <= 1e-12 * abs(sT)): sys.exit(f"STOP: slope_with_errors' slope {s2} differs from T.slope's {sT}")
    etaT = eta_rec * L0 / (L0 - delta / 2)
    r = dict(eta=f"{etaT:.6f}", L0=L0_str, c_s=f"{sT:.5f}", c_s_err=f"{errT:.6f}", c_s_err_scaled=f"{err_scaled:.6f}",
             c_s_scatter_mass=f"{scatT:.5f}", chi2_red=f"{chi2T:.3f}", n_masses=str(nm), trajectories_used=str(sum(q["n"] for q in cs)),
             trajectories_discarded=str(sum(q["nd"] for q in cs)), eta_rec=f"{eta_rec:.6f}", delta_sigma=f"{delta:.6f}",
             L_eff_true=f"{LeT:.6f}", KR="", dev_KR_pct="")
    if eta_rec <= 0.69:
        kr = T.kr_cs(etaT); r["KR"] = f"{kr:.5f}"; r["dev_KR_pct"] = f"{100 * (sT / kr - 1):+.3f}"
    return r


# ---------------------------------------------------------------- the gate on A1 v2
def gate_a1v2(jobs):
    from multiprocessing import Pool
    table = T.a1_leaf_table(); tasks = []
    for leaf in table:
        for M in T.A1_MASSES:
            runs = T.cell_runs(os.path.join(T.DROOT, leaf["leaf"], f"m_{M}"), M)
            if runs: tasks.append((leaf["eta"], float(leaf["L0"]), M, runs))
    print(f"# GATE: e0_preview_261011.py's per-cell code on A1 v2's own trajectories against the canonical table\n\n"
          f"{len(tasks)} (eta, M) cells, {sum(len(t[3]) for t in tasks)} trajectories (DROOT {os.path.relpath(T.DROOT, HS)}); "
          f"per-trajectory frequencies by paper1_populate_cs_err.cell() itself, Pool({jobs})")
    with Pool(jobs) as pool: cells = pool.map(E0.cell, tasks, chunksize=1)
    canon = {r["eta_rec"]: r for r in csv.DictReader(open(T.plot_path(E0.SRC19)))}
    print(f"canonical table: {os.path.relpath(T.plot_path(E0.SRC19), os.path.dirname(HS))} ({len(canon)} rows)\n")
    print("| eta_rec | L0 | c_s | c_s_err_scaled | chi2_red | E0 columns equal (of 13) | E0 columns differing | KR / dev_KR_pct (information) |\n"
          "|---|---|---|---|---|---|---|---|")
    bad, n, kr_rel, kr_rows = [], 0, 0.0, 0
    for leaf in table:
        mine = e0_row(leaf["L0"], leaf["eta"], [dict(M=c["M"], nu=c["nu"], sd=c["sd"], n=c["n"], nd=c["nd"]) for c in cells
                                                   if c["eta"] == leaf["eta"]])
        if mine is None: continue
        n += 1; ref = canon.get(mine["eta_rec"])
        diff = [k for k in E0_FIELDS if ref is None or str(ref.get(k, "")) != mine[k]]
        if diff: bad.append((mine["eta_rec"], diff, {k: (ref.get(k) if ref else None, mine[k]) for k in diff}))
        th = [k for k in THEORY if ref is not None and str(ref.get(k, "")) != mine[k]]
        if ref is not None and ref.get("KR") and mine["KR"]:
            kr_rel = max(kr_rel, abs(float(ref["KR"]) / float(mine["KR"]) - 1))
        kr_rows += bool(th)
        thtxt = "equal" if not th else "; ".join(f"{k} {ref.get(k)} (table) / {mine[k]} (now)" for k in th)
        print(f"| {mine['eta_rec']} | {mine['L0']} | {mine['c_s']} | {mine['c_s_err_scaled']} | {mine['chi2_red']} | "
              f"{len(E0_FIELDS) - len(diff)} | {', '.join(diff) if diff else '-'} | {thtxt} |")
    for e, d, v in bad: print(f"  {e}: {v}")
    ok = not bad and n == len(canon)
    print(f"\nINFORMATION (not the gate): KR and dev_KR_pct differ in {kr_rows} of {n} rows; largest |KR table / KR now - 1| = {kr_rel:.1e}")
    print(f"GATE: {'PASS' if ok else 'FAIL'} -- {n - len(bad)} of {len(canon)} canonical rows reproduced in all 13 E0 columns"
          + ("" if n == len(canon) else f"; {n} rows computed"))
    return ok


# ---------------------------------------------------------------- the preview
def recorded(cell_dir, H, Ns):
    """the binary's [EDMD3-RUN] lines of one (cell, M): (L0, eta, box px, exact_box) -- one value per cell -- and the checks"""
    txt = open(os.path.join(cell_dir, "run.log"), errors="replace").read()
    recs = {(m.group(1), m.group(4), m.group(5), m.group(6), m.group(7), m.group(8), m.group(9)) for m in RUN_RE.finditer(txt)}
    if len(recs) != 1: sys.exit(f"STOP: {cell_dir}: {len(recs)} distinct [EDMD3-RUN] records (one expected)")
    L0s, bw, bh, ex, eta, N, NL = next(iter(recs)); L0 = float(L0s)
    eta_nom = 2 * int(Ns) * math.pi * 0.25 / (2 * L0 * float(H))
    delta = E0.box_delta(L0)
    checks = {"eta = N pi r^2 / (2 L0 H)": abs(float(eta) / eta_nom - 1) < 1e-12,
              "box width = (int)(2 L0 24) px, i.e. box_delta(L0)": abs((2 * L0 * 24 - int(bw)) / 24 - delta) < 1e-6,
              "exact_box = 0 (the A1 v2 command)": ex == "0", "N = 2 N_L = 2 N_s": int(N) == 2 * int(NL) == 2 * int(Ns)}
    return dict(L0=L0, eta_full=float(eta), box_px=f"{bw} x {bh}", exact_box=int(ex), delta=delta), checks


def preview_table(merged, tasks_path, out):
    L = [l.split() for l in open(tasks_path) if l.strip()]
    cellmeta = {}
    for p in L:
        rel = p[1]; cid = rel.split("/")[2]; cellmeta.setdefault(cid, dict(L0s=p[5], H=p[6], Ns=p[7], rel_cell="/".join(rel.split("/")[:3])))
    rows, recs = [], []
    for cid, m in cellmeta.items():
        cells, rec0 = [], None
        for d in sorted(glob.glob(os.path.join(merged, m["rel_cell"], "m_*")), key=lambda q: int(q.rsplit("_", 1)[1])):
            M = int(d.rsplit("_", 1)[1]); red = pd.read_csv(os.path.join(d, "red_nu.csv"))
            nu = red["nu"].to_numpy(float); rec, chk = recorded(d, m["H"], m["Ns"])
            nd = len(glob.glob(os.path.join(d, ".failed_run*")))
            cells.append(dict(M=M, nu=float(nu.mean()), sd=float(nu.std(ddof=1)) if len(nu) > 1 else float("nan"), n=len(nu), nd=nd))
            recs.append(dict(cell=cid, M=M, **rec, **{k: v for k, v in chk.items()}))
            if rec0 is None: rec0 = rec
            elif (rec["L0"], rec["eta_full"], rec["box_px"]) != (rec0["L0"], rec0["eta_full"], rec0["box_px"]):
                sys.exit(f"STOP: {cid}: the masses record different lengths")
        r = e0_row(m["L0s"], round(rec0["eta_full"], 6), cells)
        if r: r.update(cell=cid, eta_recorded_full=repr(rec0["eta_full"]), box_px=rec0["box_px"]); rows.append(r)
    pd.DataFrame(rows).to_csv(out, index=False)
    pd.DataFrame(recs).to_csv(out.replace(".csv", "_recorded.csv"), index=False)
    return rows, recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-a1v2", action="store_true"); ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--preview"); ap.add_argument("--tasks"); ap.add_argument("--out")
    a = ap.parse_args()
    if a.gate_a1v2: sys.exit(0 if gate_a1v2(a.jobs) else 1)
    if a.preview and a.tasks and a.out:
        rows, recs = preview_table(a.preview, a.tasks, a.out)
        print(f"{len(rows)} cells -> {a.out}; recorded lengths -> {a.out.replace('.csv', '_recorded.csv')}")
        sys.exit(0)
    ap.error("--gate-a1v2, or --preview with --tasks and --out")


if __name__ == "__main__":
    main()
