#!/usr/bin/env python3
"""##CHRIS 2026-10-12 (plan-author decision 39 item 5 f; registration 261012 sec. 4.8.6 and 4.8.8; gate G1): THE KOA REDUCTION OF ONE
STAGE 1 CELL FOLDER m_<M>, after reduce_B.py (merge_s1_koa.py runs both): the per-trajectory tables the Mac needs, so that the fetch
carries no traces (sec. 4.8.8: "reduce_B.py per (cell, mass); nu_d and the seed-averaged ACF; one summary row per trajectory for Z and
psi6; the snapshot reduction").
  nud_runs.csv, runs_health.csv, failed_runs.csv   cluster/gen3_testR_261011/reduce_testR_koa.py's three tables, imported unchanged:
                   nu_d per trajectory by the sec. 4.4.13 item 4a fit with nu_ref = the folder's own argmax mean over its clean
                   trajectories (decision 39 item 5 d, "as Test R"); the ##RUN health table; the refused runs.
  acf_mean_fit.csv the item-4a fit of the folder's seed-averaged ACF (the ACFs cut to the shortest, the same fit_acf and nu_ref):
                   A, tau_T, B, tau_r, omega, phi; nu_d_meanfit = omega / 2 pi; Q_meanfit = pi nu_d_meanfit tau_r; fit_ok.
  s1_runs.csv      one row per trajectory:
                   - Z (sec. 4.8.6.4): the mean over the hold's blocks in its last quarter of the two compartments' mean virial Z_c
                     (virial_blocks_*_run<s>.csv.gz); Z_q3 and Z_q4, the same over the third and the fourth quarter (sec. 4.8.7.6);
                     the blocks used per quarter;
                   - psi6 (sec. 4.8.6.4): the record means of the two compartments' psi6 (psi6_t_*_run<s>.csv.gz, columns psi6_0,
                     psi6_1, record rows), psi6_rec_0, psi6_rec_1 and their mean psi6_rec; the hold's last value of each (information);
                   - the halves (sec. 4.8.7.6): nu1, nu2 = the argmax estimator on the first and the second half of the analysed
                     prefix, exactly as validation/paper1_window_aging_261009.traj2 (notes sec. 4.7.7): n = the prefix of TD = 200
                     predicted periods, h2 = n // 2, the floor k = round((TD/2) / X_EDGE); dnu = nu2 / nu1 - 1;
                   - the files present: trace, psi6, virial, snapshot.
  snap_runs.csv    one row per snapshot file (sec. 4.8.6.4, the interior density): per snapshot phase (0 hold, 1 record) the mean over
                   its snapshots of each compartment's interior packing fraction: the disk centres inside the rectangle that leaves
                   2 sigma to every outer wall and to each face of the divider (its centre +- t/2, t = 0.05, at that snapshot),
                   counted only where that rectangle is at least 4 sigma wide, times pi/4 over its area.
Prints one line: the counts. Exit 0, or 2 with a STOP line (an input missing or inconsistent).
usage: python3 cluster/gen3_stage1_261012/reduce_s1_koa.py <.../m_<M>>
"""
import glob, gzip, math, os, re, struct, sys, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_testR_261011")); sys.path.insert(0, os.path.join(HS, "validation")); sys.path.insert(0, HS)
import reduce_testR_koa as RT                                   # nud_table, health_table, failed_list (Test R's, unchanged)
import resched_testT_followup_261007 as FU                       # fit_acf: the item 4a model
import tests_20260913 as T
from paper1_populate_cs_err_20261002 import TD, X_EDGE
warnings.filterwarnings("ignore")
WALL_T, STRIP, MINW = 0.05, 2.0, 4.0


def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def one(d, pat, r):
    f = glob.glob(os.path.join(d, pat.replace("<r>", str(r))) + "") + glob.glob(os.path.join(d, pat.replace("<r>", str(r))) + ".gz")
    return f[0] if len(f) == 1 else (None if not f else "AMBIGUOUS")


def acf_mean_fit(d):
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
    z = np.load(os.path.join(d, "acf_runs.npz")); nl = min(len(z[f"run{r}"]) for r in red["run"])
    dt = float(red["dt"].iloc[0]); nu_ref = float(red["nu"].mean())
    m = np.mean([np.asarray(z[f"run{r}"][:nl], dtype=float) for r in red["run"]], axis=0)
    try:
        p = FU.fit_acf(m, dt, nu_ref); ok = 1
    except Exception:
        p = [math.nan] * 6; ok = 0
    nd = p[4] / (2 * math.pi)
    pd.DataFrame([dict(n=len(red), nl=nl, dt=dt, nu_ref=nu_ref, A=p[0], tau_T=p[1], B=p[2], tau_r=p[3], omega=p[4], phi=p[5],
                       nu_d_meanfit=nd, Q_meanfit=math.pi * nd * p[3], fit_ok=ok)]).to_csv(os.path.join(d, "acf_mean_fit.csv"), index=False,
                                                                                          float_format="%.17g")
    return ok


def quarters(vb, Teq):
    v = pd.read_csv(vb); h = v[v["phase"] == "hold"]
    zc = (h["Z_0"] + h["Z_1"]) / 2.0; t = h["t_sigma"]
    q = lambda a, b: zc[(t > a * Teq + 1e-9) & (t <= b * Teq + 1e-9)]
    q3, q4 = q(0.50, 0.75), q(0.75, 1.0)
    return float(q4.mean()) if len(q4) else math.nan, float(q3.mean()) if len(q3) else math.nan, len(q3), len(q4)


def halves(tr):
    """paper1_window_aging_261009.traj2's split, with TD = 200 (the stage 1 record is 200 predicted periods)"""
    t, x, nup = T._load(tr); dt = (t[-1] - t[0]) / (len(t) - 1)
    n = T._prefix(t, nup, TD)
    if n is None: return math.nan, math.nan

    def argmax_nu(xx, td):
        P, df = T._spectrum(xx, dt); k = int(round(td / X_EDGE)); return (k + int(np.argmax(P[k:]))) * df
    h2 = n // 2
    return argmax_nu(x[:h2], TD / 2.0), argmax_nu(x[h2:2 * h2], TD / 2.0)


def psi6_means(ps):
    p = pd.read_csv(ps); rec = p[p["phase"] == "record"]; hold = p[p["phase"] == "hold"]
    if "psi6_0" not in p.columns or "psi6_1" not in p.columns: return [math.nan] * 5
    a, b = float(rec["psi6_0"].mean()), float(rec["psi6_1"].mean())
    return a, b, (a + b) / 2.0, float(hold["psi6_0"].iloc[-1]) if len(hold) else math.nan, float(hold["psi6_1"].iloc[-1]) if len(hold) else math.nan


def snapshot_reduce(sp):
    b = (gzip.open(sp, "rb") if sp.endswith(".gz") else open(sp, "rb")).read()
    if b[:8] != b"G3SNAP1\n": raise ValueError("not a G3SNAP1 file")
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, every = struct.unpack_from("<dddd", b, 16); off = 48
    acc = {0: [[], []], 1: [[], []]}; cnt = {0: 0, 1: 0}
    while off < len(b):
        t, ph, ndv = struct.unpack_from("<dii", b, off); off += 16
        cen = struct.unpack_from(f"<{ndv}d", b, off); off += 8 * ndv
        xy = np.frombuffer(b, dtype="<f8", count=2 * N, offset=off).reshape(N, 2); off += 16 * N
        cnt[ph] += 1
        if ndv != 1: continue
        c = cen[0]
        for k, (lo, hi) in enumerate(((STRIP, c - WALL_T / 2 - STRIP), (c + WALL_T / 2 + STRIP, W - STRIP))):
            if hi - lo < MINW: acc[ph][k].append(math.nan); continue
            y0, y1 = STRIP, H - STRIP
            if y1 - y0 <= 0: acc[ph][k].append(math.nan); continue
            inside = ((xy[:, 0] > lo) & (xy[:, 0] < hi) & (xy[:, 1] > y0) & (xy[:, 1] < y1)).sum()
            acc[ph][k].append(inside * math.pi / 4.0 / ((hi - lo) * (y1 - y0)))
    out = dict(N=N, box_w=W, box_h=H, radius=R, every=every, snaps_hold=cnt[0], snaps_record=cnt[1])
    for ph, lab in ((0, "hold"), (1, "record")):
        for k in (0, 1):
            v = np.array(acc[ph][k], dtype=float)
            out[f"eta_int_{lab}_{k}"] = float(np.nanmean(v)) if np.isfinite(v).any() else math.nan
            out[f"used_{lab}_{k}"] = int(np.isfinite(v).sum())
    return out


def main():
    if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]): sys.exit("usage: reduce_s1_koa.py <.../m_<M>>")
    d = sys.argv[1]
    for f in ("red_nu.csv", "acf_runs.npz", "run.log"):
        if not os.path.exists(os.path.join(d, f)): print(f"STOP: {d} has no {f}"); sys.exit(2)
    n, nf = RT.nud_table(d); nh = RT.health_table(d); nfl = RT.failed_list(d); ok = acf_mean_fit(d)
    red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run").reset_index(drop=True)
    log = open(os.path.join(d, "run.log"), errors="replace").read()
    rows, snaps = [], []
    for r in red["run"]:
        tr = one(d, "wall_x_positions_L0_*_run<r>.csv", r); ps = one(d, "psi6_t_L0_*_run<r>.csv", r)
        vb = one(d, "virial_blocks_L0_*_run<r>.csv", r); sn = one(d, "snapshots_L0_*_run<r>.bin", r)
        if "AMBIGUOUS" in (tr, ps, vb, sn): print(f"STOP: more than one file of run {r} in {d}"); sys.exit(2)
        sec = re.search(rf"^##RUN \S+ \S+ run {r} seed \d+ .*?(?=^##RUN |\Z)", log, re.M | re.S)
        hs = re.search(r"^Wall hold steps: (\d+)$", sec.group(0), re.M) if sec else None      # the binary's info block
        Teq = int(hs.group(1)) * 0.4 / 24.0 if hs else math.nan
        row = dict(run=int(r), seed=int(red.loc[red["run"] == r, "seed"].iloc[0]), T_eq=Teq, trace=int(tr is not None), psi6=int(ps is not None),
                   virial=int(vb is not None), snapshot=int(sn is not None))
        z4, z3, n3, n4 = quarters(vb, Teq) if vb and math.isfinite(Teq) else (math.nan, math.nan, 0, 0)
        row.update(Z=z4, Z_q3=z3, Z_q4=z4, blocks_q3=n3, blocks_q4=n4)
        a, b, m, h0, h1 = psi6_means(ps) if ps else [math.nan] * 5
        row.update(psi6_rec_0=a, psi6_rec_1=b, psi6_rec=m, psi6_hold_end_0=h0, psi6_hold_end_1=h1)
        n1, n2 = halves(tr) if tr else (math.nan, math.nan)
        row.update(nu1=n1, nu2=n2, dnu=n2 / n1 - 1 if n1 and math.isfinite(n1) else math.nan)
        rows.append(row)
        if sn:
            s = snapshot_reduce(sn); s.update(run=int(r), file=os.path.basename(sn)); snaps.append(s)
    pd.DataFrame(rows).to_csv(os.path.join(d, "s1_runs.csv"), index=False, float_format="%.17g")
    if snaps: pd.DataFrame(snaps).to_csv(os.path.join(d, "snap_runs.csv"), index=False, float_format="%.17g")
    miss = sum(1 for x in rows if not (x["trace"] and x["psi6"] and x["virial"]))
    print(f"{n} trajectories: nud_runs.csv ({nf} fit failures), runs_health.csv ({nh} ##RUN sections), failed_runs.csv ({nfl} failed runs), "
          f"acf_mean_fit.csv (fit_ok {ok}), s1_runs.csv ({len(rows)} rows; {miss} without trace, psi6 or virial file), snap_runs.csv ({len(snaps)} files)")


if __name__ == "__main__":
    main()
