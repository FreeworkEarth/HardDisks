#!/usr/bin/env python3
"""##CHRIS 2026-09-12: one pass over every A2 (fixed-eta N ladder) trace, caching BOTH
frequency estimators so the cut-sensitivity study never refits.

Per trace:
  nu_damped, sigma_nu   time-domain damped-cosine fit (fit_nu_damped.fit_trace)
  nu_binned, df_binned  peak of the power spectrum of the same post-transient,
                        mean-removed displacement. Only the DC bin is excluded, so
                        this is the assumption-free "no cut" estimator; df = 1/T is
                        its quantisation and the reason it cannot resolve better
                        than a few percent on these traces.
  nu_predicted          the driver's own estimate, carried in the trace header
  health counters, T_i  from run.log (strict contract + HD_KE_TRACE audit)
Read-only on every trace.
"""
import csv, glob, math, os, re, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit_nu_damped import fit_trace, DROP

# ##CHRIS: comma-separated roots so the original A2 tree and its top-up tree merge
# into one cache. The top-up is written to a separate tree because the runner numbers
# traces from run0 in each leaf and would otherwise overwrite the originals.
CAMP = sys.argv[1]
OUT = sys.argv[2]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 8
ROOTS = [c for c in CAMP.split(",") if c]


def binned_nu(path):
    t, x = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                t.append(float(row["Time"])); x.append(float(row["Displacement(σ)"]))
            except (KeyError, ValueError):
                pass
    if len(t) < 50:
        return float("nan"), float("nan")
    t = np.asarray(t); x = np.asarray(x)
    i0 = int(DROP * len(t)); t = t[i0:] - t[i0]; x = x[i0:] - x[i0:].mean()
    dt = float(np.median(np.diff(t)))
    if not (dt > 0):
        return float("nan"), float("nan")
    P = np.abs(np.fft.rfft(x)) ** 2
    f = np.fft.rfftfreq(len(x), d=dt)
    P[0] = 0.0                      # DC only; no band, no window, no other cut
    k = int(np.argmax(P))
    return float(f[k]), float(f[1] - f[0]) if len(f) > 1 else float("nan")


def scan_cell(cell):
    log_p = os.path.join(cell, "run.log")
    health, ti_of, pend = {}, {}, None
    if os.path.exists(log_p):
        for ln in open(log_p, errors="replace"):
            k = re.search(r"2b after per-segment equalize\s+N=(\d+)\s+KE_tot=\S+\s+"
                          r"KE_left=(\S+)\s+KE_right=(\S+)", ln)
            if k:
                ns = int(k.group(1)) / 2.0
                pend = (float(k.group(2)) / ns, float(k.group(3)) / ns)
            m = re.search(r"Running: L0 = [\d.]+, M = (\d+)\*m, run = (\d+), seed = (\d+)", ln)
            if m and pend:
                ti_of[(int(m.group(1)), int(m.group(2)))] = pend; pend = None
            h = re.search(r"EDMD-HEALTH\] L0=[\d.]+ M=(\d+) run=(\d+) seed=(\d+): "
                          r"forced_advance=(\d+) wall_clamp_repairs=(\d+) "
                          r"overlap_repairs=(\d+) wall_overdue=(\d+)", ln)
            if h:
                health[(int(h.group(1)), int(h.group(2)))] = tuple(int(h.group(i)) for i in (4, 5, 6, 7))
    mm = re.search(r"eta_([\dp]+)/N(\d+)/m_(\d+)$", cell)
    if not mm:
        return []
    eta = float(mm.group(1).replace("p", ".")); N = int(mm.group(2))
    rows = []
    for f in sorted(glob.glob(os.path.join(cell, "wall_x_positions_*.csv"))):
        g = re.search(r"wallmassfactor_(\d+)_run(\d+)\.csv$", f)
        if not g:
            continue
        M, run = int(g.group(1)), int(g.group(2))
        with open(f) as fh:
            row0 = next(csv.DictReader(fh), None)
        if row0 is None:
            continue
        res = fit_trace(f)
        nb, df = binned_nu(f)
        fa, cr, orp, wo = health.get((M, run), (0, 0, 0, 0))
        ti = ti_of.get((M, run))
        rows.append(dict(
            eta=f"{eta:.6g}", N=N, M=M, run=run, L0=row0.get("L0", ""),
            nu_damped=f"{res['nu']:.10g}" if res else "",
            sigma_nu=f"{res['sigma_nu']:.6g}" if res else "",
            gamma=f"{res['gamma']:.6g}" if res else "",
            fit_rms=f"{res['rms']:.6g}" if res else "",
            nu_binned=f"{nb:.10g}" if nb == nb else "",
            df_binned=f"{df:.6g}" if df == df else "",
            nu_predicted=row0.get("Predicted_Frequency", ""),
            forced_advance=fa, clamp_repair=cr, overlap_repair=orp, wall_overdue=wo,
            T_i_left=f"{ti[0]:.9g}" if ti else "", T_i_right=f"{ti[1]:.9g}" if ti else "",
            path=os.path.relpath(f), tree=os.path.basename(os.path.dirname(
                os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(f))))))))
    return rows


if __name__ == "__main__":
    cells = []
    for root in ROOTS:
        got = sorted(d for d in glob.glob(f"{root}/eta_*/N*/m_*") if os.path.isdir(d))
        print(f"  {root}: {len(got)} cells")
        cells += got
    print(f"cells {len(cells)} from {len(ROOTS)} tree(s)", flush=True)
    with Pool(NW) as p:
        out = p.map(scan_cell, cells, chunksize=1)
    rows = [r for sub in out for r in sub]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    nd = sum(1 for r in rows if r["nu_damped"])
    nb = sum(1 for r in rows if r["nu_binned"])
    hl = sum(1 for r in rows if any((r["forced_advance"], r["clamp_repair"],
                                     r["overlap_repair"], r["wall_overdue"])))
    print(f"traces {len(rows)}  damped-fit ok {nd}  binned ok {nb}  any health event {hl}")
    print(f"-> {OUT}")
