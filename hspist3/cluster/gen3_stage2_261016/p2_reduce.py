#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 57; notes 261012 sec. 4.7.135): THE P2 REDUCTION, run by p2_merge.sbatch on KOA after the stage 1
merge and its reduction (merge_s1_koa.py: reduce_B.py and reduce_s1_koa.py per cell folder, unchanged). It reads the merged tree and adds,
per (cell, mass) folder merged/stage2/gen3/<cell>/m_<M>, two tables (nothing else is written or changed):
  p2_runs.csv    one row per trajectory (red_nu.csv's runs):
    - the run: seed, node, seconds, executed events, clean (runs_health.csv; run.log);
    - the estimators: nu_arg, nu_d and fit_ok as registered (red_nu.csv, nud_runs.csv); the item-4a fit of the run's own ACF (FU.fit_acf,
      unchanged; nu_ref the folder's argmax mean, as nud_runs.csv): nu_d, tau_r, Q = pi nu_d tau_r, A, tau_T, B, the rms of
      (normalised ACF - fit), and "stall" when A and B end within 1e-3 of the initial guess 0.3 and 0.7 (Report D56-1 (d));
    - the trace over the analysed prefix (200 predicted periods, T._prefix): mean, sd, drift (least squares slope, per 1000 sigma-time),
      the slow part (running mean over one predicted period) and the fast part, the share of the periodogram's power within 0.8-1.2 nu_d
      (the spectrum share), the Welch peak / nu_pred;
    - psi6 per compartment (psi6_t): the record means of |psi6| and the hold's last quarter; the orientation's circular mean over the
      record, mod 60 deg in [0, 60) (0: close-packed lines along x; 30: along y);
    - Z (s1_runs.csv, unchanged);
    - the strain: e_lattice from the run header's arrangement (a_s, a_r: the start's lattice before the jitter; p2_lattice.strain);
      per compartment e, theta0, the line direction and the ordered fraction (p2_lattice.config_strain) of the initial state
      (.done_runs/run<r>_chunk<j>/initial_state.txt), of the last snapshot of the hold and of the last snapshot (the record's end), with
      the divider's distance from L0 in the dump and in the last snapshot.
  p2_frames.csv  per trajectory and snapshot: t, phase, the divider - L0, and per compartment e, theta0, the mean local |psi6| of the bulk,
                 the ordered fraction, the disks used.
Prints one line per folder and a summary; writes merged/p2_reduce_summary.txt. A trajectory whose reduction raises gets its error text
in the column 'error' (the rest of the folder goes on); exit 0, or 2 when a folder lacks the stage 1 tables (nothing written there).
usage: python3 cluster/gen3_stage2_261016/p2_reduce.py <.../P2/merged> --tasks <tasks_p2_koa_261016.txt>
"""
import glob, gzip, io, math, os, re, sys, time, traceback, warnings
import numpy as np, pandas as pd
from scipy.signal import welch
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.dirname(HERE); HS = os.path.dirname(CL)
for p in (HERE, os.path.join(HS, "validation"), HS): sys.path.insert(0, p)
import p2_lattice as LT
import tests_20260913 as T
import resched_testT_followup_261007 as FU
from paper1_populate_cs_err_20261002 import TD
warnings.filterwarnings("ignore")
SIDE = re.compile(r"left=\[orient (\w+) rows (\w+) n_s (\d+) n_r (\d+) a_s (\S+) a_r (\S+) g (\S+) s (\S+) n_vac (\d+) jitter (\S+) px\]")


def circ60(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    if not len(a): return math.nan
    v = float(np.mod(np.degrees(np.angle(np.mean(np.exp(6j * np.radians(a))))) / 6.0, 60.0)); return 0.0 if v > 60.0 - 1e-6 else v


def one(d, pat):
    f = sorted(glob.glob(os.path.join(d, pat)))
    return f[0] if len(f) == 1 else None


def load_trace(p):
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
        hdr = fh.readline().rstrip("\n").split(","); first = fh.readline().rstrip("\n").split(",")
    nup = float(first[hdr.index("Predicted_Frequency")])
    d = pd.read_csv(p, usecols=["Time", "Displacement(σ)"])
    return d["Time"].to_numpy(float), d["Displacement(σ)"].to_numpy(float), nup


def run_sections(log):
    out = {}
    for sec in re.split(r"^##RUN ", log, flags=re.M)[1:]:
        h = re.match(r"\S+ \S+ run (\d+) seed (\d+) node (\S+) \((\d+) s\)", sec)
        if not h: continue
        ev = re.search(r"executed events (\d+)", sec); hd = re.search(r"^\[EDMD3-RUN\].*$", sec, re.M)
        out[int(h.group(1))] = dict(seed=int(h.group(2)), node=h.group(3), seconds=int(h.group(4)), events=int(ev.group(1)) if ev else -1,
                                    header=hd.group(0) if hd else "", m4c=(re.search(r"^\[EDMD3-M4C\].*$", sec, re.M) or [""])[0])
    return out


def strain_cols(prefix, res, cen, L0):
    o = {f"{prefix}_div_minus_L0": cen - L0}
    for k in (0, 1):
        r = res[k]
        o.update({f"{prefix}_e_{k}": r["e"], f"{prefix}_theta0_{k}": r["theta0"], f"{prefix}_line_{k}": r["line"],
                  f"{prefix}_psi6loc_{k}": r["psi6_local_bulk"], f"{prefix}_ordered_{k}": r["ordered_bulk"], f"{prefix}_used_{k}": r["n_used"]})
    return o


def reduce_run(d, r, red, nud, s1, hl, z, nu_ref, dt, sec, task, frames_out):
    row = dict(run=int(r)); tr = one(d, f"wall_x_positions_L0_*_run{r}.csv*")
    t, x, nup = load_trace(tr); n = T._prefix(t, nup, TD)
    if n is None: raise ValueError("the trace is shorter than 200 predicted periods")
    t, x = t[:n], x[:n]; ddt = (t[-1] - t[0]) / (len(t) - 1)
    rn = nud[nud["run"] == r].iloc[0]; rr = red[red["run"] == r].iloc[0]
    row.update(seed=int(rn["seed"]), nu_pred=nup, nu_arg=float(rr["nu"]), nu_d=float(rn["nu_d"]), fit_ok=int(rn["fit_ok"]))
    try:
        p = FU.fit_acf(np.asarray(z[f"run{r}"], float), dt, nu_ref); c = np.asarray(z[f"run{r}"], float); c = c / c[0]
        g = FU._g(np.arange(len(c)) * dt, *p); nd = p[4] / (2 * math.pi)
        row.update(fit_nu_d=nd, fit_tau_r=p[3], fit_Q=math.pi * nd * p[3], fit_A=p[0], fit_tau_T=p[1], fit_B=p[2],
                   fit_rms=float(np.sqrt(np.mean((c - g) ** 2))), stall=int(abs(p[0] - 0.3) < 1e-3 and abs(p[2] - 0.7) < 1e-3), fit_raised=0)
    except Exception:
        row.update(fit_nu_d=math.nan, fit_tau_r=math.nan, fit_Q=math.nan, fit_A=math.nan, fit_tau_T=math.nan, fit_B=math.nan, fit_rms=math.nan,
                   stall=0, fit_raised=1)
    w = max(1, int(1.0 / (nup * ddt))); sl = np.convolve(x, np.ones(w) / w, mode="same"); fast = x - sl
    P, df = T._spectrum(x, ddt); f = np.arange(len(P)) * df; ref = row["nu_d"] if math.isfinite(row["nu_d"]) else row["nu_arg"]
    band = (f > 0.8 * ref) & (f < 1.2 * ref); fw, Pw = welch(x - x.mean(), fs=1 / ddt, nperseg=max(16, n // 8), window="hann", detrend="constant")
    row.update(x_mean=float(x.mean()), x_sd=float(x.std()), drift_per_1000=1000 * float(np.polyfit(t, x, 1)[0]),
               slow_sd=float(sl[w:-w].std()) if n > 2 * w else math.nan, fast_sd=float(fast[w:-w].std()) if n > 2 * w else math.nan,
               spectrum_share=float(P[band].sum() / P[1:].sum()), welch_peak_over_nup=float(fw[1 + np.argmax(Pw[1:])] / nup))
    ps = pd.read_csv(one(d, f"psi6_t_L0_*_run{r}.csv*")); rec = ps[ps["phase"] == "record"]; hold = ps[ps["phase"] == "hold"]
    teq = float(s1[s1["run"] == r]["T_eq"].iloc[0]) if "T_eq" in s1 else 5000.0; h4 = hold[hold["t_sigma"] > 0.75 * teq]
    for k in (0, 1):
        row.update({f"psi6_rec_{k}": float(rec[f"psi6_{k}"].mean()), f"psi6_holdq4_{k}": float(h4[f"psi6_{k}"].mean()),
                    f"orient_rec_{k}": circ60(rec[f"psi6_orient_deg_{k}"]), f"orient_holdq4_{k}": circ60(h4[f"psi6_orient_deg_{k}"])})
    srow = s1[s1["run"] == r].iloc[0]; row.update(Z=float(srow["Z"]), Z_q3=float(srow["Z_q3"]))
    hrow = hl[hl["run"] == r]; row.update(clean=int(hrow["gen3_clean"].iloc[0]) if len(hrow) else -1)
    row.update(node=sec.get("node", ""), seconds=sec.get("seconds", -1), events=sec.get("events", -1),
               m4c=sec.get("m4c", "").replace("[EDMD3-M4C] lattice orientation requested: ", ""))
    m = SIDE.search(sec.get("header", ""))
    if m:
        a = dict(ns=int(m.group(3)), nr=int(m.group(4)), a_s=float(m.group(5)), a_r=float(m.group(6)))
        row.update(arrangement=f"{m.group(1)[:4]} {m.group(2)[:3]} {a['ns']}x{a['nr']}", n_vac=int(m.group(9)), gap_sigma=float(m.group(8)) / LT.PX,
                   e_lattice=LT.strain(a))
    L0 = float(task[5]); dumps = glob.glob(os.path.join(d, ".done_runs", f"run{r}_chunk*", "initial_state.txt")) + \
        glob.glob(os.path.join(d, ".done_runs", f"run{r}", "initial_state.txt"))
    if len(dumps) == 1:
        D = LT.read_dump(dumps[0]); res = LT.config_strain(D["xy"], D["cen"], D["W"], D["H"], th=D["th"])
        row.update(strain_cols("start", res, D["cen"], L0)); row.update(start_NL=int((D["xy"][:, 0] < D["cen"]).sum()), start_N=D["N"])
    sn = one(d, f"snapshots_L0_*_run{r}.bin*")
    if sn:
        S = LT.read_snapshots(sn); last_hold = None
        for i, fr in enumerate(S["frames"]):
            res = LT.config_strain(fr["xy"], fr["cen"], S["W"], S["H"])
            fo = dict(run=int(r), frame=i, t=fr["t"], phase=fr["phase"], div_minus_L0=fr["cen"] - L0)
            for k in (0, 1):
                fo.update({f"e_{k}": res[k]["e"], f"theta0_{k}": res[k]["theta0"], f"psi6loc_{k}": res[k]["psi6_local_bulk"],
                           f"ordered_{k}": res[k]["ordered_bulk"], f"used_{k}": res[k]["n_used"]})
            frames_out.append(fo)
            if fr["phase"] == 0: last_hold = (res, fr["cen"])
        row.update(strain_cols("hold", last_hold[0], last_hold[1], L0) if last_hold else {})
        if S["frames"]: row.update(strain_cols("end", res, S["frames"][-1]["cen"], L0))
        row.update(snapshots=len(S["frames"]))
    return row


def main():
    if len(sys.argv) < 2 or "--tasks" not in sys.argv: sys.exit(__doc__)
    M = os.path.abspath(sys.argv[1]); TL = [l.split() for l in open(sys.argv[sys.argv.index("--tasks") + 1]) if l.strip()]
    tasks = {(l[1], int(l[3])): l for l in TL}
    folders = sorted(glob.glob(os.path.join(M, "stage2", "gen3", "*", "m_*")))
    if not folders: print(f"STOP: no folder merged/stage2/gen3/<cell>/m_<M> under {M}"); sys.exit(2)
    t0 = time.time(); summ, nerr, nrun, bad = [], 0, 0, 0
    for d in folders:
        rel = os.path.relpath(d, M); need = ("red_nu.csv", "nud_runs.csv", "acf_runs.npz", "s1_runs.csv", "runs_health.csv", "run.log")
        miss = [f for f in need if not os.path.exists(os.path.join(d, f))]
        if miss: print(f"{rel}: MISSING {miss} -- not reduced"); bad += 1; continue
        red = pd.read_csv(os.path.join(d, "red_nu.csv")).sort_values("run"); nud = pd.read_csv(os.path.join(d, "nud_runs.csv"))
        s1 = pd.read_csv(os.path.join(d, "s1_runs.csv")); hl = pd.read_csv(os.path.join(d, "runs_health.csv")); z = np.load(os.path.join(d, "acf_runs.npz"))
        nl = min(len(z[f"run{r}"]) for r in red["run"]); zz = {f"run{r}": z[f"run{r}"][:nl] for r in red["run"]}
        nu_ref = float(red["nu"].mean()); dt = float(red["dt"].iloc[0]); secs = run_sections(open(os.path.join(d, "run.log"), errors="replace").read())
        rows, frames = [], []
        for r in red["run"]:
            task = tasks.get((rel, int(r)))
            try:
                if task is None: raise ValueError(f"no task line for {rel} s={r}")
                row = reduce_run(d, int(r), red, nud, s1, hl, zz, nu_ref, dt, secs.get(int(r), {}), task, frames); row["error"] = ""
            except Exception as e:
                row = dict(run=int(r), error=f"{type(e).__name__}: {e}".replace("\n", " ")[:300]); nerr += 1
            cid = rel.split(os.sep)[2]; f = cid.split("_")
            row.update(cell=cid, eta=float(f[2][3:]), G=f[3], pin=f[4][3:], M=int(rel.split("_")[-1]),
                       H=float(task[6]) if task else math.nan, L0=float(task[5]) if task else math.nan)
            rows.append(row); nrun += 1
        R = pd.DataFrame(rows); lead = ["cell", "eta", "G", "pin", "H", "L0", "M", "run"]
        R = R[lead + [c for c in R.columns if c not in lead]]; R.to_csv(os.path.join(d, "p2_runs.csv"), index=False, float_format="%.10g")
        if frames: pd.DataFrame(frames).to_csv(os.path.join(d, "p2_frames.csv"), index=False, float_format="%.8g")
        g = lambda c: float(R[c].mean()) if c in R else math.nan
        line = (f"{rel}: {len(R)} runs, errors {int((R['error'] != '').sum())}; e start {100 * g('start_e_0'):+.2f} / {100 * g('start_e_1'):+.2f} %, "
                f"end {100 * g('end_e_0'):+.2f} / {100 * g('end_e_1'):+.2f} %; Q (own fits, median) {R['fit_Q'].median() if 'fit_Q' in R else math.nan:.1f}; "
                f"stalls {int(R['stall'].sum()) if 'stall' in R else 0}; psi6 rec {g('psi6_rec_0'):.3f} / {g('psi6_rec_1'):.3f}")
        print(line, flush=True); summ.append(line)
    tail = (f"p2_reduce: {len(folders)} folders ({bad} without the stage 1 tables), {nrun} trajectories, {nerr} with an error; "
            f"{time.time() - t0:.0f} s")
    open(os.path.join(M, "p2_reduce_summary.txt"), "w").write("\n".join(summ + [tail]) + "\n"); print(tail)
    sys.exit(2 if bad else 0)


if __name__ == "__main__":
    main()
