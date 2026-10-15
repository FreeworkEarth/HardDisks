#!/usr/bin/env python3
"""##CHRIS 2026-10-15 (plan-author decision 56, task 2 c; notes 261012 sec. 4.7.133): DOES THE DAMPED-COSINE ESTIMATOR SURVIVE IN THE
SOLID? The five Mac runs of d56_pilot.py: (N 900, eta 0.69, M 4500) s = 0, 1, 2 (the registered stage 1 lines, rerun on the Mac: same
protocol and seeds, other trajectories than KOA's) and the Stage 2 pilot cell (N 900, eta 0.76, M 4500) s = 0, 1.
  1. The registered reduction, unchanged: cluster/confinement_20261013/reduce_B.py per cell (red_nu.csv: the argmax estimator; acf_runs.npz:
     the position ACF to 20 predicted periods), then cluster/gen3_stage1_261012/reduce_s1_koa.py per m_<M> (nud_runs.csv: nu_d per
     trajectory by the item-4a fit, nu_ref = the folder's argmax mean; acf_mean_fit.csv; s1_runs.csv: Z, psi6, the halves; snap_runs.csv).
     Run once: a folder that already holds red_nu.csv is not reduced again.
  2. Per trajectory: the item-4a fit of its own ACF (FU.fit_acf, all six parameters): nu_d, tau_r, Q = pi nu_d tau_r, the deviation of
     the fitted curve from the normalised ACF (rms and largest, over the fitted lags); the periodogram (T._spectrum, the argmax estimator's
     own) and a Welch estimate (8 segments, Hann); the trace's mean, drift (least squares slope), its slow part (the running mean over
     one predicted period) against its fast part; the amplitude check: the trace's sd against the sd that the cell's own E1 speed
     implies for a divider between two columns of that stiffness, sd_E1 = sqrt(kT L_eff / (2 H rho m c^2)) (a locked divider would sit
     far below it; information, not an estimator).
  3. Per cell, E1 exactly as s1_verdict_261012.per_mass on the trajectories present (nu_d-bar, tau_r of the seed-averaged ACF,
     nu_0 = sqrt(nu_d-bar^2 + (1 / (2 pi tau_r))^2), Q = pi nu_d-bar tau_r, c = nu_0 / x, x = K(alpha) / (2 pi (L0 - 1.025))); the
     jackknife needs 3 trajectories, so the pilot's 2 give no SE. For (900, 0.69) the same from KOA's fetched tables: the 3 seeds and
     all 25 (information).
  4. psi6 per compartment over time (psi6_t, hold and record) and the crystal per snapshot (local psi6 per disk from the snapshots,
     neighbours within 1.35 lattice spacings on the disk's own side of the divider; colour = the local orientation, mod 60 deg).
Writes into experiments_gen3_d56_261015/c_analysis/ (a new folder): c_runs.csv, c_cells.csv, c_summary.txt, and the figures
fig_c1_N0900_eta0.69_traces_fits_spectra_psi6.png, fig_c2_N0900_eta0.76_traces_fits_spectra_psi6.png, fig_c3_snapshots_local_psi6.png.
usage (from hspist3/): python3 experiments_gen3_d56_261015/scripts/d56_c_analysis.py [--out <a new folder instead>]
"""
import glob, gzip, math, os, re, struct, subprocess, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import welch
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d56_261015")
for p in (os.path.join(HS, "cluster", "gen3_stage1_261012"), os.path.join(HS, "cluster", "gen3_testR_261011"), os.path.join(HS, "validation"), HS):
    sys.path.insert(0, p)
import tests_20260913 as T
import resched_testT_followup_261007 as FU
import plot_speed_of_sound_edmd as SOS
from paper1_populate_cs_err_20261002 import TD
warnings.filterwarnings("ignore")
OUT = os.path.join(E, "c_analysis")
KOA69 = os.path.join(HS, "experiments_gen3_stage1_koa_261012", "B900", "merged", "stage1", "gen3", "B_N0900_eta0.690000", "m_4500")
CELLS = [dict(key="69", eta=0.69, N=900, H=30.0, L0=17.073873117335832, alpha=5, M=4500, Teq=5000.0,
              cell=os.path.join(E, "data", "stage1", "gen3", "B_N0900_eta0.690000"),
              title="N 900, eta 0.69, M 4500 (alpha 5): the registered stage 1 lines s = 0, 1, 2 rerun on the Mac"),
         dict(key="76", eta=0.76, N=900, H=30.0, L0=15.501279540739111, alpha=5, M=4500, Teq=5000.0,
              cell=os.path.join(E, "data", "pilot", "P2_N0900_eta0.760000"),
              title="N 900, eta 0.76, M 4500 (alpha 5): the Stage 2 pilot cell, s = 0, 1")]
C_L, C_R = "#2a6fbb", "#d9822b"            # the two compartments (left, right)
C_TR, C_FIT = "#1f3b5a", "#c0392b"


def reduce_once(c):
    m = os.path.join(c["cell"], f"m_{c['M']}")
    if os.path.exists(os.path.join(m, "red_nu.csv")): return "already reduced"
    a = subprocess.run([sys.executable, os.path.join(HS, "cluster", "confinement_20261013", "reduce_B.py"), c["cell"]], capture_output=True, text=True)
    b = subprocess.run([sys.executable, os.path.join(HS, "cluster", "gen3_stage1_261012", "reduce_s1_koa.py"), m], capture_output=True, text=True)
    return f"reduce_B exit {a.returncode}: {a.stdout.strip()} {a.stderr.strip()[-200:]} | reduce_s1_koa exit {b.returncode}: {b.stdout.strip()} {b.stderr.strip()[-200:]}"


def fit6(a, dt, nu_ref):
    p = FU.fit_acf(np.asarray(a, float), dt, nu_ref); c = np.asarray(a, float) / a[0]; lag = np.arange(len(c)) * dt
    g = FU._g(lag, *p); return p, c, lag, g


def e1_cell(nd, acfs, dt, nu_ref, L0, alpha):
    """s1_verdict_261012.per_mass on these trajectories (the seed-averaged fit for tau_r; the jackknife from 3 trajectories)"""
    nl = min(len(a) for a in acfs); m = np.mean([np.asarray(a[:nl], float) for a in acfs], axis=0)
    p = FU.fit_acf(m, dt, nu_ref); ndb = float(np.mean(nd)); tau_r = p[3]
    nu0 = math.sqrt(ndb ** 2 + (1 / (2 * math.pi * tau_r)) ** 2); x = T.k_root(alpha) / (2 * math.pi * (L0 - 1.025))
    jk = []
    if len(nd) >= 3:
        for i in range(len(nd)):
            keep = [k for k in range(len(nd)) if k != i]; mi = np.mean([np.asarray(acfs[k][:nl], float) for k in keep], axis=0)
            fi = FU.fit_acf(mi, dt, nu_ref); jk.append(math.sqrt(np.mean([nd[k] for k in keep]) ** 2 + (1 / (2 * math.pi * fi[3])) ** 2))
    n = len(jk); se = math.sqrt((n - 1) / n * sum((v - np.mean(jk)) ** 2 for v in jk)) if n >= 3 else math.nan
    return dict(n=len(nd), nud=ndb, tau_r=tau_r, Q=math.pi * ndb * tau_r, nu0=nu0, nu0_se=se, x=x, c=nu0 / x, sig_c=se / x, A=p[0], tau_T=p[1], B=p[2])


def koa_cell(seeds=None):
    nud = pd.read_csv(os.path.join(KOA69, "nud_runs.csv")); red = pd.read_csv(os.path.join(KOA69, "red_nu.csv")); z = np.load(os.path.join(KOA69, "acf_runs.npz"))
    if seeds is not None: nud = nud[nud["seed"].isin(seeds)]
    nud = nud[nud["fit_ok"] == 1].sort_values("run")
    return e1_cell(nud["nu_d"].to_numpy(float), [z[f"run{r}"] for r in nud["run"]], float(red["dt"].iloc[0]), float(red["nu"].mean()), 17.073873117335832, 5)


def cs_refs(eta):
    a = np.array([eta]); out = {}
    if eta <= 0.69 + 1e-12:
        Z, dZ = SOS.Z_kolafa_rottner_2006(a)[0], SOS.dZ_kolafa_rottner_2006(a)[0]; out["ref"] = "KR"
    else:
        Z, dZ = SOS.Z_liu_solid(a)[0], SOS.dZ_liu_solid(a)[0]; out["ref"] = "Liu solid branch (bulk only: a lower bound in the solid)"
    out.update(Z=Z, cs=float(SOS.cs_adiabatic_2d_monatomic(np.array([Z]), np.array([dZ]), a, kbt=1.0, m=1.0)[0]), cT=math.sqrt(Z + eta * dZ))
    Zh, dZh = SOS.Z_henderson_eos(a)[0], SOS.dZ_henderson_eos(a)[0]
    out["cs_henderson"] = float(SOS.cs_adiabatic_2d_monatomic(np.array([Zh]), np.array([dZh]), a, kbt=1.0, m=1.0)[0])
    return out


def circ60(a):
    """the circular mean of a psi6 orientation series [deg], mod 60 in [0, 60) (0: close-packed lines along x; 30: along y)"""
    return float(np.mod(np.degrees(np.angle(np.mean(np.exp(6j * np.radians(np.asarray(a, float)))))) / 6.0, 60.0))


def running_mean(x, w):
    w = max(1, int(w)); k = np.ones(w) / w; return np.convolve(x, k, mode="same"), w


def snap_read(p):
    b = (gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")).read()
    if b[:8] != b"G3SNAP1\n": raise ValueError("not a G3SNAP1 file")
    N, nd = struct.unpack_from("<ii", b, 8); W, H, R, every = struct.unpack_from("<dddd", b, 16); off = 48; out = []
    while off < len(b):
        t, ph, ndv = struct.unpack_from("<dii", b, off); off += 16
        cen = struct.unpack_from(f"<{ndv}d", b, off); off += 8 * ndv
        xy = np.frombuffer(b, dtype="<f8", count=2 * N, offset=off).reshape(N, 2).copy(); off += 16 * N
        out.append((t, ph, cen[0] if ndv else math.nan, xy))
    return dict(N=N, W=W, H=H, R=R), out


def local_psi6(xy, c, rc):
    side = xy[:, 0] > c; d = xy[None, :, :] - xy[:, None, :]; r = np.hypot(d[..., 0], d[..., 1])
    nb = (r < rc) & (r > 0) & (side[None, :] == side[:, None]); th = np.arctan2(d[..., 1], d[..., 0])
    z = np.where(nb, np.exp(6j * th), 0).sum(axis=1); k = nb.sum(axis=1)
    psi = np.where(k > 0, z / np.maximum(k, 1), 0); return np.abs(psi), np.degrees(np.angle(psi)) / 6.0, side


def per_run(c, red, nud, s1, hl, z, nu_ref, c_e1):
    m = os.path.join(c["cell"], f"m_{c['M']}"); rows, curves = [], []
    refs = cs_refs(c["eta"]); rho = 4 * c["eta"] / math.pi; Leff = c["L0"] - 1.025
    for r in red.sort_values("run")["run"]:
        tr = glob.glob(os.path.join(m, f"wall_x_positions_L0_*_run{r}.csv"))[0]
        t, x, nup = T._load(tr); dt = (t[-1] - t[0]) / (len(t) - 1); n = T._prefix(t, nup, TD); t, x = t[:n], x[:n]
        p, cn, lag, g = fit6(z[f"run{r}"], dt, nu_ref); nd = p[4] / (2 * math.pi)
        P, df = T._spectrum(x, dt); fw, Pw = welch(x - x.mean(), fs=1 / dt, nperseg=n // 8, window="hann", detrend="constant")
        sl, w = running_mean(x, 1 / (nup * dt)); fast = x - sl; drift = np.polyfit(t, x, 1)[0]
        band = (np.arange(len(P)) * df > 0.8 * nd) & (np.arange(len(P)) * df < 1.2 * nd)
        sd_e1 = math.sqrt(Leff / (2 * c["H"] * rho * c_e1 ** 2))
        ps = pd.read_csv(glob.glob(os.path.join(m, f"psi6_t_L0_*_run{r}.csv.gz"))[0])
        rec = ps[ps["phase"] == "record"]; hold = ps[ps["phase"] == "hold"]; h4 = hold[hold["t_sigma"] > 0.75 * c["Teq"]]
        srow = s1[s1["run"] == r].iloc[0]; hrow = hl[hl["run"] == r].iloc[0]; nrow = nud[nud["run"] == r].iloc[0]
        ev = re.search(rf"^##RUN \S+ \S+ run {r} seed \d+ .*?executed events (\d+)", open(os.path.join(m, "run.log")).read(), re.M | re.S)
        rows.append(dict(cell=c["key"], eta=c["eta"], run=int(r), seed=int(nrow["seed"]), mac_seconds=int(hrow["seconds"]), events=int(ev.group(1)),
                         gen3_clean=int(hrow["gen3_clean"]), nu_pred=nup, nu_arg=float(nrow["nu_arg"]), nu_d=float(nrow["nu_d"]), fit_ok=int(nrow["fit_ok"]),
                         nu_d_ownfit=nd, tau_r=p[3], Q_run=math.pi * nd * p[3], A=p[0], tau_T=p[1], B=p[2], phi=p[5],
                         fit_rms=float(np.sqrt(np.mean((cn - g) ** 2))), fit_maxdev=float(np.max(np.abs(cn - g))),
                         x_mean=float(x.mean()), x_sd=float(x.std()), drift_per_1000=1000 * drift, slow_sd=float(sl[w:-w].std()), fast_sd=float(fast[w:-w].std()),
                         band_power_frac=float(P[band].sum() / P[1:].sum()), welch_peak=float(fw[1 + np.argmax(Pw[1:])]),
                         x_sd_E1=sd_e1, x_sd_over_E1=float(x.std()) / sd_e1,
                         Z=float(srow["Z"]), psi6_rec_0=float(srow["psi6_rec_0"]), psi6_rec_1=float(srow["psi6_rec_1"]),
                         psi6_hold_q4_0=float(h4["psi6_0"].mean()), psi6_hold_q4_1=float(h4["psi6_1"].mean()),
                         orient_rec_0=circ60(rec["psi6_orient_deg_0"]), orient_rec_1=circ60(rec["psi6_orient_deg_1"]),
                         nu1=float(srow["nu1"]), nu2=float(srow["nu2"]), dnu=float(srow["dnu"])))
        curves.append(dict(r=int(r), t=t, x=x, sl=sl, w=w, nup=nup, lag=lag, cn=cn, g=g, P=P, df=df, fw=fw, Pw=Pw, ps=ps, nd=nd))
    return rows, curves


def orient_unwrap(a):
    """the representation without a wrap at the series' own mean orientation: [0, 60) near 30 deg, (-30, 30] near 0"""
    a = np.asarray(a, float); ref = np.degrees(np.angle(np.mean(np.exp(6j * np.radians(a))))) / 6.0
    return np.mod(a, 60.0) if abs(ref) > 15 else (np.mod(a + 30.0, 60.0) - 30.0)


def draw_cell(c, rows, curves, cell, koa, path):
    nr = len(curves); fig = plt.figure(figsize=(25, 4.3 * nr + 1.2))
    gs = fig.add_gridspec(nr, 5, width_ratios=[1.35, 1.0, 1.15, 1.15, 1.2], wspace=0.28, hspace=0.42)
    for i, (cu, rw) in enumerate(zip(curves, rows)):
        a0 = fig.add_subplot(gs[i, 0]); a0.plot(cu["t"], cu["x"], color=C_TR, lw=0.35, alpha=0.55, label="divider displacement")
        w = cu["w"]; a0.plot(cu["t"][w:-w], cu["sl"][w:-w], color=C_FIT, lw=1.6, label="running mean over one predicted period")
        a0.axhline(0, color="0.5", lw=0.6, ls=":")
        a0.set_xlabel("record time  [sigma sqrt(m/kT)]", fontsize=8); a0.set_ylabel("displacement  [sigma]", fontsize=8)
        a0.set_title(f"s = {rw['run']}  (seed {rw['seed']}): the whole record, {TD} predicted periods\nmean {rw['x_mean']:+.3f}, sd {rw['x_sd']:.3f}, "
                     f"drift {rw['drift_per_1000']:+.3f} per 1000; slow part sd {rw['slow_sd']:.3f}, fast {rw['fast_sd']:.3f}; sd / sd_E1 = {rw['x_sd_over_E1']:.2f}", fontsize=8)
        if i == 0: a0.legend(fontsize=7, loc="upper right", framealpha=0.9)
        a1 = fig.add_subplot(gs[i, 1]); T10 = 10 / cu["nup"]; k = cu["t"] <= cu["t"][0] + T10
        a1.plot(cu["t"][k], cu["x"][k], color=C_TR, lw=0.9); a1.axhline(0, color="0.5", lw=0.6, ls=":")
        a1.set_title("the first 10 predicted periods", fontsize=8); a1.set_xlabel("record time", fontsize=8)
        a2 = fig.add_subplot(gs[i, 2]); a2.plot(cu["lag"] * cu["nup"], cu["cn"], color=C_TR, lw=1.0, label="ACF of the displacement (normalised)")
        a2.plot(cu["lag"] * cu["nup"], cu["g"], color=C_FIT, lw=1.2, ls="--", label="item-4a fit: A e^(-t/tau_T) + B e^(-t/tau_r) cos(wt + phi)")
        a2.axhline(0, color="0.5", lw=0.6, ls=":"); a2.set_xlabel("lag  [predicted periods]", fontsize=8)
        a2.set_title(f"nu_d {rw['nu_d_ownfit']:.5f}, tau_r {rw['tau_r']:.1f}, Q = {rw['Q_run']:.2f};  A {rw['A']:.2f}, tau_T {rw['tau_T']:.0f}, B {rw['B']:.2f}\n"
                     f"fit minus ACF: rms {rw['fit_rms']:.3f}, largest {rw['fit_maxdev']:.3f}", fontsize=8)
        if i == 0: a2.legend(fontsize=7, loc="upper right", framealpha=0.9)
        a3 = fig.add_subplot(gs[i, 3]); f = np.arange(len(cu["P"])) * cu["df"]
        a3.semilogy(f[1:] / cu["nup"], cu["P"][1:] * cu["df"], color="0.75", lw=0.5, label="periodogram (the argmax estimator's)")
        a3.semilogy(cu["fw"][1:] / cu["nup"], cu["Pw"][1:] / cu["Pw"][1:].max() * (cu["P"][1:] * cu["df"]).max(),
                    color=C_TR, lw=1.3, label="Welch, 8 segments (scaled to the periodogram's peak)")
        for v, ls, lab in ((1.0, ":", "nu_pred (binary, Henderson)"), (rw["nu_d_ownfit"] / cu["nup"], "--", "nu_d (fit)"), (rw["nu_arg"] / cu["nup"], "-.", "nu_arg")):
            a3.axvline(v, color=C_FIT if ls == "--" else "0.35", ls=ls, lw=1.0, label=lab)
        a3.set_xlim(0, 4); a3.set_xlabel("frequency / nu_pred", fontsize=8); pk = (cu["P"][1:] * cu["df"]).max(); a3.set_ylim(pk * 1e-6, pk * 3)
        a3.set_title(f"power within 0.8-1.2 nu_d: {100 * rw['band_power_frac']:.0f} %;  Welch peak {rw['welch_peak'] / cu['nup']:.2f} nu_pred\n"
                     f"halves: nu1 {rw['nu1']:.5f}, nu2 {rw['nu2']:.5f} (dnu {100 * rw['dnu']:+.1f} %)", fontsize=8)
        if i == 0: a3.legend(fontsize=6.5, loc="lower left", framealpha=0.9)
        sub = gs[i, 4].subgridspec(2, 1, hspace=0.08); b0 = fig.add_subplot(sub[0]); b1 = fig.add_subplot(sub[1], sharex=b0)
        ps = cu["ps"]
        for k_, col, lab in ((0, C_L, "left"), (1, C_R, "right")):
            b0.plot(ps["t_sigma"], ps[f"psi6_{k_}"], color=col, lw=0.8, label=f"{lab} compartment")
            b1.plot(ps["t_sigma"], orient_unwrap(ps[f"psi6_orient_deg_{k_}"]), ".", color=col, ms=1.2, alpha=0.5)
        for b in (b0, b1): b.axvline(c["Teq"], color="0.3", lw=0.8, ls="--")
        b0.set_ylim(0, 1); b0.set_ylabel("|psi6|", fontsize=8); b1.set_ylabel("orientation [deg]", fontsize=8); b1.set_xlabel("time from the start (the hold ends at the dashed line)", fontsize=8)
        b0.set_title(f"psi6 per compartment: record means {rw['psi6_rec_0']:.3f} / {rw['psi6_rec_1']:.3f}; Z (hold, last quarter) {rw['Z']:.3f}", fontsize=8)
        plt.setp(b0.get_xticklabels(), visible=False)
        if i == 0: b0.legend(fontsize=7, loc="lower left", framealpha=0.9)
        for a in (a0, a1, a2, a3, b0, b1): a.tick_params(labelsize=7); a.grid(alpha=0.25)
    ref = cs_refs(c["eta"])
    e1 = f"E1 on these {cell['n']} trajectories: nu_d-bar {cell['nud']:.5f}, tau_r {cell['tau_r']:.1f} (the seed-averaged ACF), Q = {cell['Q']:.2f}, nu_0 {cell['nu0']:.5f}, c = {cell['c']:.3f}"
    e1 += (f" +- {cell['sig_c']:.3f}" if math.isfinite(cell["sig_c"]) else " (no jackknife with 2 trajectories)")
    kt = "".join(f"\n{k['source']} (the fetched tables, the same fit): Q = {k['Q']:.2f}, c = {k['c']:.3f} +- {k['sig_c']:.3f}" for k in koa)
    fig.suptitle(f"{c['title']}\n{e1};  {ref['ref']}: c_s {ref['cs']:.3f} (Henderson {ref['cs_henderson']:.3f}){kt}", fontsize=11, y=1.0)
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)


def draw_snaps(path):
    sets = [("69", 0), ("76", 0), ("76", 1)]
    fig, axs = plt.subplots(len(sets), 3, figsize=(17, 5.6 * len(sets))); fig.subplots_adjust(top=0.92, hspace=0.32, wspace=0.10)
    out = []
    for i, (key, r) in enumerate(sets):
        c = [x for x in CELLS if x["key"] == key][0]; m = os.path.join(c["cell"], f"m_{c['M']}")
        hdr, sn = snap_read(glob.glob(os.path.join(m, f"snapshots_L0_*_run{r}.bin.gz"))[0])
        rho = hdr["N"] / (hdr["W"] * hdr["H"]); a_lat = math.sqrt(2 / (math.sqrt(3) * rho)); rc = 1.35 * a_lat
        hold = [s for s in sn if s[1] == 0]; recd = [s for s in sn if s[1] == 1]
        for j, (s, lab) in enumerate(((hold[0], "the first snapshot (hold)"), (hold[-1], "the hold's last"), (recd[-1], "the record's last"))):
            t, ph, cen, xy = s; mag, ang, side = local_psi6(xy, cen, rc); ax = axs[i, j]
            col = plt.cm.hsv(np.mod(ang, 60.0) / 60.0); col[mag < 0.7] = (0.85, 0.85, 0.85, 1.0)
            from matplotlib.collections import EllipseCollection
            ax.add_collection(EllipseCollection(np.full(len(xy), 2 * hdr["R"]), np.full(len(xy), 2 * hdr["R"]), np.zeros(len(xy)), units="xy",
                                                offsets=xy, transOffset=ax.transData, facecolors=col, edgecolors="0.25", linewidths=0.2))
            ax.axvline(cen, color="k", lw=2.0); ax.set_xlim(0, hdr["W"]); ax.set_ylim(0, hdr["H"]); ax.set_aspect("equal")
            fo = [float(np.mean(mag[side == k_])) for k_ in (False, True)]; ordered = [float(np.mean(mag[side == k_] >= 0.7)) for k_ in (False, True)]
            ax.set_title(f"eta {c['eta']}, s = {r}: {lab}, t = {t:.0f}\nlocal |psi6| {fo[0]:.2f} / {fo[1]:.2f} (left / right); ordered {100 * ordered[0]:.0f} % / "
                         f"{100 * ordered[1]:.0f} %\ndivider {cen - c['L0']:+.3f} sigma from L0", fontsize=8)
            ax.tick_params(labelsize=7)
            out.append(dict(cell=key, run=r, which=lab, t=t, divider_minus_L0=cen - c["L0"], local_psi6_left=fo[0], local_psi6_right=fo[1], ordered_left=ordered[0], ordered_right=ordered[1]))
    sm = plt.cm.ScalarMappable(cmap=plt.cm.hsv, norm=plt.Normalize(0, 60)); cb = fig.colorbar(sm, ax=axs, shrink=0.30, pad=0.01)
    cb.set_label("local orientation mod 60 deg (0: lines along x; 30: along y)", fontsize=8)
    fig.suptitle("the crystal per snapshot: each disk coloured by its local bond orientation; grey: local |psi6| < 0.7 (\"ordered\": >= 0.7)\n"
                 "neighbours within 1.35 lattice spacings on the disk's own side of the divider (black line)", fontsize=11, y=0.965)
    fig.savefig(path, dpi=100, bbox_inches="tight"); plt.close(fig)
    return out


def main():
    global OUT
    if "--out" in sys.argv: OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    os.makedirs(OUT); lines = [f"# d56_c_analysis.py: task 2 c of decision 56", ""]
    allrows, cells, koa = [], [], []
    for lab, seeds in (("KOA, the same 3 seeds", [3227069107, 1470420657, 3149991795]), ("KOA, all 25 seeds", None)):
        k = koa_cell(seeds); ref = cs_refs(0.69)
        k.update(cell="69", eta=0.69, source=lab, ref=ref["ref"], cs_ref=ref["cs"], cs_henderson=ref["cs_henderson"], Z_ref=ref["Z"], cT_ref=ref["cT"], c_over_ref=k["c"] / ref["cs"] - 1)
        koa.append(k)
    for c in CELLS:
        lines.append(f"- {c['key']}: {reduce_once(c)}")
        m = os.path.join(c["cell"], f"m_{c['M']}")
        red = pd.read_csv(os.path.join(m, "red_nu.csv")); nud = pd.read_csv(os.path.join(m, "nud_runs.csv")); s1 = pd.read_csv(os.path.join(m, "s1_runs.csv"))
        hl = pd.read_csv(os.path.join(m, "runs_health.csv")); z = np.load(os.path.join(m, "acf_runs.npz")); nu_ref = float(red["nu"].mean())
        ok = nud[nud["fit_ok"] == 1].sort_values("run")
        cell = e1_cell(ok["nu_d"].to_numpy(float), [z[f"run{r}"] for r in ok["run"]], float(red["dt"].iloc[0]), nu_ref, c["L0"], c["alpha"])
        rows, curves = per_run(c, red, nud, s1, hl, z, nu_ref, cell["c"])
        ref = cs_refs(c["eta"]); cell.update(cell=c["key"], eta=c["eta"], source="Mac", ref=ref["ref"], cs_ref=ref["cs"], cs_henderson=ref["cs_henderson"],
                                             Z_ref=ref["Z"], cT_ref=ref["cT"], c_over_ref=cell["c"] / ref["cs"] - 1)
        cells.append(cell); allrows += rows
        draw_cell(c, rows, curves, cell, koa if c["key"] == "69" else [], os.path.join(OUT, f"fig_c{1 if c['key'] == '69' else 2}_N0900_eta{c['eta']:.2f}_traces_fits_spectra_psi6.png"))
    assert sorted(r["seed"] for r in allrows if r["cell"] == "69") == sorted([3227069107, 1470420657, 3149991795])
    cells += koa
    snaps = draw_snaps(os.path.join(OUT, "fig_c3_snapshots_local_psi6.png"))
    R = pd.DataFrame(allrows); Cc = pd.DataFrame(cells); S = pd.DataFrame(snaps)
    R.to_csv(os.path.join(OUT, "c_runs.csv"), index=False, float_format="%.10g"); Cc.to_csv(os.path.join(OUT, "c_cells.csv"), index=False, float_format="%.10g")
    S.to_csv(os.path.join(OUT, "c_snapshots.csv"), index=False, float_format="%.10g")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
    lines += ["", "## per trajectory", "", R.to_string(index=False, float_format=lambda v: f"{v:.5g}"), "", "## per cell (E1 as s1_verdict_261012.per_mass)", "",
              Cc.to_string(index=False, float_format=lambda v: f"{v:.5g}"), "", "## snapshots", "", S.to_string(index=False, float_format=lambda v: f"{v:.4g}")]
    open(os.path.join(OUT, "c_summary.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))


if __name__ == "__main__":
    main()
