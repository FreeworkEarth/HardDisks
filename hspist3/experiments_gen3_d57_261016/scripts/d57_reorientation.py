#!/usr/bin/env python3
"""##CHRIS 2026-10-16 (plan-author decision 57, addendum 3, task 4; notes 261012 sec. 4.7.140): THE REORIENTATION AT eta 0.72.
(After the independent verification of Report D57-3: the snapshot interval in the text is read from the data, 67-117 at 0.72, not the
61-117 of both etas; the Mac-run table keeps its t = 150 row.)
  1. P2 (KOA): the two pin-0 cells at eta 0.72 (G1 H 30, G2 H 29.708), all 50 runs each: the psi6 orientation per compartment against time
     (psi6_t, every Delta_psi = 1-3.4 sigma-time; drawn in [-15, 45) deg so that 0 = lines along x and 30 = along y lie away from the wrap),
     the hold's end marked; and per compartment the time of the first SNAPSHOT with the lines along y (p2_frames.csv; snapshots every
     67-117 sigma-time at 0.72: the direction from theta0, y when theta0 mod 60 lies within 15 deg of 30) and of the first psi6_t sample within
     7.5 deg of 30; the number of compartments with that time before 200, 300 and 500 sigma-time.
  2. The ONE Mac run (experiments_gen3_d57_261016/data: P2's line 21, (0.72, G1, pin 0, alpha 5, s 0), snapshots every 10 sigma-time):
     the frames of its first 600 sigma-time, every disk coloured by its local psi6 orientation (neighbours within 1.35 lattice spacings on
     its own side of the divider; mod 60; grey when |psi6_local| < 0.5), with each compartment's theta0, line direction and ordered
     fraction (p2_lattice.config_strain) per frame; and its psi6_t orientation per compartment.
Writes reorientation/: reorientation_p2.csv (per compartment), reorientation_mac_frames.csv (per frame), reorientation.md,
fig_d57_reorientation_p2.png, fig_d57_reorientation_frames.png.
usage (from hspist3/): python3 experiments_gen3_d57_261016/scripts/d57_reorientation.py [--out <a new folder>]
"""
import glob, math, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); E = os.path.join(HS, "experiments_gen3_d57_261016")
sys.path.insert(0, os.path.join(HS, "cluster", "gen3_stage2_261016"))
import p2_lattice as LT
P2 = os.path.join(HS, "experiments_gen3_stage2_koa_261016", "P2", "merged", "stage2", "gen3")
CELLS = (("P2_N0900_eta0.720000_G1_pin0", "G1, H 30"), ("P2_N0900_eta0.720000_G2_pin0", "G2, H 29.708"))
MAC = os.path.join(E, "data", "stage2", "gen3", "P2_N0900_eta0.720000_G1_pin0", "m_4500")
FRAMES_T = (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 300, 600)
OUT = os.path.abspath(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else os.path.join(E, "reorientation")
TEQ, LIMITS = 5000.0, (200, 300, 500)


def shown(a): return np.mod(np.asarray(a, float) + 15.0, 60.0) - 15.0          # [-15, 45): 0 = x, 30 = y


def line_of(theta):
    t = np.mod(np.asarray(theta, float), 60.0); return np.where(np.isnan(t), np.nan, np.where(np.abs(t - 30.0) < 15.0, 90.0, 0.0))


def local_psi6(xy, cen, W, H, th=LT.TH):
    """per disk: |psi6| and orientation (deg, mod 60) over its neighbours within 1.35 a on its own side (a from the side's density)"""
    out_abs, out_or = np.full(len(xy), np.nan), np.full(len(xy), np.nan)
    for side, width in ((xy[:, 0] < cen, cen - th / 2), (xy[:, 0] > cen, W - cen - th / 2)):
        q = xy[side]
        a = math.sqrt(2.0 / (math.sqrt(3.0) * len(q) / (width * H))); d = q[None, :, :] - q[:, None, :]; r = np.hypot(d[..., 0], d[..., 1])
        nb = (r < 1.35 * a) & (r > 0); z = np.where(nb, np.exp(6j * np.arctan2(d[..., 1], d[..., 0])), 0).sum(axis=1)
        psi = np.where(nb.sum(axis=1) > 0, z / np.maximum(nb.sum(axis=1), 1), 0)
        out_abs[side], out_or[side] = np.abs(psi), np.mod(np.degrees(np.angle(psi)) / 6.0, 60.0)
    return out_abs, out_or


def main():
    if os.path.exists(OUT): sys.exit(f"REFUSED: {OUT} exists (a new folder only)")
    os.makedirs(OUT); rows, series = [], {}
    for cid, lab in CELLS:
        for d in sorted(glob.glob(os.path.join(P2, cid, "m_*"))):
            M = int(d.split("_")[-1]); F = pd.read_csv(os.path.join(d, "p2_frames.csv")); R = pd.read_csv(os.path.join(d, "p2_runs.csv"))
            for r in R.to_dict("records"):
                ps = pd.read_csv(glob.glob(os.path.join(d, f"psi6_t_L0_*_run{int(r['run'])}.csv*"))[0]); f = F[F["run"] == r["run"]].sort_values("t")
                for k in (0, 1):
                    o = shown(ps[f"psi6_orient_deg_{k}"]); series.setdefault(cid, []).append((ps["t_sigma"].to_numpy(), o))
                    ln = line_of(f[f"theta0_{k}"]); fy = f["t"].to_numpy()[np.where(ln == 90)[0]]
                    near = np.where(np.abs(np.asarray(o) - 30.0) <= 7.5)[0]
                    rows.append(dict(cell=cid, M=M, run=int(r["run"]), side=k, start_line=r[f"start_line_{k}"], end_line=r[f"end_line_{k}"],
                                     first_y_snapshot=float(fy[0]) if len(fy) else math.nan, snapshot_every=float(f["t"].diff().median()),
                                     first_psi6t_near30=float(ps["t_sigma"].iloc[near[0]]) if len(near) else math.nan, psi6t_every=float(ps["t_sigma"].diff().median())))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "reorientation_p2.csv"), index=False, float_format="%.6g")
    # the Mac run's frames
    S = LT.read_snapshots(glob.glob(os.path.join(MAC, "snapshots_*_run0.bin*"))[0]); fr_rows, frames = [], {}
    for fr in S["frames"]:
        if fr["t"] > 600.5: break
        res = LT.config_strain(fr["xy"], fr["cen"], S["W"], S["H"])
        fr_rows.append(dict(t=fr["t"], phase=fr["phase"], **{f"{k}_{s}": res[i][k] for i, s in ((0, "left"), (1, "right")) for k in ("theta0", "line", "ordered_bulk", "e")}))
        for t in FRAMES_T:
            if abs(fr["t"] - t) < 0.5 and t not in frames: frames[t] = fr       # the snapshot times are k x 10 to rounding
    MF = pd.DataFrame(fr_rows); MF.to_csv(os.path.join(OUT, "reorientation_mac_frames.csv"), index=False, float_format="%.6g")
    mps = pd.read_csv(glob.glob(os.path.join(MAC, "psi6_t_L0_*_run0.csv*"))[0])
    # the report
    L = ["# The reorientation at eta 0.72 (decision 57 addendum 3, task 4; d57_reorientation.py)", "",
         "## P2 (KOA): the two pin-0 cells, 50 runs x 2 compartments each", "",
         "| cell | compartments | start x / y | end x / y | snapshots every [sigma] | psi6_t every [sigma] | first snapshot with y: median (range) | "
         + " | ".join(f"first y snapshot before {t}" for t in LIMITS) + " | first psi6_t within 7.5 deg of 30: median (range) | "
         + " | ".join(f"before {t}" for t in LIMITS) + " |", "|---|---|---|---|---|---|---|" + "---|" * 3 + "---|" + "---|" * 3]
    for cid, lab in CELLS:
        g = D[D["cell"] == cid]; fy = g["first_y_snapshot"].dropna(); fp = g["first_psi6t_near30"].dropna()
        L.append(f"| {lab} | {len(g)} | {int((g['start_line'] == 0).sum())} / {int((g['start_line'] == 90).sum())} | {int((g['end_line'] == 0).sum())} / "
                 f"{int((g['end_line'] == 90).sum())} | {g['snapshot_every'].min():.0f}-{g['snapshot_every'].max():.0f} | {g['psi6t_every'].min():.1f}-{g['psi6t_every'].max():.1f} | "
                 f"{fy.median():.0f} ({fy.min():.0f} to {fy.max():.0f}) | " + " | ".join(f"{int((fy < t).sum())} of {len(g)}" for t in LIMITS)
                 + f" | {fp.median():.0f} ({fp.min():.0f} to {fp.max():.0f}) | " + " | ".join(f"{int((fp < t).sum())} of {len(g)}" for t in LIMITS) + " |")
    first_snap = D.groupby("cell").apply(lambda g: int((g["first_y_snapshot"] <= g["snapshot_every"] * 1.5).sum()))
    L += ["", "- The first snapshot with y is the first snapshot after the start (t = one snapshot interval) in "
              + ", ".join(f"{lab}: {first_snap[cid]} of 100" for cid, lab in CELLS)
              + f" compartments: the snapshots every {D['snapshot_every'].min():.0f}-{D['snapshot_every'].max():.0f} sigma-time do not show the rotation itself, "
              f"only that it is over; psi6_t (every {D['psi6t_every'].min():.1f}-{D['psi6t_every'].max():.1f} sigma-time) shows the orientation's path but not "
              "the configurations.", "",
          "## The Mac run (snapshots every 10 sigma-time; P2's line 21 with T_snap 10)", "",
          "| t [sigma] | left: theta0 [deg], lines, ordered fraction, e [%] | right: theta0 [deg], lines, ordered fraction, e [%] |", "|---|---|---|"]
    for r in MF[MF["t"] <= 150.5].to_dict("records"):
        L.append(f"| {r['t']:.0f} | {r['theta0_left']:.1f}, {'y' if r['line_left'] == 90 else 'x'}, {r['ordered_bulk_left']:.2f}, {100 * r['e_left']:+.2f} | "
                 f"{r['theta0_right']:.1f}, {'y' if r['line_right'] == 90 else 'x'}, {r['ordered_bulk_right']:.2f}, {100 * r['e_right']:+.2f} |")
    for k, s in ((0, "left"), (1, "right")):
        ys = MF[MF[f"line_{s}"] == 90]["t"]; o = shown(mps[f"psi6_orient_deg_{k}"]); near = np.where(np.abs(o - 30) <= 7.5)[0]
        back = MF[(MF["t"] > ys.iloc[0]) & (MF[f"line_{s}"] == 0)] if len(ys) else MF.iloc[0:0]
        L.append(f"\n- the {s} compartment: the first 10-sigma snapshot with lines along y at t = {ys.iloc[0]:.0f}; the first psi6_t sample within 7.5 deg of 30 at "
                 f"t = {mps['t_sigma'].iloc[near[0]]:.1f}; snapshots with lines along x after that, up to 600: {len(back)}")
    open(os.path.join(OUT, "reorientation.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    # figure 1: P2's 200 compartments
    fig, axs = plt.subplots(2, 2, figsize=(16, 8.6), gridspec_kw=dict(width_ratios=[1.6, 1.0]))
    for i, (cid, lab) in enumerate(CELLS):
        for t, o in series[cid]:
            for ax, lim in ((axs[i, 0], None), (axs[i, 1], 600)):
                k = slice(None) if lim is None else (t <= lim); ax.plot(t[k], np.asarray(o)[k], color="#1f3b5a", lw=0.25, alpha=0.25)
        for ax in axs[i]:
            ax.axhline(0, color="#2a6fbb", lw=0.9, ls=":"); ax.axhline(30, color="#d9822b", lw=0.9, ls=":"); ax.set_ylim(-15, 45); ax.grid(alpha=0.2)
            ax.set_ylabel("psi6 orientation [deg]  (0: lines along x, 30: along y)", fontsize=8); ax.tick_params(labelsize=8)
        axs[i, 0].axvline(TEQ, color="0.2", ls="--", lw=1.0); axs[i, 0].text(TEQ, 42, " the hold's end", fontsize=8, va="top")
        axs[i, 0].set_title(f"eta 0.72, {lab}, start pin 0: 50 runs x 2 compartments, the whole run", fontsize=9)
        axs[i, 1].set_title(f"{lab}: the first 600 sigma-time", fontsize=9)
    for ax in axs[1]: ax.set_xlabel("time from the start  [sigma-time]", fontsize=9)
    fig.suptitle("P2 (KOA): the psi6 orientation per compartment of the two eta 0.72 pin-0 cells (psi6_t)", fontsize=10.5)
    fig.savefig(os.path.join(OUT, "fig_d57_reorientation_p2.png"), dpi=110, bbox_inches="tight"); plt.close(fig)
    # figure 2: the Mac run frame by frame
    ts = [t for t in FRAMES_T if t in frames]; nc = 5; nr_ = math.ceil(len(ts) / nc)
    fig = plt.figure(figsize=(3.4 * nc, 3.3 * nr_ + 3.0)); gs = fig.add_gridspec(nr_ + 1, nc, height_ratios=[0.9] + [1.0] * nr_, hspace=0.32, wspace=0.08)
    top = fig.add_subplot(gs[0, :]); cmap = plt.get_cmap("twilight")
    for k, s, col in ((0, "left", "#2a6fbb"), (1, "right", "#d9822b")):
        m = mps["t_sigma"] <= 600; top.plot(mps["t_sigma"][m], shown(mps[f"psi6_orient_deg_{k}"])[m], color=col, lw=0.8, label=f"{s} compartment (psi6_t, every 1.7 sigma)")
        top.plot(MF["t"], shown(MF[f"theta0_{s}"]), "o", color=col, ms=3, mfc="white")
    for t in ts: top.axvline(t, color="0.75", lw=0.6)
    top.set_ylim(-15, 45); top.set_xlim(0, 600); top.axhline(0, color="0.5", lw=0.6, ls=":"); top.axhline(30, color="0.5", lw=0.6, ls=":")
    top.set_ylabel("orientation [deg]", fontsize=8); top.set_xlabel("time from the start [sigma-time] (open circles: theta0 of the 10-sigma snapshots; grey lines: "
                                                                  "the frames below)", fontsize=8)
    top.legend(fontsize=7, loc="center right"); top.tick_params(labelsize=7)
    for j, t in enumerate(ts):
        fr = frames[t]; ax = fig.add_subplot(gs[1 + j // nc, j % nc]); pa, po = local_psi6(fr["xy"], fr["cen"], S["W"], S["H"])
        ok = pa >= 0.5
        ax.scatter(fr["xy"][~ok, 0], fr["xy"][~ok, 1], s=7, color="0.8", lw=0)
        sc = ax.scatter(fr["xy"][ok, 0], fr["xy"][ok, 1], s=7, c=po[ok], cmap=cmap, vmin=0, vmax=60, lw=0)
        ax.axvline(fr["cen"], color="k", lw=0.8); ax.set_xlim(0, S["W"]); ax.set_ylim(0, S["H"]); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        r = MF.iloc[(MF["t"] - t).abs().idxmin()]
        ax.set_title(f"t = {t}: L {'y' if r['line_left'] == 90 else 'x'} ({r['ordered_bulk_left']:.2f}), R {'y' if r['line_right'] == 90 else 'x'} ({r['ordered_bulk_right']:.2f})",
                     fontsize=8)
    cb = fig.colorbar(sc, ax=fig.axes[1:], fraction=0.015, pad=0.01); cb.set_label("local psi6 orientation [deg] (0 = 60: lines along x; 30: along y); grey: |psi6| < 0.5", fontsize=8)
    fig.suptitle("The Mac run: P2's line 21 (eta 0.72, G1 H 30, pin 0, alpha 5, s 0) with snapshots every 10 sigma-time; each frame's title: the line "
                 "direction per compartment (its ordered fraction)", fontsize=10)
    fig.savefig(os.path.join(OUT, "fig_d57_reorientation_frames.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
