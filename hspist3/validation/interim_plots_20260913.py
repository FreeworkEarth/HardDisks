#!/usr/bin/env python3
"""##CHRIS 2026-09-13: interim look at TEST A while it runs. Read-only on the traces, never part of
the overnight pipeline, and limited to 2 worker processes so the running simulations keep their cores.

Figure 1 (per chosen cell): (a) divider position over the first 20 oscillations, (b) the same run over
its whole record, (c) that run's FFT power spectrum with its largest bin marked, (d) every finished
seed's spectrum overlaid with their average, (e) the average as seeds are added, (f) each seed's
largest bin on a 25- versus a 1000-oscillation record.
Figure 2: c_s against record length from the finished trajectories only.
"""
import glob
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tests_20260913 as T  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.figure as mfig  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, BLUE_LIGHT, ORANGE = "#2a78d6", "#b7d3f6", "#eb6834"
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # sequential blue, ordinal steps 250..700
STAMP = time.strftime("%Y-%m-%d %H:%M %Z")


def style(ax, title):
    ax.set_facecolor(SURF)
    for sp in ax.spines.values():
        sp.set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", color=INK, fontsize=10)


def labels(ax, xl, yl):
    ax.set_xlabel(xl, color=INK2, fontsize=9)
    ax.set_ylabel(yl, color=INK2, fontsize=9)


def cell_dir(L0s, M):
    return os.path.join(T.TROOT, "A_length", f"L0_{T.tag(L0s)}", f"m_{M}")


def finished(L0s, M):
    return [q for q in T.cell_runs(cell_dir(L0s, M), M) if not q[2]]


def spec_at(t, x, nup, N):
    n = T._prefix(t, nup, N)
    dt = (t[-1] - t[0]) / (len(t) - 1)
    P, df = T._spectrum(x[:n], dt)
    k, nu, _ = T._peak(P, df)
    return P, df, k, nu


def plot_cell(L0s, L0, M, out):
    runs = finished(L0s, M)
    data = [(r, *T._load(p)) for r, p, _ in runs]
    R = len(data)
    eta = T.eta_of(L0)
    fig = plt.figure(figsize=(13.5, 13), facecolor=SURF)
    gs = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.24, left=0.07, right=0.97, top=0.91, bottom=0.05)
    r0, t, x, nup = data[0]
    ph = t * nup
    # frequency scale for every zoom: the 25-oscillation median peak, which no heavy divider takes at the drift
    f_res = float(np.median([spec_at(tt, xx, nn, 25)[3] for _, tt, xx, nn in data]))

    ax = fig.add_subplot(gs[0, 0])
    style(ax, "(a) divider position, first 20 oscillations after release")
    m = ph <= 20
    ax.plot(ph[m], x[m], color=BLUE, lw=1.0)
    ax.axhline(0, color=AXIS, lw=0.8)
    labels(ax, "time, in predicted oscillation periods", "displacement from centre  [σ]")

    ax = fig.add_subplot(gs[0, 1])
    style(ax, f"(b) the same run, whole record of {ph[-1]:.0f} oscillations")
    step = max(1, len(x) // 60000)
    ax.plot(ph[::step], x[::step], color=BLUE, lw=0.3)
    ax.axhline(0, color=AXIS, lw=0.8)
    labels(ax, "time, in predicted oscillation periods", "displacement from centre  [σ]")

    specs = [spec_at(tt, xx, nn, 1000) for _, tt, xx, nn in data]
    nus = np.array([s[3] for s in specs])
    nubar = float(nus.mean())
    P0, df0, k0, nu0 = specs[0]
    f0 = np.arange(len(P0)) * df0

    ax = fig.add_subplot(gs[1, 0])
    style(ax, "(c) FFT power spectrum of that run, whole record")
    band = (f0 > 0) & (f0 <= 2.2 * f_res)
    ax.semilogy(f0[band], P0[band], color=BLUE, lw=0.6)
    ax.plot([nu0], [P0[k0]], "o", color=ORANGE, ms=8, mec=SURF, mew=1.3, zorder=5)
    ax.annotate(f"largest bin  ν = {nu0:.5f}\n(bin {k0}, bin width {df0:.2e})", xy=(nu0, P0[k0]),
                xytext=(0.62, 0.9), textcoords="axes fraction", color=INK2, fontsize=8.5, va="top",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8))
    labels(ax, "frequency f", "power |FFT|²")

    L = min(len(s[0]) for s in specs)
    f = np.arange(L) * specs[0][1]
    stack = np.array([s[0][:L] for s in specs])
    Pm = stack.mean(axis=0)
    nures = float(np.median(nus))
    n_drift = int(sum(1 for s in specs if s[2] <= 3))
    zb =(f >= 0.55 * f_res) & (f <= 1.6 * f_res)

    ax = fig.add_subplot(gs[1, 1])
    style(ax, f"(d) all {R} finished seeds overlaid, and their average")
    for row in stack:
        ax.semilogy(f[zb], row[zb], color=BLUE_LIGHT, lw=0.5, alpha=0.85, zorder=1)
    ax.semilogy(f[zb], Pm[zb], color=BLUE, lw=2.2, zorder=3)
    trans = ax.get_xaxis_transform()
    for v in nus:
        ax.plot([v, v], [0.0, 0.06], color=INK, lw=0.9, transform=trans, zorder=4)
    if f[zb][0] <= nures <= f[zb][-1]:
        ax.axvline(nures, color=INK, lw=1.0, ls=":", zorder=2)
    if f[zb][0] <= nubar <= f[zb][-1]:
        ax.axvline(nubar, color=INK, lw=1.0, ls="--", zorder=2)
    ax.set_xlim(f[zb][0], f[zb][-1])
    ax.text(0.02, 0.97, "window centred on the 25-oscillation median peak\nthin: each seed\nbold: average of all seeds\n"
            "ticks at the bottom: each seed's largest bin\n"
            f"dotted: median peak {nures:.5f};  dashed: mean peak {nubar:.5f}\n"
            f"{n_drift} of {R} seeds had their largest bin at the slow drift (bins 1-3, off the left edge)",
            transform=ax.transAxes, color=INK2, fontsize=8.3, va="top", zorder=6,
            bbox=dict(fc=SURF, ec="none", alpha=0.85, pad=2))
    labels(ax, "frequency f", "power |FFT|²")

    ax = fig.add_subplot(gs[2, 0])
    style(ax, "(e) the average spectrum as more seeds are added")
    counts = [c for c in (1, 5, 10, 25) if c <= R]
    if R not in counts:
        counts.append(R)
    counts = counts[-4:]
    for c, col in zip(counts, RAMP[-len(counts):]):
        ax.semilogy(f[zb], stack[:c].mean(axis=0)[zb], color=col, lw=1.4 if c < R else 2.2, label=f"{c} seed{'s' if c > 1 else ''}")
    if f[zb][0] <= nubar <= f[zb][-1]:
        ax.axvline(nubar, color=INK, lw=1.0, ls="--")
    ax.set_xlim(f[zb][0], f[zb][-1])
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper right")
    ax.text(0.02, 0.97, "random bin-to-bin noise shrinks as 1/√(seeds);\nthe resonance peak stays where it is",
            transform=ax.transAxes, color=INK2, fontsize=8.3, va="top")
    labels(ax, "frequency f", "power |FFT|², averaged")

    ax = fig.add_subplot(gs[2, 1])
    style(ax, "(f) each seed's largest bin: 25 versus 1000 oscillations")
    allv = []
    for row, N, col in ((1, 25, ORANGE), (0, 1000, BLUE)):
        sp = [spec_at(tt, xx, nn, N) for _, tt, xx, nn in data]
        v = np.array([s[3] for s in sp])
        dfN = sp[0][1]
        allv += list(v)
        jitter = ((np.arange(len(v)) % 7) - 3) * 0.045
        ax.scatter(v, row + jitter, s=30, color=col, edgecolor=SURF, lw=0.8, zorder=3)
        ax.plot([v.mean()] * 2, [row - 0.3, row + 0.3], color=INK, lw=1.6, ls="--", zorder=4)
        ax.plot([np.median(v)] * 2, [row - 0.3, row + 0.3], color=INK, lw=1.6, ls=":", zorder=4)
        ax.text(0.01, row + 0.42, f"{N} oscillations: mean {v.mean():.5f}, median {np.median(v):.5f}, "
                f"one bin = {100 * dfN / f_res:.1f} % of the oscillation frequency",
                transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData),
                ha="left", va="bottom", color=INK2, fontsize=8.3)
        if N == 25:
            lo, hi = min(v), max(v)
            kk = np.arange(max(1, int(lo / dfN) - 1), int(hi / dfN) + 2)
            for kb in kk:
                ax.axvline(kb * dfN, color=GRID, lw=0.8, zorder=1)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["1000 oscillations", "25 oscillations"], color=INK2, fontsize=9)
    ax.set_ylim(-0.6, 1.9)
    pad = 0.05 * (max(allv) - min(allv) + 1e-12)
    ax.set_xlim(min(allv) - pad, max(allv) + pad)
    ax.set_xlabel("frequency of the largest bin   (dashed = mean, dotted = median, grey lines = 25-oscillation bins)",
                  color=INK2, fontsize=9)

    fig.suptitle(f"TEST A, interim {STAMP}:  L0 = {L0:g} σ (η = {eta:.3f}),  divider mass M = {M},  {R} finished seeds\n"
                 "Román 2002 estimator: FFT of the divider position from release, largest nonzero bin per seed, mean over seeds",
                 color=INK, fontsize=11)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=170, facecolor=SURF)
    plt.close(fig)
    return R, nubar


def interim_cs(out):
    tasks, used = [], 0
    masses = []
    for L0s, L0 in T.L0S:
        for M in T.ROMAN_MASSES:
            runs = T.cell_runs(cell_dir(L0s, M), M)
            if len([q for q in runs if not q[2]]) >= 3:
                tasks.append(((L0, M), L0, M, runs, T.CUTS))
                used += len(runs)
                masses.append((L0, M))
    if not tasks:
        return None
    results = T.analyse(tasks, workers=2)
    rows = []
    for _, L0 in T.L0S:
        rr, _ = T.summarize(results, L0, T.CUTS, T.ROMAN_MASSES)
        rows += rr
    T.write_csv(f"{out}.csv", rows)
    per_L0 = {L0: sorted(M for l, M in masses if l == L0) for _, L0 in T.L0S}
    title = (f"TEST A, interim {STAMP}: c_s against record length, finished trajectories only ({used} of 600)\n"
             f"masses with ≥ 3 finished seeds: L0 = 20: {per_L0[20.0]};  L0 = 7.5: {per_L0[7.5]}.  "
             "Error bar on the Román mean = 1σ scatter of per-mass c_s.")
    orig = mfig.Figure.suptitle
    mfig.Figure.suptitle = lambda self, t, **kw: orig(self, title, **kw)
    try:
        T.figure_A(rows, T.CUTS, out)
    finally:
        mfig.Figure.suptitle = orig
    return rows, used, per_L0


if __name__ == "__main__":
    P = T.PLOTS
    done = []
    for L0s, L0 in T.L0S:
        for M in T.ROMAN_MASSES:
            if M > 100 and len(finished(L0s, M)) >= 5:
                R, nb = plot_cell(L0s, L0, M, os.path.join(P, f"260913_testA_interim_divider_spectra_L0_{L0:g}_M{M}".replace(".", "p")))
                done.append(f"L0={L0:g} M={M}: {R} seeds, mean peak {nb:.5f}")
    light = None
    for L0s, L0 in T.L0S:
        for M in (20, 50, 100):
            if len(finished(L0s, M)) >= 5:
                light = (L0s, L0, M)
                break
        if light:
            break
    if light:
        R, nb = plot_cell(light[0], light[1], light[2],
                          os.path.join(P, f"260913_testA_interim_divider_spectra_L0_{light[1]:g}_M{light[2]}".replace(".", "p")))
        done.append(f"L0={light[1]:g} M={light[2]}: {R} seeds, mean peak {nb:.5f}")
    else:
        done.append("no light divider (M <= 100) has 5 finished seeds yet")
    ready = {L0: sum(1 for M in T.ROMAN_MASSES if len(finished(L0s, M)) >= 3) for L0s, L0 in T.L0S}
    res = None
    if max(ready.values()) >= 3:
        res = interim_cs(os.path.join(P, "260913_testA_interim_cs_vs_record_length"))
    else:
        done.append(f"interim c_s skipped: masses with >= 3 finished seeds per L0 = {ready} (the slope needs 3)")
    if res:
        rows, used, per_L0 = res
        done.append(f"interim c_s from {used} trajectories, masses {per_L0}")
        for r in rows:
            done.append(f"   L0={r['L0']:g} cut={r['cut']:>4}: Román mean {T.fmt(r['roman_mean'])} ± {T.fmt(r['roman_mean_scatter'])}  "
                        f"median {T.fmt(r['roman_median'])}  parabolic {T.fmt(r['parabolic'])}  fit {T.fmt(r['spectrum_fit'])}  masses {r['n_masses']}")
    print("\n".join(done))
