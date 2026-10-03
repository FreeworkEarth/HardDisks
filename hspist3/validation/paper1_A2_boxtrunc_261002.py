#!/usr/bin/env python3
"""##CHRIS 2026-10-02: methods sec. 14.3 -- the box-truncation correction applied to the A2 size ladder, by the
pre-registered rule (commit 0ddefa9). No trajectory is re-analysed: the per-mass frequencies are data, the correction
changes only x_M = K(alpha)/(2 pi L_eff) by one factor per (eta, N) cell, and the eta at which KR is evaluated.

Geometries: P = as published (L_eff = L0 - 2r), B = before (L0 - 2r - t/2, the canonical geometry), A = after
(L0 - 2r - t/2 - delta/2, eta_true). Verdict: B vs A only.
  per point : D = (c_s - KR(eta))/sigma, sigma = the plotted error (slope_with_errors, sd/sqrt(n), sqrt(chi2_red));
              REGENERATE if |D_A - D_B| > 0.5 for any point of either table.
  fit       : c_s = c_inf + b/sqrt(N), weights = SD of nu/x (analyze_A2_X2p5_20260914.py:120), on the draft's source
              table 260919_A2_cs_per_mass.csv; in A on the deviations from KR(eta_true,N); REGENERATE if c_inf or b
              moves by more than its own (B) 1 sigma.
usage: python3 hspist3/validation/paper1_A2_boxtrunc_261002.py [--write]   (--write: only after REGENERATE, see sec. 14.3)
"""
import csv, glob, math, os, re, sys
from collections import defaultdict
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos
from paper1_populate_cs_err_20261002 import slope_with_errors
R, TW = 0.5, T.WALL_T
TABLES = ("260919_A2_cs_per_mass.csv", "260919_A2_cs_per_mass_famB.csv")
PUB = "260917_A2_cs_vs_N_extrapolation.csv"
A2ROOT = os.path.join(T.ROOT)          # .../00_eta_sweep_ROMAN, holds the A2_* campaign roots

def pre(name):                          # the uncorrected input: the dated copy once it exists (re-runs never double-correct)
    p = name.replace(".csv", "_pre_boxtrunc_261002.csv")
    return p if os.path.exists(T.plot_path(p)) else name

def delta(L0):                          # the binary's float32 SIM_WIDTH (00ALLINONE.c:323), as in sec. 14.2
    w = np.float32(2) * np.float32(L0) * np.float32(24)
    return (float(w) - math.floor(float(w))) / 24.0

ROOTS = ("A2_*", "famB_20260911")      # famB_20260911 holds the N = 100/400 ladder cells (not the _ABANDONED_seedpad copy)

def height(eta, N):
    """H from --height= in the runs' 00_COMMAND.md; else read back from a trace's recorded eta and L0."""
    tag = f"eta_{eta:.2f}".replace(".", "p")
    cmds = [c for r in ROOTS for c in glob.glob(os.path.join(A2ROOT, r, tag, f"N{N}", "m_*", "**", "00_COMMAND.md"), recursive=True)]
    for c in sorted(cmds):
        m = re.search(r"--height=([\d.]+)", open(c, errors="ignore").read())
        if m: return float(m.group(1)), "command"
    trs = [t for r in ROOTS for t in glob.glob(os.path.join(A2ROOT, r, tag, f"N{N}", "m_*", "**", "wall_x_positions_*.csv"), recursive=True)]
    for tr in sorted(trs)[:1]:
        import pandas as pd
        h = pd.read_csv(tr, nrows=1)
        if "eta" in h and "L0" in h:
            return N * math.pi * R * R / (2 * float(h["L0"].iloc[0]) * float(h["eta"].iloc[0])), "read back"
    return float("nan"), "not found"

def x_of(M, N, Le): return T.k_root(M / float(N)) / (2 * math.pi * Le)

def cells(table):
    by = defaultdict(list)
    for r in csv.DictReader(open(T.plot_path(pre(table)))):
        by[(float(r["eta"]), int(r["N"]))].append(dict(L0=float(r["L0"]), M=int(r["M"]), n=int(r["n_runs"]),
                                                         nu=float(r["nu_mean"]), sd=float(r["nu_sd"])))
    return by

def point(rows, eta, N, geom):
    L0 = rows[0]["L0"]; d = delta(L0)
    Le = {"P": L0 - 2 * R, "B": L0 - 2 * R - TW / 2, "A": L0 - 2 * R - TW / 2 - d / 2}[geom]
    e = eta * L0 / (L0 - d / 2) if geom == "A" else eta
    x = np.array([x_of(q["M"], N, Le) for q in rows]); nu = np.array([q["nu"] for q in rows])
    sy = np.array([q["sd"] / math.sqrt(max(1, q["n"])) for q in rows])
    cs = float((x * nu).sum() / (x * x).sum()); sc = float(np.std(nu / x, ddof=1))
    _s, _e, es, _x2 = slope_with_errors(x, nu, sy)
    kr = float(T.kr_cs(e))
    return dict(L0=L0, d=d, Le=Le, eta=e, cs=cs, sc=sc, es=float(es), kr=kr, D=(cs - kr) / float(es))

def fit(pts, eta, geom):
    Ns = np.array(sorted(N for (e, N) in pts if e == eta), float)
    p = [pts[(eta, int(N))][geom] for N in Ns]
    dv = np.array([q["cs"] - q["kr"] for q in p]); s = np.array([max(q["sc"], 1e-12) for q in p])
    sl, sl_err, ic, ic_err = sos.weighted_linreg(1 / np.sqrt(Ns), dv, s, force_zero_intercept=False)[:4]
    chi2 = float((((dv - (ic + sl / np.sqrt(Ns))) / s) ** 2).sum()); kr0 = float(T.kr_cs(eta))
    return dict(Ns=Ns, cinf=kr0 + ic, cinf_err=ic_err, b=sl, b_err=sl_err, chi2=chi2, dof=len(Ns) - 2, kr=kr0)

def main():
    allpts = {}
    for tab in TABLES:
        by = cells(tab); allpts[tab] = {k: {g: point(v, k[0], k[1], g) for g in "PBA"} for k, v in by.items()}
    # ---- estimator gate: the published fit, reproduced from the draft's source table in geometry P
    pub = {round(float(r["eta"]), 2): r for r in csv.DictReader(open(T.plot_path(pre(PUB))))}
    print(f"inputs (uncorrected): {', '.join(pre(t) for t in TABLES)}, {pre(PUB)}\n")
    src = allpts[TABLES[0]]; fits = {}
    print(f"### Estimator gate: {PUB} reproduced from {TABLES[0]} (geometry P, analyze_A2_X2p5_20260914.py definitions)\n")
    print("| eta | N values | c_inf (pub) | c_inf (here) | ± (pub/here) | b (pub/here) | ± (pub/here) | chi2 (pub/here) | reproduced |")
    print("|---|---|---|---|---|---|---|---|---|")
    gate = True
    for eta in sorted({e for e, _ in src}):
        if sum(1 for (e, _) in src if e == eta) < 3: continue
        fits[eta] = {g: fit(src, eta, g) for g in "PBA"}; f = fits[eta]["P"]; p = pub.get(round(eta, 2))
        ok = p is not None and all(f"{a:.5f}" == p[k] for a, k in ((f["cinf"], "c_inf"), (f["cinf_err"], "c_inf_err"), (f["b"], "slope"), (f["b_err"], "slope_err"))) \
             and f"{f['chi2']:.3f}" == p["chi2"]
        gate &= ok
        print(f"| {eta:.2f} | {' '.join(str(int(n)) for n in f['Ns'])} | {p['c_inf'] if p else '-'} | {f['cinf']:.5f} | {p['c_inf_err'] if p else '-'} / {f['cinf_err']:.5f} | "
              f"{p['slope'] if p else '-'} / {f['b']:.5f} | {p['slope_err'] if p else '-'} / {f['b_err']:.5f} | {p['chi2'] if p else '-'} / {f['chi2']:.3f} | {'yes' if ok else '**NO**'} |")
    print(f"\nestimator gate: {'PASS -- every published row reproduced to its printed digits' if gate else 'FAIL -- VOID, stop'}")
    if not gate: sys.exit(1)

    # ---- geometry and per-point table
    print("\n### Per ladder point (sigma = plotted error; KR at the point's own eta)\n")
    print("| table | eta_rec | N | L_0 | H (source) | H*24 | delta | eta_true | L_eff,rec (B) | L_eff,true (A) | c_s before ± σ | c_s after ± σ | KR(eta_rec) | KR(eta_true) | D before | D after | change |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    moved, Hbad, maxd = [], [], 0.0
    for tab in TABLES:
        for (eta, N), g in sorted(allpts[tab].items()):
            b, a = g["B"], g["A"]; H, hs = height(eta, N); hp = H * 24
            if not (abs(hp - round(hp)) < 1e-6): Hbad.append((tab, eta, N, H))
            ch = a["D"] - b["D"]; maxd = max(maxd, abs(ch))
            if abs(ch) > 0.5: moved.append((tab, eta, N, ch))
            print(f"| {'A2' if tab == TABLES[0] else 'famB'} | {eta:.2f} | {N} | {b['L0']:.6f} | {H:.5f} ({hs}) | {hp:.3f} | {a['d']:.6f} | {a['eta']:.6f} | {b['Le']:.6f} | "
                  f"{a['Le']:.6f} | {b['cs']:.5f} ± {b['es']:.5f} | {a['cs']:.5f} ± {a['es']:.5f} | {b['kr']:.5f} | {a['kr']:.5f} | {b['D']:+.2f} | {a['D']:+.2f} | {ch:+.2f} |")
    print(f"\nheight cast: H*24 integer in {'every cell' if not Hbad else 'all but ' + str(Hbad)}")
    print(f"largest |change in D|: {maxd:.2f}; points with |change| > 0.5: {len(moved)} {[(('A2' if t == TABLES[0] else 'famB'), e, n, round(c, 2)) for t, e, n, c in moved]}")

    # ---- fit table
    print(f"\n### Finite-size fit c_s = c_inf + b/sqrt(N), weights = mass scatter, source {TABLES[0]}\n")
    print("| eta | N values | P: c_inf ± | P: b ± | B: c_inf ± | B: b ± | A: c_inf ± | A: b ± | (A-B) c_inf / σ_B | (A-B) b / σ_B | KR(eta) | A: (c_inf - KR)/σ | chi2 B / A |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fitmoved = []
    for eta, f in sorted(fits.items()):
        P, B, A = f["P"], f["B"], f["A"]; dc = (A["cinf"] - B["cinf"]) / B["cinf_err"]; db = (A["b"] - B["b"]) / B["b_err"]
        if abs(dc) > 1 or abs(db) > 1: fitmoved.append((eta, round(dc, 2), round(db, 2)))
        print(f"| {eta:.2f} | {' '.join(str(int(n)) for n in A['Ns'])} | {P['cinf']:.5f} ± {P['cinf_err']:.5f} | {P['b']:+.4f} ± {P['b_err']:.4f} | "
              f"{B['cinf']:.5f} ± {B['cinf_err']:.5f} | {B['b']:+.4f} ± {B['b_err']:.4f} | {A['cinf']:.5f} ± {A['cinf_err']:.5f} | {A['b']:+.4f} ± {A['b_err']:.4f} | "
              f"{dc:+.2f} | {db:+.2f} | {A['kr']:.5f} | {(A['cinf'] - A['kr']) / A['cinf_err']:+.2f} | {B['chi2']:.3f} / {A['chi2']:.3f} |")
    print(f"\nfit parameters moving by more than 1 sigma (B): {fitmoved or 'none'}")
    reg = bool(moved or fitmoved)
    print(f"\n**VERDICT (pre-registered rule, sec. 14.3): {'REGENERATE' if reg else 'KEEP'}** -- "
          f"{len(moved)} point(s) with |change in D| > 0.5 (largest {maxd:.2f}); {len(fitmoved)} fit parameter(s) beyond 1 sigma.")
    print("P -> B (divider thickness, not in the verdict): c_inf moves by "
          + ", ".join(f"{e:.2f}: {(f['B']['cinf'] - f['P']['cinf']) / f['P']['cinf_err']:+.2f} σ" for e, f in sorted(fits.items())))
    return reg, allpts, fits

def write(allpts, fits):
    """sec. 14.3, REGENERATE branch: corrected per-mass tables and the corrected extrapolation table."""
    print("\n### Writing the corrected A2 tables (inputs: the _pre_boxtrunc_261002 copies)\n")
    for tab in TABLES:
        src = pre(tab); assert src != tab, "dated copy missing -- make it first (cp -n, cmp)"
        rows = list(csv.DictReader(open(T.plot_path(src)))); worst = 0.0
        for r in rows:
            eta, N, M = float(r["eta"]), int(r["N"]), int(r["M"]); g = allpts[tab][(eta, N)]
            nu = float(r["nu_mean"]); cB = nu / x_of(M, N, g["B"]["Le"])
            worst = max(worst, abs(cB - float(r["c_s_mass"])))          # gate: B geometry reproduces the table's column
            r.update(eta=f"{g['A']['eta']:.6f}", c_s_mass=f"{nu / x_of(M, N, g['A']['Le']):.6f}",
                     eta_rec=r["eta"], delta_sigma=f"{g['A']['d']:.6f}", L_eff_true=f"{g['A']['Le']:.6f}")
        if worst > 1.5e-6: sys.exit(f"{src}: B geometry does not reproduce c_s_mass (max diff {worst:.2e}) -- STOP")
        with open(T.plot_path(tab), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"  {tab}: {len(rows)} rows; B geometry reproduces the old c_s_mass to {worst:.1e}; "
              f"c_s_mass now at L_eff,true; eta = eta_true; new columns eta_rec, delta_sigma, L_eff_true")
    out = []
    for eta, f in sorted(fits.items()):
        A = f["A"]
        out.append(dict(eta=f"{eta:.2f}", N_values=" ".join(str(int(n)) for n in A["Ns"]), c_inf=f"{A['cinf']:.5f}",
                        c_inf_err=f"{A['cinf_err']:.5f}", slope=f"{A['b']:.5f}", slope_err=f"{A['b_err']:.5f}",
                        chi2=f"{A['chi2']:.3f}", dof=A["dof"], KR=f"{A['kr']:.5f}",
                        dev_KR_pct=f"{100 * (A['cinf'] - A['kr']) / A['kr']:+.3f}", dev_KR_sigma=f"{(A['cinf'] - A['kr']) / A['cinf_err']:+.2f}"))
    T.write_csv(T.plot_path(PUB), out)
    print(f"  {PUB}: {len(out)} rows, geometry A (L_eff = L0 - 2r - t/2 - delta/2, KR at eta_true per point, "
          f"c_inf = KR(eta) + intercept of the deviation fit); same columns as before")

def figure(allpts, fits):
    """##CHRIS 2026-10-02 (Task M4): 260917_A2_cs_vs_N redrawn in geometry A, the layout and colours of
    analyze_A2_X2p5_20260914.py:128-147 (which cannot be rerun for this: it recomputes everything at L_eff = L0 - 2r and
    would overwrite the corrected extrapolation table). Each point is c_s at L_eff,true, referred to the nominal eta by
    KR, y = c_s - KR(eta_true,N) + KR(eta) [DERIVATION: exactly the quantity the sec. 14.3 deviation fit uses, shifted by
    the constant KR(eta)], so the fitted line c_inf + b/sqrt(N) passes through the points as fitted. Gate before saving:
    the plotted c_inf, its error, b and its error equal the regenerated 260917 extrapolation table to its printed digits."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for ext in ("png", "pdf"):
        assert os.path.exists(T.plot_path(f"260917_A2_cs_vs_N_pre_boxtrunc_261002.{ext}")), "dated copy missing -- cp -n first"
    tab = {r["eta"]: r for r in csv.DictReader(open(T.plot_path(PUB)))}
    ok = all(tab[f"{e:.2f}"][k] == f"{f['A'][v]:.5f}" for e, f in fits.items()
             for k, v in (("c_inf", "cinf"), ("c_inf_err", "cinf_err"), ("slope", "b"), ("slope_err", "b_err")))
    print(f"\nfigure gate: plotted c_inf, b and errors equal {PUB} to its printed digits -> {'PASS' if ok else 'FAIL'}")
    if not ok: sys.exit("figure not drawn")
    src = allpts[TABLES[0]]; etas = sorted(fits); nc = min(3, len(etas)); nr = math.ceil(len(etas) / nc)
    fig, axs = plt.subplots(nr, nc, figsize=(4.6 * nc, 3.9 * nr), squeeze=False)
    for ax, eta in zip(axs.flat, etas):
        A = fits[eta]["A"]; Ns = A["Ns"]; u = 1 / np.sqrt(Ns)
        p = [src[(eta, int(N))]["A"] for N in Ns]
        y = np.array([q["cs"] - q["kr"] + A["kr"] for q in p]); s = np.array([max(q["sc"], 1e-12) for q in p])
        ax.errorbar(u, y, yerr=s, fmt="o", color="#2a78d6", capsize=3, ms=6, label="EDMD A2, per-N c_s ± mass scatter")
        uu = np.linspace(0, u.max() * 1.08, 50)
        ax.plot(uu, A["cinf"] + A["b"] * uu, "-", color="#52514e", lw=1.2, label=f"c_∞ + b/√N, χ² = {A['chi2']:.2f}")
        ax.errorbar([0], [A["cinf"]], yerr=[A["cinf_err"]], fmt="s", color="#eb6834", capsize=3, ms=7,
                    label=f"c_∞ = {A['cinf']:.4f} ± {A['cinf_err']:.4f}")
        ax.axhline(A["kr"], color="#e34948", lw=1.4, label=f"KR 2006 = {A['kr']:.4f}")
        for n, ui, ci in zip(Ns, u, y):
            ax.annotate(f"N={int(n)}", (ui, ci), textcoords="offset points", xytext=(4, 6), fontsize=7, color="#52514e")
        ax.set_title(f"η = {eta:.2f}  ({100*(A['cinf']-A['kr'])/A['kr']:+.2f} % vs KR, {(A['cinf']-A['kr'])/A['cinf_err']:+.1f}σ)", fontsize=10)
        ax.set_xlabel("1/√N"); ax.set_ylabel("c_s"); ax.grid(True, ls=":", alpha=0.6); ax.legend(fontsize=6.8, loc="best")
    for ax in list(axs.flat)[len(etas):]:
        ax.axis("off")
    fig.suptitle("A2, fixed η, N ladder: largest FFT bin at f ≥ ν_pred/2.5, mean over seeds, through-origin fit over masses\n"
                 "corrected for box truncation (methods §14)", fontsize=10.5)
    fig.text(0.99, 0.004, "each point: c_s at L_eff = L0 − 2r − t/2 − δ/2 and η_true, referred to the nominal η by KR: "
             "c_s − KR(η_true) + KR(η)\nfit on the deviations from KR (methods §14.3) · error bars = mass scatter (the fit weights) · "
             "data: 260919_A2_cs_per_mass.csv, 260917_A2_cs_vs_N_extrapolation.csv", ha="right", va="bottom", fontsize=7, color="0.4")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(T.plot_path(f"260917_A2_cs_vs_N.{ext}"), dpi=170)
    print("wrote 260917_A2_cs_vs_N.png/.pdf (geometry A)")

if __name__ == "__main__":
    reg, allpts, fits = main()
    if "--write" in sys.argv:
        if not reg: sys.exit("verdict KEEP -- nothing is written (sec. 14.3)")
        write(allpts, fits)
    if "--figure" in sys.argv:
        figure(allpts, fits)
