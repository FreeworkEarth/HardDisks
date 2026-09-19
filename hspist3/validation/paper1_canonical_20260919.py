#!/usr/bin/env python3
"""##CHRIS 2026-09-19: make the divider thickness canonical in Paper 1's figures.

L_eff = L0 - 2r - t/2 with t = 0.05 (tests_20260913.WALL_T). The estimator is linear in L_eff, so
every published c_s rescales by (L0 - 2r - t/2)/(L0 - 2r) exactly -- verified against the production
estimator on one cell recomputed from raw traces (eta = 0.5236, M = 1000, 25 seeds): the code path
and the closed form agree to 1.1e-16. That is why this rescales the CSVs instead of re-reading 7875
trajectories.

Writes 260919_* CSVs and figures. The 260914/260916/260917 versions are left alone.
"""
import csv, math, os, shutil, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tests_20260913 as T

P = T.PLOTS
RDISK, WALL_T = T.RDISK, T.WALL_T


def factor(L0):
    lo = L0 - 2 * RDISK
    return (lo - 0.5 * WALL_T) / lo


def rescale(src, dst, l0_col="L0", cols=("c_s", "c_s_err", "c_s_scatter_mass", "c_s_mass", "nu_mean", "nu_sd")):
    """Rescale every c_s-like column by the thickness factor of that row's own L0."""
    rows = list(csv.DictReader(open(os.path.join(P, src))))
    hit = [c for c in rows[0] if c in cols]
    for r in rows:
        f = factor(float(r[l0_col]))
        for c in hit:
            if c.startswith("nu"):
                continue          # frequencies are data, not derived -- leave them
            try:
                r[c] = f"{float(r[c]) * f:.6f}"
            except (ValueError, TypeError):
                pass
    with open(os.path.join(P, dst), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    return len(rows), hit


def draw_main():
    """The headline c_s(eta) figure, 260914 style, drawn from the corrected CSV."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plot_speed_of_sound_edmd as sos

    def cs_of(Z, e):
        h = 1e-5
        return sos.cs_adiabatic_2d_monatomic(Z(e), (Z(e + h) - Z(e - h)) / (2 * h), e, kbt=1, m=1)

    rows = list(csv.DictReader(open(os.path.join(P, "260919_A1v2_final_cs_vs_eta.csv"))))
    E = np.array([float(r["eta"]) for r in rows])
    C = np.array([float(r["c_s"]) for r in rows])
    S = np.array([float(r["c_s_scatter_mass"]) for r in rows])
    fig, ax = plt.subplots(figsize=(10.0, 6.4))
    sos.add_eta_regime_shading(ax, x_max=0.78)
    for t_ in list(ax.texts):
        if t_.get_rotation() == 90 and t_.get_position()[0] < 0.05:
            t_.set_position((t_.get_position()[0], 0.42))
    e = np.linspace(0.001, 0.78, 600); ek = e[e <= 0.69]
    ax.plot(e, cs_of(sos.Z_spt_eos, e), "--", color="#eda100", lw=1.4, label="SPT equation of state")
    ax.plot(e, cs_of(sos.Z_henderson_eos, e), "-.", color="#1baf7a", lw=1.4, label="Henderson (a = 0.125)")
    ax.plot(ek, cs_of(sos.Z_kolafa_rottner_2006, ek), "-", color="#e34948", lw=2.4,
            label="Kolafa-Rottner 2006 (valid to eta = 0.69)")
    ax.axhline(math.sqrt(2), color="0.45", lw=1.0, ls=":", label="ideal-gas limit: c_s = sqrt(2)")
    ax.errorbar(E, C, yerr=S, fmt="o-", color="#2a78d6", ms=4.5, lw=1.1, capsize=2.5, zorder=5,
                label="A1 v2, N = 100, 9 masses x 25 seeds, 200 periods,\nlargest FFT bin at f >= nu_pred/2.5, "
                      "L_eff = L0 - 2r - t/2")
    ax.set_ylim(0, 22); ax.set_xlim(0, 0.78)
    ax.set_xlabel("Packing fraction  eta", fontsize=11.5)
    ax.set_ylabel("Speed of sound  c_s  [sqrt(k_BT/m)]", fontsize=11.5)
    ax.set_title("Speed of sound against packing fraction (canonical, divider thickness in L_eff)", fontsize=12.5)
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=8.5, framealpha=0.95)
    ax.text(0.555, 2.2, "eta >= 0.65: 6 sigma compartment, structure\nchanges during measurement - not a\nfluid-branch value",
            fontsize=8.5, color="0.3")
    fig.text(0.99, 0.004, "L_eff = L0 - 2r - t/2 with t = 0.05 sigma - data: 260919_A1v2_final_cs_vs_eta.csv",
             ha="right", va="bottom", fontsize=7, color="0.4")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = os.path.join(P, "260919_cs_vs_eta")
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    print(f"  main figure -> {os.path.basename(out)}.png/.pdf")


def main():
    print("### Making L_eff = L0 - 2r - t/2 canonical (t = %.3f)\n" % WALL_T)
    for src, dst in (("260914_A1v2_final_cs_vs_eta.csv", "260919_A1v2_final_cs_vs_eta.csv"),
                     ("260917_A2_cs_per_mass.csv", "260919_A2_cs_per_mass.csv"),
                     ("260916_A2_cs_per_mass.csv", "260919_A2_cs_per_mass_famB.csv")):
        if not os.path.exists(os.path.join(P, src)):
            print(f"  {src}: MISSING, skipped"); continue
        n, hit = rescale(src, dst)
        print(f"  {src} -> {dst}: {n} rows, rescaled columns {hit}")

    # the deviation column in the A1 CSV is derived; recompute it rather than leave it stale
    import numpy as np
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import plot_speed_of_sound_edmd as sos
    p = os.path.join(P, "260919_A1v2_final_cs_vs_eta.csv")
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        e = float(r["eta"])
        if e <= 0.69:
            h = 1e-5; Z = sos.Z_kolafa_rottner_2006
            kr = float(sos.cs_adiabatic_2d_monatomic(Z(e), (Z(e + h) - Z(e - h)) / (2 * h), e, kbt=1, m=1))
            r["KR"] = f"{kr:.5f}"
            r["dev_KR_pct"] = f"{100 * (float(r['c_s']) - kr) / kr:+.3f}"
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"  deviation column recomputed in {os.path.basename(p)}")

    print("\n### Regenerating the three circulated figures from the corrected CSVs\n")
    here = os.path.dirname(os.path.abspath(__file__))
    env = dict(os.environ, HD_A1_CSV="260919_A1v2_final_cs_vs_eta.csv")
    jobs = [
        ("dilute zoom  -> 260919_cs_vs_eta_lowdensity_zoom",
         [sys.executable, os.path.join(here, "lowdensity_zoom_20260917.py"),
          "260919_A2_cs_per_mass.csv", "260919_cs_vs_eta_lowdensity_zoom"], env),
        ("N100 vs A2   -> 260919_cs_vs_eta_N100_vs_A2",
         [sys.executable, os.path.join(here, "overlay_N100_vs_A2_20260915.py")],
         dict(env, HD_A2_CSV="260919_A2_cs_per_mass_famB.csv", HD_OUT="260919_cs_vs_eta_N100_vs_A2")),
    ]
    for name, cmd, e in jobs:
        r = subprocess.run(cmd, cwd=here, capture_output=True, text=True, env=e)
        tail = (r.stdout or r.stderr).strip().splitlines()[-1:] or ["(no output)"]
        print(f"  {name}: {'ok' if r.returncode == 0 else 'FAILED'} -- {tail[0][-90:]}")
    draw_main()


if __name__ == "__main__":
    main()
