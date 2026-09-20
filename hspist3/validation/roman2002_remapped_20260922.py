#!/usr/bin/env python3
"""##CHRIS 2026-09-22: B1 -- Roman 2002 Table I re-mapped, against KR and against our N = 100.

Roman compares his measured c_s to SPT and Henderson through the gamma = 2 mapping
c_s^2 = 2 (kT/m) Z_something (his Eq. 26), and finds an excess that GROWS with density: +1.3 % at
eta = 0.112 rising to +12.8 % at 0.524. That growth is the mapping, not his data. The correct
adiabatic mapping for a 2D monatomic hard-disk fluid is

    c_s^2 = (kT/m) [ Z + eta Z' + Z^2 ] ,

and re-mapping his own numbers through it flattens the excess to about +1 % at every density --
the same flat offset our N = 100 divider box shows against the same equation of state.

Two independent codes, 24 years apart, same box (N_s = 50 per side, H = 10, same observable): that
is a much stronger statement for the paper than "we agree with Roman". It says the +1 % is a
property of the N = 100 divider box, not of either implementation.

Analysis only -- no runs, nothing overwritten. Roman's numbers are typed from his Table I (p. 850)
and his eta is RECOMPUTED here from his own L_0 rather than taken from his table, so a typo in
either column would show up as a mismatch.
"""
import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos

OUT = os.path.join(
    "/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL",
    "paper1_speedofsound/experiments/final")

# Roman 2002 Table I, p. 850. (L_0, c_s, sigma_c_s). His box: N_0 = 100 total, N_s = 50 per side,
# A = H = 10, sigma = 1 so r = 0.5. eta is recomputed below from L_0, not copied.
ROMAN = [(7.5, 5.99, 0.09), (10.0, 3.78, 0.08), (15.0, 2.61, 0.03), (20.0, 2.20, 0.02),
         (25.0, 2.01, 0.02), (30.0, 1.89, 0.02), (35.0, 1.81, 0.02)]
# Roman's own SPT / Henderson columns, for the transcription check
ROMAN_SPT = [5.31, 3.53, 2.50, 2.15, 1.97, 1.86, 1.78]
ROMAN_HEN = [5.45, 3.59, 2.53, 2.16, 1.97, 1.86, 1.79]
NS, H, R = 50, 10.0, 0.5


def eta_of(L0):
    """Roman's packing fraction: N_s disks of radius r in one compartment L_0 x H."""
    return NS * math.pi * R * R / (L0 * H)


def adiabatic(Zf, dZf, eta):
    e = np.array([eta], float)
    return float(sos.cs_adiabatic_2d_monatomic(Zf(e), dZf(e), e, kbt=1.0, m=1.0)[0])


def gamma2(Zf, dZf, eta):
    """Roman's Eq. (26): c_s^2 = 2 (kT/m) (Z + eta Z'), the constant-gamma=2 assumption."""
    e = np.array([eta], float)
    return float(sos.cs_roman_gamma2_from_eos(Zf(e), dZf(e), e, kbt=1.0, m=1.0)[0])


def main():
    print("## B1 -- Roman 2002 Table I re-mapped\n")
    print("Transcription check: his own SPT / Henderson columns against this repo's EOS, gamma = 2\n")
    print("| L0 | eta (recomputed) | his SPT | ours | his Henderson | ours |")
    print("|---|---|---|---|---|---|")
    for (L0, _, _), sp, he in zip(ROMAN, ROMAN_SPT, ROMAN_HEN):
        e = eta_of(L0)
        print(f"| {L0:g} | {e:.4f} | {sp:.2f} | {gamma2(sos.Z_spt_eos, sos.dZ_spt_eos, e):.2f} | "
              f"{he:.2f} | {gamma2(sos.Z_henderson_eos, sos.dZ_henderson_eos, e):.2f} |")

    print("\n### Roman against Kolafa-Rottner, both mappings\n")
    print("| L0 | eta | c_s Roman | his +- [%] | KR gamma=2 | dev | KR adiabatic | dev |")
    print("|---|---|---|---|---|---|---|---|")
    rows = []
    for L0, cs, dcs in ROMAN:
        e = eta_of(L0)
        g2 = gamma2(sos.Z_kolafa_rottner_2006, sos.dZ_kolafa_rottner_2006, e)
        ad = adiabatic(sos.Z_kolafa_rottner_2006, sos.dZ_kolafa_rottner_2006, e)
        rows.append((e, L0, cs, dcs, g2, ad))
        print(f"| {L0:g} | {e:.4f} | {cs:.2f} | {100*dcs/cs:.1f} | {g2:.3f} | "
              f"{100*(cs-g2)/g2:+.1f} % | {ad:.3f} | {100*(cs-ad)/ad:+.1f} % |")
    dev_ad = [100 * (c - a) / a for _, _, c, _, _, a in rows]
    dev_g2 = [100 * (c - g) / g for _, _, c, _, g, _ in rows]
    print(f"\n  adiabatic mapping: mean {np.mean(dev_ad):+.2f} %, spread "
          f"{np.min(dev_ad):+.2f} to {np.max(dev_ad):+.2f} %  -- FLAT")
    print(f"  his gamma=2      : {np.min(dev_g2):+.2f} to {np.max(dev_g2):+.2f} %  -- grows with eta")

    # ---- our canonical N = 100 curve, same EOS, same mapping
    ours = []
    path = T.plot_path("260919_A1v2_final_cs_vs_eta.csv")
    for r in csv.DictReader(open(path)):
        e, c = float(r["eta"]), float(r["c_s"])
        s = float(r.get("c_s_scatter_mass", 0.0) or 0.0)
        if e <= 0.69:
            ad = adiabatic(sos.Z_kolafa_rottner_2006, sos.dZ_kolafa_rottner_2006, e)
            ours.append((e, c, s, ad, 100 * (c - ad) / ad))
    lo = [d for e, _, _, _, d in ours if e <= 0.4]
    print(f"\n  ours, N = 100, eta <= 0.4: mean {np.mean(lo):+.2f} % above KR adiabatic "
          f"({len(lo)} densities)  [{os.path.basename(path)}]")

    # the one place the two disagree
    e52 = min(rows, key=lambda r: abs(r[0] - 0.5236))
    o52 = min(ours, key=lambda o: abs(o[0] - 0.5236))
    dr = 100 * (e52[2] - e52[5]) / e52[5]
    sr = 100 * e52[3] / e52[2]
    print(f"\n  eta ~ 0.52: Roman {dr:+.1f} +- {sr:.1f} %, ours {o52[4]:+.1f} % "
          f"-> {abs(o52[4]-dr)/sr:.1f} sigma apart on his error bar alone")

    # ---- figure
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(8.6, 7.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2.1, 1]})
    e = np.linspace(0.05, 0.60, 400)
    ax.plot(e, sos.cs_adiabatic_2d_monatomic(sos.Z_kolafa_rottner_2006(e),
                                             sos.dZ_kolafa_rottner_2006(e), e, kbt=1, m=1),
            "-", color="#e34948", lw=2.3, label="Kolafa–Rottner, adiabatic mapping")
    ax.plot(e, sos.cs_roman_gamma2_from_eos(sos.Z_kolafa_rottner_2006(e), sos.dZ_kolafa_rottner_2006(e), e, kbt=1, m=1),
            "--", color="#e34948", lw=1.5, alpha=0.85, label="Kolafa–Rottner, Román's $\\gamma=2$ mapping")
    ax.errorbar([r[0] for r in rows], [r[2] for r in rows], yerr=[r[3] for r in rows],
                fmt="s", color="black", ms=6, capsize=3, lw=0, elinewidth=1.2,
                label="Román 2002, Table I ($N=100$, 100 traj.)")
    E = np.array([o[0] for o in ours]); C = np.array([o[1] for o in ours]); S = np.array([o[2] for o in ours])
    k = E <= 0.60
    ax.errorbar(E[k], C[k], yerr=S[k], fmt="o", color="#2a78d6", ms=4, capsize=2, lw=0,
                elinewidth=1.0, label="this work, $N=100$ (25 seeds, 9 masses)")
    ax.set_ylabel("Speed of sound  $c_s$  [$\\sqrt{k_BT/m}$]", fontsize=11)
    ax.set_ylim(1.5, 6.6); ax.grid(True, ls=":", alpha=0.6)
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.95)
    ax.set_title("Román 2002 re-mapped: the density-dependent excess was the mapping", fontsize=12)

    axd.axhline(0, color="#e34948", lw=1.6)
    axd.axhspan(-1, 1, color="#e1e0d9", alpha=0.7, zorder=0)
    axd.errorbar([r[0] for r in rows], dev_ad,
                 yerr=[100 * r[3] / r[5] for r in rows], fmt="s", color="black", ms=6,
                 capsize=3, lw=0, elinewidth=1.2, label="Román, adiabatic mapping")
    axd.plot([r[0] for r in rows], dev_g2, "s--", color="0.55", ms=5, lw=1.0,
             label="Román, his own $\\gamma=2$ mapping")
    KRa = np.array([adiabatic(sos.Z_kolafa_rottner_2006, sos.dZ_kolafa_rottner_2006, x) for x in E[k]])
    axd.errorbar(E[k], 100 * (C[k] - KRa) / KRa, yerr=100 * S[k] / KRa, fmt="o", color="#2a78d6",
                 ms=4, capsize=2, lw=0, elinewidth=1.0, label="this work, $N=100$")
    axd.set_xlabel("Packing fraction  $\\eta$", fontsize=11)
    axd.set_ylabel("deviation from KR [%]", fontsize=10)
    axd.set_xlim(0.05, 0.60); axd.set_ylim(-3, 14)
    axd.grid(True, ls=":", alpha=0.6); axd.legend(fontsize=8.5, loc="upper left", framealpha=0.95)
    axd.text(0.055, -2.4, "±1 % band", fontsize=8, color="0.35")
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    base = os.path.join(OUT, "260922_roman2002_remapped_vs_KR")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", dpi=200)
    with open(f"{base}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["L0", "eta", "cs_roman", "err_roman", "KR_gamma2", "dev_gamma2_pct",
                    "KR_adiabatic", "dev_adiabatic_pct"])
        for (e_, L0, cs, dcs, g2, ad) in rows:
            w.writerow([L0, f"{e_:.6f}", cs, dcs, f"{g2:.4f}", f"{100*(cs-g2)/g2:.3f}",
                        f"{ad:.4f}", f"{100*(cs-ad)/ad:.3f}"])
    print(f"\nwritten {base}.png/.pdf/.csv")


if __name__ == "__main__":
    main()
