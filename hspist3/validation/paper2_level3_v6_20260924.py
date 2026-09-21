#!/usr/bin/env python3
"""##CHRIS 2026-09-24: Level 3 v6 -- close the settled-state comparison with the BOX's own inputs.

v5 established that the wall's settled position is the right observable (the "transient excess" was
a peak statistic, reproduced by a no-push control) and that the bulk-EOS model predicts it to a few
per cent. This closes the comparison two ways:

  1. WINDOW VALIDITY. The settled average needs the wall to have completed several free periods,
     T_w = 2 pi sqrt(M_s/(k + k_gas)) = 60 sigma-time at M_s = 50 and 120 at 200. A cell whose
     record leaves less than 3 periods after t > tau + 3 T_w is not a settled measurement, it is a
     phase sample; v5 reported one such cell (u = 0.1, M_s = 200, 0.5 period) as a 1.7 sigma result
     and it should have been flagged. Valid cells are pooled by inverse variance.

  2. BOX INPUTS INSTEAD OF BULK INPUTS. The model is parameter-free but its inputs are bulk: the
     Kolafa-Rottner equation of state and the adiabat it implies. The box is N = 100, and Paper 1
     measures a 100-disk divider box as ~1 % stiffer in c_s than bulk KR, i.e. ~2 % in c_s^2 = the
     gas spring. Level 1 measures the gas ending 0.5 % hotter than the quasi-static adiabat. Both
     are MEASUREMENTS with error bars, not fits, and both feed the same fixed point. Nothing here
     is adjusted to improve agreement.

A logical point that v5 got wrong and is corrected here: the settled displacement solves
k s = F(L - dx + s) - F(L), which depends only on the DIFFERENCE of the force. A constant offset in
F cancels exactly. So the +1.3 % excess standing force measured on the held wall CANNOT explain the
settled residual -- it explains the wall's pre-push offset, which the control already removes. v5
offered it as "a third to a half" of the residual; that was wrong.
"""
import glob, math, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plot_speed_of_sound_edmd as sos

ET = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/hspist3/"
      "experiments_energy_transfer")
V3, V4, V6 = f"{ET}/level3_master_preload_20260921", f"{ET}/level3_v4_20260922", f"{ET}/level3_v6_20260924"
P1CSV = ("/Users/chrisharing/Desktop/CCS_complex_coupled_systems/Repo/HardDisks/0000_PLAN_OVERALL/"
         "paper1_speedofsound/experiments/final/260919_A1v2_final_cs_vs_eta.csv")

N, L, H, R, DX, K = 100, 78.5, 10.0, 0.5, 7.96, 0.5
ETA = N * math.pi * R * R / (L * H)
Z = float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0])
DZ = float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
CS_BULK = math.sqrt(Z + ETA * DZ + Z * Z)
KGAS = N * CS_BULK ** 2 / L ** 2
F_EOS = N * Z / L
F_BOX_NOM = 1.5961                      # A1 held-wall, pre-signal window
TF_ADIABAT, TF_LEVEL1 = 1.1432, 1.149   # adiabat vs Level 1's measurement


def adiabat_force(cs_scale=1.0, tf_target=None, f_scale=1.0):
    """Tabulate F(L) on the adiabat. cs_scale scales the gas STIFFNESS via c_s (k_gas ~ c_s^2);
    tf_target rescales the adiabat exponent so T(L-dx)/T(L) hits a measured value.

    NOTE on the held-wall force. F_box exceeds the bulk-EOS F by 1.35 %, but that excess is applied
    as an ADDITIVE offset, not a multiplicative scaling, and therefore cancels from the fixed point.
    The reason is Paper 1's own result: the wall pressure Z_wall is a SURFACE CONTACT value, not the
    bulk stiffness -- measured at Paper 1's densities it over-predicts the c_s offset five-fold, and
    the claim that it explained the N = 100 offset was withdrawn on 2026-09-19. Scaling the whole
    force curve by F_box/F_EOS (and hence its slope, hence k_gas) would repeat exactly that error.
    The additive excess still matters in the ENERGY, dE = F s + k s^2/2, where F is the actual
    standing force at the operating point."""
    LG = np.linspace(L - DX - 10.0, L + 10.0, 4001)
    Zs = np.array([float(sos.Z_kolafa_rottner_2006(np.array([N * math.pi * R * R / (x * H)]))[0])
                   for x in LG])
    lnT = np.concatenate([[0.0], np.cumsum(-0.5 * (Zs[1:] + Zs[:-1]) * np.diff(np.log(LG / L)))])
    lnT = lnT - np.interp(L, LG, lnT)
    if tf_target is not None:
        cur = float(np.interp(L - DX, LG, lnT))
        lnT = lnT * (math.log(tf_target) / cur)
    FG = N * np.exp(lnT) * Zs / LG + (f_scale - 1.0) * F_EOS   # ADDITIVE, see docstring
    # stiffness scaling: steepen the force curve about L without moving F(L)
    if cs_scale != 1.0:
        F0 = float(np.interp(L, LG, FG))
        FG = F0 + (FG - F0) * cs_scale ** 2
    return LG, FG


def fixed_point(LG, FG):
    F0 = float(np.interp(L, LG, FG))
    lo, hi = 0.0, 6.0
    for _ in range(80):
        s = 0.5 * (lo + hi)
        if K * s < float(np.interp(L - DX + s, LG, FG)) - F0: lo = s
        else: hi = s
    return s


def dE(s, F=F_EOS):
    return F * s + 0.5 * K * s * s


def seeds(pat):
    out, tt = [], None
    for f in sorted(glob.glob(pat)):
        d = pd.read_csv(f, low_memory=False, usecols=["Time", "W0_x_sigma", "PistonR_v"])
        t = d["Time"].to_numpy(float); x = d["W0_x_sigma"].to_numpy(float)
        mv = np.nonzero(np.abs(d["PistonR_v"].to_numpy(float)) > 1e-12)[0]
        i0 = mv[0] if len(mv) else 0
        out.append(-(x[i0:] - x[i0])); ti = t[i0:] - t[i0]
        if tt is None or len(ti) < len(tt): tt = ti
    if not out: return None, None
    n = min(len(a) for a in out)
    return tt[:n], np.array([a[:n] for a in out])


def main():
    TW = {m: 2 * math.pi / math.sqrt((K + KGAS) / m) for m in (50, 200)}
    print("## Level 3 v6 -- window validity and the box-input model\n")
    print(f"eta = {ETA:.4f}, c_s(bulk KR) = {CS_BULK:.4f}, k_gas = {KGAS:.5f}, "
          f"k + k_gas = {K+KGAS:.4f}")
    print(f"T_w = 2 pi sqrt(M/(k+k_gas)): {TW[50]:.0f} sigma-time at M_s = 50, "
          f"{TW[200]:.0f} at 200\n")

    ctl = {m: seeds(f"{V4}/A4_free_M{m}/tr_*.csv") for m in (50, 200)}

    print("### 1. Window validity and pooling\n")
    print("| u | M_s | T_w | window start | record end | periods in window | verdict | s̄ corrected |")
    print("|---|---|---|---|---|---|---|---|")
    cells = [(0.02, 50, V4, "A3_M50_u0.02"), (0.02, 200, V4, "A3_M200_u0.02"),
             (0.05, 50, V3, "k0.5_M50_u0.05"), (0.05, 200, V6, "k0.5_M200_u0.05"),
             (0.1, 50, V3, "k0.5_M50_u0.1"), (0.1, 200, V6, "k0.5_M200_u0.1")]
    valid = []
    for u, ms, src, pat in cells:
        t, S = seeds(f"{src}/{pat}/tr_*.csv")
        if t is None:
            print(f"| {u} | {ms} | {TW[ms]:.0f} | — | — | — | NO DATA | — |"); continue
        tau = DX / u; t0 = tau + 3 * TW[ms]
        per_win = (t[-1] - t0) / TW[ms]
        if per_win < 3:
            print(f"| {u} | {ms} | {TW[ms]:.0f} | {t0:.0f} | {t[-1]:.0f} | {per_win:.1f} | "
                  f"**INVALID** (<3) | — |")
            continue
        w = t >= t0
        perseed = S[:, w].mean(axis=1)
        tc, Sc = ctl[ms]
        wc = (tc >= 0) & (tc <= t[-1] - t0)
        perc = Sc[:, wc].mean(axis=1)
        sb, sbe = perseed.mean(), perseed.std(ddof=1) / math.sqrt(len(perseed))
        sc, sce = perc.mean(), perc.std(ddof=1) / math.sqrt(len(perc))
        s_corr, e_corr = sb - sc, math.hypot(sbe, sce)
        valid.append((u, ms, s_corr, e_corr))
        print(f"| {u} | {ms} | {TW[ms]:.0f} | {t0:.0f} | {t[-1]:.0f} | {per_win:.1f} | valid | "
              f"{s_corr:.3f} ± {e_corr:.3f} |")

    w_ = np.array([1 / e ** 2 for _, _, _, e in valid])
    sp = float(np.sum(w_ * np.array([s for _, _, s, _ in valid])) / np.sum(w_))
    spe = float(1 / math.sqrt(np.sum(w_)))
    print(f"\n  pooled over {len(valid)} valid cells (inverse variance): "
          f"s̄ = **{sp:.4f} ± {spe:.4f} sigma**")

    print("\n### 2. The fixed point with bulk inputs, then with the box's own\n")
    # sensitivities
    dlnS_dlnK = K / (K + KGAS)
    print(f"  sensitivities:  dln s/dln k_gas = k/(k+k_gas) = {dlnS_dlnK:.3f};")
    print(f"                  dln s/dln F for an ADDITIVE offset = 0 exactly -- the fixed point")
    print(f"                  k s = F(L-dx+s) - F(L) depends only on the DIFFERENCE. A")
    print(f"                  MULTIPLICATIVE scaling is a different thing: it scales the slope")
    print(f"                  too, and would give dln s = 0.910 x dln F. The held-wall excess is")
    print(f"                  treated as additive because Z_wall is a contact value (Paper 1).\n")

    # Paper 1 box stiffening, interpolated to this eta
    rows = [(float(r["eta"]), float(r["c_s"]), float(r["c_s_scatter_mass"])) for r in
            __import__("csv").DictReader(open(P1CSV))]
    rows.sort()
    es = np.array([r[0] for r in rows]); cs = np.array([r[1] for r in rows])
    ce = np.array([r[2] for r in rows])          # c_s_scatter_mass; c_s_err is all zeros in the CSV
    krs = np.array([float(sos.cs_adiabatic_2d_monatomic(
        sos.Z_kolafa_rottner_2006(np.array([e])), sos.dZ_kolafa_rottner_2006(np.array([e])),
        np.array([e]), kbt=1, m=1)[0]) for e in es])
    dev = cs / krs
    r_box = float(np.interp(ETA, es, dev)); r_err = float(np.interp(ETA, es, ce / cs))
    print(f"  Paper 1 at eta = {ETA:.4f} (interpolated between the two bracketing densities):")
    print(f"    c_s(box)/c_s(bulk KR) = {r_box:.4f} +- {r_err:.4f}  -> k_gas ratio "
          f"{r_box**2:.4f} +- {2*r_err:.4f}")
    print(f"  A1 held wall: F_box = {F_BOX_NOM:.4f} vs F_EOS = {F_EOS:.4f} "
          f"({100*(F_BOX_NOM/F_EOS-1):+.2f} %)")
    print(f"  Level 1: T_f/T_i = {TF_LEVEL1} vs adiabat {TF_ADIABAT} "
          f"({100*(TF_LEVEL1/TF_ADIABAT-1):+.2f} %)\n")

    print("| model inputs | s_qs | vs pooled s̄ | σ | ΔE(s_qs) [kT] |")
    print("|---|---|---|---|---|")
    variants = [
        ("bulk KR (v3/v4/v5 model)", dict()),
        ("(a) + uniform F_box scaling", dict(f_scale=F_BOX_NOM / F_EOS)),
        ("(a)+(b) + Paper 1 c_s,box", dict(f_scale=F_BOX_NOM / F_EOS, cs_scale=r_box)),
        ("(a)+(b)+(c) + Level 1 adiabat", dict(f_scale=F_BOX_NOM / F_EOS, cs_scale=r_box,
                                               tf_target=TF_LEVEL1)),
    ]
    out = {}
    for lab, kw in variants:
        LG, FG = adiabat_force(**kw)
        s = fixed_point(LG, FG)
        out[lab] = s
        d = 100 * (sp - s) / s
        nsig = abs(sp - s) / spe
        Fuse = F_BOX_NOM if "(a)" in lab else F_EOS
        print(f"| {lab} | {s:.4f} | {d:+.2f} % | {nsig:.1f} | {dE(s, Fuse):.4f} |")

    # error propagation on the full box model
    LG, FG = adiabat_force(f_scale=F_BOX_NOM / F_EOS, cs_scale=r_box * (1 + r_err),
                           tf_target=TF_LEVEL1)
    s_hi = fixed_point(LG, FG)
    LG, FG = adiabat_force(f_scale=F_BOX_NOM / F_EOS, cs_scale=r_box * (1 - r_err),
                           tf_target=TF_LEVEL1)
    s_lo = fixed_point(LG, FG)
    s_full = out["(a)+(b)+(c) + Level 1 adiabat"]
    mod_err = 0.5 * abs(s_hi - s_lo)
    tot = math.hypot(spe, mod_err)
    print(f"\n  box-input model: s_qs = {s_full:.4f} +- {mod_err:.4f} (from Paper 1's c_s error)")
    print(f"  measurement    : s̄   = {sp:.4f} +- {spe:.4f}")
    print(f"  residual       : {100*(sp-s_full)/s_full:+.2f} +- {100*tot/s_full:.2f} %  "
          f"-> **{abs(sp-s_full)/tot:.1f} sigma**")

    print("\n### 3. Consistency check on the pre-push offset\n")
    off = (F_BOX_NOM - F_EOS) / (K + KGAS)
    cm = []
    for m in (50, 200):
        tc, Sc = ctl[m]
        w = (tc >= 0) & (tc <= 800)
        cm.append(Sc[:, w].mean())
    print(f"  predicted from the held-wall force: (F_box - F_EOS)/(k + k_gas) = "
          f"{(F_BOX_NOM-F_EOS):.4f}/{K+KGAS:.4f} = **{off:.4f} sigma**")
    print(f"  measured, A4 no-push control mean : **{np.mean(cm):.4f} sigma** "
          f"(M_s = 50: {cm[0]:.4f}, 200: {cm[1]:.4f})")
    print(f"  two independent measurements of the same +1.3 % excess standing force.")
    print(f"  NOTE this offset does NOT enter the settled residual: it is removed by the control,")
    print(f"  and a uniform F offset cancels from the fixed point anyway (section 2).")


if __name__ == "__main__":
    main()
