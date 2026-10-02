#!/usr/bin/env python3
"""##CHRIS 2026-10-12: linewidth decomposition of the divider mode on the 261006 mode ladder.

    Gamma_meas = gamma_M / M_eff(alpha) + c * omega^2 + Gamma_res

Analysis only, on 261006_mode_ladder.json (seven cells, v1 binary, contraction ON -- no post-flag
data enters). Two stages, run in this order and never the other way round:

    --prereg   prints every fixed input and prediction (no fit). Its output is pasted into the
               pre-registration in 260913_tests_REPORT.md, which is committed BEFORE --fit runs.
    --fit      refuses to run unless the committed HEAD copy of that report already contains the
               pre-registration marker; then fits ONCE and prints the verdict by the written rule.

CONVENTION (260912_paper1_methods.md sec. 13): Gamma is the ENERGY-decay rate, Gamma = 2/tau_r,
where tau_r is the AMPLITUDE decay time of C(t) = A e^{-t/tau_T} + B e^{-t/tau_r} cos(omega t).
Gamma is also the FWHM of the line in angular frequency; Delta f_FWHM = Gamma/(2 pi).

SOURCES, read from the PDF (Mansour, Garcia & Baras, PRE 73, 016121 (2006)):
  p. 3, Eqs. 8-11: Enskog eta_s, zeta, kappa, g2 for hard disks (they cite ref. [20]: Gass,
       J. Chem. Phys. 54, 1898 (1971), and Barker & Henderson, RMP 48, 587 (1976)). The Gass PDF
       is NOT in the repo, so no Gass page is cited; the coefficients are Mansour's transcription.
  p. 5, Eq. 15: linear velocity profile; Eq. 17: M_hat Xp'' = L_y(P_L - P_R)
       - L_y v_p (Gamma_L/X_p + Gamma_R/(L_x - X_p)), Gamma = zeta + eta;  Eq. 18: M_hat = M + mN/3.
  => friction coefficient gamma_M = 2 L_y (eta_s + zeta) / X_p at the midpoint; the oscillator's
     energy-decay rate is gamma_M / M_hat. M_hat is exactly the K -> 0 limit of M_eff below.
"""
import json, math, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import plot_speed_of_sound_edmd as sos
REPO = os.path.dirname(os.path.dirname(HERE))
LADDER = os.path.join(HERE, "261006_mode_ladder.json")
REPORT_REL = "0000_PLAN_OVERALL/ALL_MARKDOWNS/260913_tests_REPORT.md"
MARKER = "PRE-REGISTRATION 1b: linewidth decomposition"
FIGDIR = os.path.join(REPO, "0000_PLAN_OVERALL", "paper1_speedofsound", "experiments", "final")

ETA = 0.10134170                       # audited, free compartment length, both boxes (261006)
H = 10.0                               # L_y
# published in 261006 sec. 3 (Mansour/measured, 2 decimals) -- the reproduction target
PUBLISHED_RATIO = [2.04, 1.72, 1.33, 1.14, 1.03, 1.87, 1.52]

def enskog(eta, T=1.0):
    """Mansour Eqs. 8-11 (p. 3), sigma = m = k_B = 1; n is the number density."""
    n = 4.0 * eta / math.pi
    g2 = (1 - 7 * n * math.pi / 64) / (1 - n * math.pi / 4) ** 2
    png = math.pi * n * g2
    eta_s = 0.2555 * n * math.sqrt(math.pi) * (1 + 2 / png + 0.4365 * png) * math.sqrt(T)
    zeta = 0.1592 * math.pi ** 1.5 * n * n * g2 * math.sqrt(T)
    kappa = 1.029 * n * math.sqrt(math.pi) * (1.5 + 2 / png + 0.4359 * png) * math.sqrt(T)
    return n, g2, eta_s, zeta, kappa

def meff_gas(K, NS):
    """Gas inertia of the standing wave, both compartments: 2 N_s m [1/2 - sin2K/(4K)] / sin^2 K."""
    return 2.0 * NS * (0.5 - math.sin(2 * K) / (4 * K)) / math.sin(K) ** 2

def g_visc(K):
    """Standing-wave viscous dissipation relative to the linear profile of Eq. 15 (secondary only)."""
    return K * K * (0.5 + math.sin(2 * K) / (4 * K)) / math.sin(K) ** 2

def inputs():
    J = json.load(open(LADDER)); cells = J["cells"]
    Z = float(sos.Z_kolafa_rottner_2006(np.array([ETA]))[0])
    dZ = float(sos.dZ_kolafa_rottner_2006(np.array([ETA]))[0])
    cs = math.sqrt(Z + ETA * dZ + Z * Z); gam = 1 + Z * Z / (Z + ETA * dZ)
    n, g2, es, ze, ka = enskog(ETA)
    rho, cv = n, 1.0; cp = gam * cv
    visc = (es + ze) / rho; therm = (gam - 1) * ka / (rho * cp)
    pred = dict(c_th=therm / cs ** 2, c_full=(visc + therm) / cs ** 2, c_visc=visc / cs ** 2)
    rows = []
    for r in cells:
        K = r["K_pred"]; NS = r["NS"]; LE = r["LEFF"]; LC = r["LC"]
        gM = 2 * H * (es + ze) / LE
        Me = r["M"] + meff_gas(K, NS)
        G = 2.0 / r["taur_x"]; sG = 2.0 * r["taur_err_x"] / r["taur_x"] ** 2
        om = 2 * math.pi / r["per_x"]; som = 2 * math.pi * r["per_err_x"] / r["per_x"] ** 2
        rows.append(dict(box=r["box"], NS=NS, LC=LC, LE=LE, M=r["M"], al=r["alpha"], K=K,
                         Me=Me, Mhat=r["Mhat"], gM=gM, fixed=gM / Me, G=G, sG=sG, om=om, som=som,
                         om2=om * om, taur=r["taur_x"], staur=r["taur_err_x"],
                         ratio=r["taur_man"] / r["taur_x"],
                         fixed_LC=2 * H * (es + ze) / LC / Me, fixed_gK=gM * g_visc(K) / Me))
    return dict(Z=Z, dZ=dZ, cs=cs, gam=gam, n=n, g2=g2, es=es, ze=ze, ka=ka, rho=rho, cp=cp,
                visc=visc, therm=therm, pred=pred, rows=rows)

def prereg():
    I = inputs(); P = I["pred"]
    a = enskog(0.100)
    print("### R3. Enskog transcription check against the 260918 audit (eta = 0.100)\n")
    print(f"eta_s = {a[2]:.4f}, zeta = {a[3]:.4f}, eta_s + zeta = {a[2]+a[3]:.4f}  (audit: 0.314, 0.017, 0.331)\n")
    print(f"### Fixed inputs at the box packing fraction eta = {ETA:.8f}\n")
    print("| quantity | value | source |\n|---|---|---|")
    print(f"| n = rho (m = 1) | {I['n']:.6f} | 4 eta / pi |")
    print(f"| g2 | {I['g2']:.5f} | Mansour Eq. 11, p. 3 |")
    print(f"| eta_s | {I['es']:.5f} | Mansour Eq. 8, p. 3 |")
    print(f"| zeta | {I['ze']:.5f} | Mansour Eq. 9, p. 3 |")
    print(f"| kappa | {I['ka']:.5f} | Mansour Eq. 10, p. 3 |")
    print(f"| Z, eta Z' | {I['Z']:.6f}, {ETA*I['dZ']:.6f} | Kolafa-Rottner 2006 |")
    print(f"| c_s | {I['cs']:.6f} | (Z + eta Z' + Z^2)^(1/2) |")
    print(f"| gamma = c_P/c_V | {I['gam']:.5f} | 1 + Z^2/(Z + eta Z') |")
    print(f"| c_v, c_p per unit mass | 1, {I['cp']:.5f} | 2D hard disks: c_v = k_B/m exactly |")
    print(f"| (eta_s + zeta)/rho | {I['visc']:.5f} | viscous diffusivity |")
    print(f"| (gamma-1) kappa/(rho c_p) | {I['therm']:.5f} | thermal part |")
    print(f"| **c_th** = (gamma-1) kappa/(rho c_p c_s^2) | **{P['c_th']:.4f}** | primary prediction |")
    print(f"| c_full = [(eta_s+zeta)/rho + (gamma-1) kappa/(rho c_p)]/c_s^2 | {P['c_full']:.4f} | literal free-wave form, secondary |")
    print(f"| c_visc = (eta_s+zeta)/(rho c_s^2) | {P['c_visc']:.4f} | the part already inside gamma_M |")
    print("\n### Per-cell inputs (tau_r, period: 261006_mode_ladder.json, x-trace fit; errors: delete-8 block jackknife)\n")
    print("| box | N_s | L_eff | M | alpha | K | M_hat | M_eff | gamma_M | gamma_M/M_eff | tau_r | Gamma_meas = 2/tau_r | omega^2 | c_th omega^2 | Mansour/measured (261006 form) | published |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    ok = True
    for r, pub in zip(I["rows"], PUBLISHED_RATIO):
        good = abs(round(r["ratio"], 2) - pub) < 1e-9; ok &= good
        print(f"| {r['box']} | {r['NS']} | {r['LE']:.2f} | {r['M']} | {r['al']:.3f} | {r['K']:.4f} | {r['Mhat']:.2f} | "
              f"{r['Me']:.2f} | {r['gM']:.5f} | {r['fixed']:.3e} | {r['taur']:.0f} ± {r['staur']:.0f} | "
              f"{r['G']:.4e} ± {r['sG']:.1e} | {r['om2']:.4e} | {P['c_th']*r['om2']:.3e} | {r['ratio']:.2f} | {pub:.2f} {'OK' if good else 'MISMATCH'} |")
    print(f"\n**R1 (Mansour/measured reproduced from the JSON to 2 decimals): {'PASS' if ok else 'FAIL'}**\n")
    A = {r["M"]: r for r in I["rows"] if r["box"] == "A"}; B = {r["M"]: r for r in I["rows"] if r["box"] == "B"}
    print("**R2 (tau_r proportional to L at fixed M):**\n")
    for M in (50, 100):
        q = B[M]["taur"] / A[M]["taur"]
        sq = q * math.hypot(B[M]["staur"] / B[M]["taur"], A[M]["staur"] / A[M]["taur"])
        print(f"- M = {M}: tau_r(B)/tau_r(A) = {q:.3f} ± {sq:.3f}; L_c ratio {B[M]['LC']/A[M]['LC']:.3f}, "
              f"L_eff ratio {B[M]['LE']/A[M]['LE']:.3f}; (q - 2)/sigma = {(q-2)/sq:+.1f}")
    print("\nsigma_omega contributes at most "
          f"{max(2*r['som']/r['om']*r['om2']*P['c_th']/r['sG'] for r in I['rows']):.3f} of sigma_Gamma per cell, so omega^2 is treated as exact.")
    return ok

def head_has_marker():
    try:
        txt = subprocess.run(["git", "-C", REPO, "show", f"HEAD:{REPORT_REL}"], capture_output=True,
                             text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return False
    return MARKER in txt

def wls(x, y, s):
    w = 1 / s ** 2; X = np.vstack([x, np.ones_like(x)]).T
    C = np.linalg.inv(X.T @ (X * w[:, None])); p = C @ (X.T @ (w * y))
    chi2 = float(np.sum(w * (y - X @ p) ** 2)); return p, C, chi2

def fit():
    if not head_has_marker():
        sys.exit(f"REFUSED: the committed HEAD copy of {REPORT_REL} has no '{MARKER}'. Commit the pre-registration first.")
    from scipy.stats import chi2 as chi2d
    I = inputs(); P = I["pred"]; R = I["rows"]
    x = np.array([r["om2"] for r in R]); s = np.array([r["sG"] for r in R])
    G = np.array([r["G"] for r in R]); fx = np.array([r["fixed"] for r in R])
    (c, Gres), C, chi2 = wls(x, G - fx, s); dof = len(R) - 2
    sc, sr = math.sqrt(C[0, 0]), math.sqrt(C[1, 1]); scale = max(1.0, math.sqrt(chi2 / dof))
    pval = float(chi2d.sf(chi2, dof))
    print(f"### Fit (run once): Gamma_meas - gamma_M/M_eff = c omega^2 + Gamma_res, weighted, {len(R)} cells, {dof} dof\n")
    print(f"c = {c:.4f} ± {sc:.4f} (raw) ± {sc*scale:.4f} (x sqrt(chi2_red) = {scale:.3f})")
    print(f"Gamma_res = {Gres:.3e} ± {sr:.1e} (raw) ± {sr*scale:.1e} (scaled)")
    print(f"chi2 = {chi2:.2f} on {dof} dof, chi2_red = {chi2/dof:.2f}, p = {pval:.3g}"
          f"{'  -> MODEL FORM REJECTED (p < 0.01) by the pre-registered flag' if pval < 0.01 else ''}")
    corr = C[0, 1] / (sc * sr); print(f"corr(c, Gamma_res) = {corr:+.3f}\n")
    se = sc * scale
    dE = (c - P["c_th"]) / se; d0 = c / se
    print(f"Verdict quantities (scaled sigma_c = {se:.4f}): (c - c_th)/sigma = {dE:+.2f};  c/sigma = {d0:+.2f};  "
          f"(c - c_full)/sigma = {(c-P['c_full'])/se:+.2f} [secondary]\n")
    nearE, near0 = abs(dE) <= 2, abs(d0) <= 2
    if nearE and not near0: verdict = "OUTCOME 1 -- c consistent with Enskog (c_th), and distinguishable from 0"
    elif near0 and not nearE: verdict = "OUTCOME 2 -- c = 0 within errors, and inconsistent with Enskog (c_th)"
    elif nearE and near0: verdict = "OUTCOME 4 (overlap) -- data cannot separate c_th from 0"
    else: verdict = "OUTCOME 3 -- neither: report Gamma_res vs L"
    print(f"**{verdict}**\n")
    print("Per-cell residual against the PREDICTION (no fitted parameter): r = Gamma_meas - gamma_M/M_eff - c_th omega^2\n")
    print("| box | L_c | M | alpha | Gamma_meas | gamma_M/M_eff | c_th omega^2 | r | r/sigma | fit residual/sigma |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in R:
        rr = r["G"] - r["fixed"] - P["c_th"] * r["om2"]
        fr = (r["G"] - r["fixed"] - c * r["om2"] - Gres) / r["sG"]
        print(f"| {r['box']} | {r['LC']:.2f} | {r['M']} | {r['al']:.3f} | {r['G']:.4e} | {r['fixed']:.3e} | "
              f"{P['c_th']*r['om2']:.3e} | {rr:+.3e} | {rr/r['sG']:+.1f} | {fr:+.1f} |")
    print("\nGamma_res vs L (weighted mean of r per box, scaled by sqrt(chi2_red) of the box mean where > 1):\n")
    for b in ("A", "B"):
        rr = np.array([r["G"] - r["fixed"] - P["c_th"] * r["om2"] for r in R if r["box"] == b])
        ss = np.array([r["sG"] for r in R if r["box"] == b]); w = 1 / ss ** 2
        m = float(np.sum(w * rr) / np.sum(w)); sm = math.sqrt(1 / np.sum(w))
        ch = float(np.sum(w * (rr - m) ** 2)); dd = len(rr) - 1
        L = [r["LC"] for r in R if r["box"] == b][0]
        print(f"- box {b} (L_c = {L}): Gamma_res = {m:+.3e} ± {sm*max(1,math.sqrt(ch/max(dd,1))):.1e}  "
              f"(chi2 = {ch:.1f} on {dd} dof; Gamma_res * L_c = {m*L:+.4f}, Gamma_res * L_c^2 = {m*L*L:+.3f})")
    print("\n### Secondary (reported, no verdict)\n")
    for lab, key in (("gamma_M with L_c instead of L_eff", "fixed_LC"),
                     ("standing-wave viscous profile, gamma_M g(K)/M_eff", "fixed_gK")):
        f2 = np.array([r[key] for r in R]); (c2, g2r), C2, ch2 = wls(x, G - f2, s)
        sc2 = math.sqrt(C2[0, 0]) * max(1, math.sqrt(ch2 / dof))
        print(f"- {lab}: c = {c2:.4f} ± {sc2:.4f}, Gamma_res = {g2r:.3e}, chi2 = {ch2:.2f}/{dof}")
    figure(R, P, c, Gres, se, s)
    return dict(c=c, sc=sc, se=se, Gres=Gres, sr=sr, chi2=chi2, dof=dof, p=pval, verdict=verdict)

def figure(R, P, c, Gres, se, s):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for b, mk in (("A", "o"), ("B", "s")):
        rs = [r for r in R if r["box"] == b]
        ax.errorbar([r["om2"] * 1e3 for r in rs], [(r["G"] - r["fixed"]) * 1e3 for r in rs],
                    yerr=[r["sG"] * 1e3 for r in rs], fmt=mk, color="tab:blue", ms=6, capsize=3,
                    label=f"box {b}, L_c = {rs[0]['LC']}", mfc="tab:blue" if b == "A" else "white")
    xx = np.linspace(0, max(r["om2"] for r in R) * 1.08, 50)
    ax.plot(xx * 1e3, (c * xx + Gres) * 1e3, color="tab:blue", lw=1.5, label=f"fit: c = {c:.2f} ± {se:.2f}")
    ax.plot(xx * 1e3, P["c_th"] * xx * 1e3, color="red", lw=1.5, label=f"Enskog thermal, c_th = {P['c_th']:.2f}")
    ax.plot(xx * 1e3, P["c_full"] * xx * 1e3, color="red", lw=1.2, ls="--", label=f"Enskog free-wave, c_full = {P['c_full']:.2f}")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel(r"$\omega^2$  [$10^{-3}\,\tau^{-2}$]"); ax.set_ylabel(r"$\Gamma_{\rm meas}-\gamma_M/M_{\rm eff}$  [$10^{-3}\,\tau^{-1}$]")
    ax.set_title("Divider-mode linewidth beyond Mansour friction (261006 ladder, v1)", fontsize=10)
    ax.legend(fontsize=8, frameon=False); fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"261012_paper1_linewidth_decomp.{ext}"), dpi=200)
    print(f"\nfigure: {os.path.relpath(os.path.join(FIGDIR, '261012_paper1_linewidth_decomp.png'), REPO)} (+ .pdf)")

if __name__ == "__main__":
    if "--prereg" in sys.argv: sys.exit(0 if prereg() else 1)
    elif "--fit" in sys.argv: fit()
    else: sys.exit("usage: --prereg | --fit")
