#!/usr/bin/env python3
"""##CHRIS 2026-10-09 (261012 sec. 4.7.4, decision 3 item 2): the Kolafa-Rottner range table for Paper 1, printed by script.
For each packing-fraction range: which reference applies, the fit accuracy AS STATED IN THE SOURCE (quoted with page), the
source's own data points, the number of our data points per Paper 1 table, and what Paper 1 does there (quoted from the draft).

[SOURCE] J. Kolafa and M. Rottner, Mol. Phys. 104, 3435-3441 (2006), DOI 10.1080/00268970600967963; PDF on disk at
ZZZ_PAPER/Kolafa and Rottner - 2006 - Simulation-based equation of state of the hard disk fluid ... .pdf (journal pages 3435-3441
are PDF pages 2-8; the quotes were read from the page images and from pdftotext). rho = N sigma^2 / A, eta = (pi/4) rho.
[DECISION, methods sec. 15/15.1, 2026-10-02] Paper 1 uses the rho_max = 0.90 fit (Eq. 7, sec. 3.2), fitted to eta <= 0.7069
(= pi 0.90 / 4), and compares data with it only for eta <= 0.69 (plot_speed_of_sound_edmd.KR2006_PLOT_ETA_MAX).
[DERIVATION] the fit-to-fit spread of c_s among the paper's three published fits (rho_max 0.88, 0.89, 0.90), c_s^2 = Z + eta Z'
+ Z^2 (kT = m = 1), with the coefficients typed from the page image of p. 3438-3439 (below). The project's sanity script
validation/paper1_kr_sanity_261002.py types the x^57 coefficient of the rho_max = 0.89 fit as 5.77730095e-23; the paper prints
5.57730095e-23 (p. 3439). Its effect is printed below; Paper 1 uses the 0.90 fit, whose coefficients agree.
usage (from hspist3/):  python3 validation/paper1_kr_ranges_261009.py
"""
import csv, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); HS = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, HS)
import plot_speed_of_sound_edmd as SOS
import paper1_kr_sanity_261002 as SAN

ROOT = os.path.dirname(HS)
FINAL = os.path.join(ROOT, "0000_PLAN_OVERALL", "paper1_speedofsound", "experiments", "final")
DRAFT = os.path.join(ROOT, "0000_PLAN_OVERALL", "paper1_speedofsound", "writeup", "paper1_draft.tex")
ETA_FIT, ETA_CMP = SOS.KR2006_ETA_MAX, SOS.KR2006_PLOT_ETA_MAX

# the three published fits, typed from the page image (pp. 3438-3439): power of x -> coefficient
KR_PDF = {
    0.88: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504951668, 6: -0.0325433846, 7: 0.0133946531,
           8: 0.00174265604, 9: -0.00944632202, 10: 0.00851111768, 11: -0.0035963525, 12: 0.000577345106, 19: -1.06399127e-7},
    0.89: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504963915, 6: -0.0325578581, 7: 0.0134816028,
           8: 0.00129187484, 9: -0.00808881628, 10: 0.00669011963, 11: -0.00250795961, 12: 0.000336036442, 22: -5.15282664e-9,
           57: 5.57730095e-23},
    0.90: {0: 1.0, 1: 2.0, 2: 1.12801775, 3: 0.00181895291, 4: -0.0526134737, 5: 0.0504960168, 6: -0.0325537792, 7: 0.0134578632,
           8: 0.00140888182, 9: -0.00834273601, 10: 0.00694127367, 11: -0.00262254723, 12: 0.000355746352, 22: -5.24672938e-9,
           57: 5.88054639e-23},
}
# Table 1 of the paper (p. 3436): rho, Z, sigma(Z); finite-size corrected "for rho < 0.88 by equation (5) and for rho > 0.88 by
# Z(1/N) linear extrapolation"; footnote a (rho 0.88): "A compromise between 10.23095(26) by (5) and 10.23055(54) by Z(1/N)
# extrapolation."; footnote b (rho 0.90): "This value may be affected by finite-size effects."
KR_TABLE1 = [(0.40, 2.1514393, 0.0000052), (0.45, 2.4276680, 0.0000065), (0.50, 2.7601235, 0.0000079), (0.55, 3.1647878, 0.0000098),
             (0.60, 3.663691, 0.000012), (0.65, 4.287926, 0.000015), (0.70, 5.082362, 0.000019), (0.75, 6.113391, 0.000026),
             (0.80, 7.476491, 0.000036), (0.83, 8.494891, 0.000050), (0.84, 8.866011, 0.000059), (0.85, 9.245785, 0.000072),
             (0.86, 9.621609, 0.000082), (0.87, 9.96782, 0.00013), (0.88, 10.2309, 0.0003), (0.89, 10.3176, 0.0011), (0.90, 10.2059, 0.0011)]

QUOTE = {
    "s": ("p. 3438", "\"The value of s for an optimum fit is around unity provided that the input standard errors sigma are reliable, "
                     "which is the case for our simulations where sigma is determined with an accuracy (error of the error) of a few percent [21].\""),
    "fit90": ("p. 3439", "\"rho_max = 0.90, s = 0.927; region rho in [0.89, 0.90] of this equation may be affected by finite-size effects\""),
    "corr": ("p. 3438", "\"Correction term (5) is not applicable for rho >= 0.89 because it contains the second derivative of the EOS and "
                        "therefore the correction term is large and not available with sufficient precision. Linear extrapolation of Z(1/N) was "
                        "used instead; the final results thus lack precision.\""),
    "mc": ("p. 3438", "\"The data in the 'difficult' region close to the phase transition agree well with recent extensive Monte Carlo data [3, 4] "
                      "with the exception of density rho = 0.9 closest to rho_c where the N = 1024^2 result Z = 10.212 [3] is significantly "
                      "larger than our Z = 10.206\""),
    "extra": ("p. 3439", "\"Any extrapolation to rho > rho_c should be done with caution because function p(rho) is likely to be non-analytical at rho_c.\""),
    "loop": ("p. 3439", "\"both equations with rho_max >= 0.89 predict to some extent the loop (with the 'classical' critical exponent alpha' = 3) "
                        "at the critical (fluid/hexatic) point, even if this is not the aim of the present work which focuses rather on the "
                        "low-density region.\""),
}


def Z(eta, c):
    x = eta / (1 - eta); return sum(a * x ** i for i, a in c.items())


def dZ(eta, c):
    x = eta / (1 - eta); return sum(i * a * x ** (i - 1) for i, a in c.items() if i) / (1 - eta) ** 2


def cs(eta, c):
    z, d = Z(eta, c), dZ(eta, c); v = z + eta * d + z * z
    return math.sqrt(v) if v > 0 else float("nan")


def rows(name):
    p = os.path.join(FINAL, name)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def bucket(e):
    return 0 if e <= ETA_CMP else (1 if e <= ETA_FIT else 2)


def draft_line(n):
    return open(DRAFT).read().splitlines()[n - 1].strip()


def main():
    print("# The Kolafa-Rottner range table for Paper 1 (261012 sec. 4.7.4 decision 3 item 2), printed by validation/paper1_kr_ranges_261009.py\n")
    # ---------------------------------------------------------------- the fit itself: the module against the page
    worst = max(abs(SOS.Z_kolafa_rottner_2006([e])[0] / Z(e, KR_PDF[0.90]) - 1) for e in (0.1, 0.3, 0.5, 0.65, 0.69, 0.7069))
    print(f"The module's rho_max = 0.90 fit against the coefficients typed here from the page: largest relative difference in Z {worst:.1e} "
          f"(eta 0.1-0.7069). Fit range eta <= pi 0.90/4 = {ETA_FIT:.4f}; comparison range eta <= {ETA_CMP:.2f} (rho <= {4 * ETA_CMP / math.pi:.4f}).\n")
    # ---------------------------------------------------------------- the source's data per range
    names = ["eta <= 0.69 (compared)", "0.69 < eta <= 0.7069 (fit range, not compared)", "eta > 0.7069 (beyond the fit)"]
    src = {0: [], 1: [], 2: []}
    for rho, z, s in KR_TABLE1:
        src[bucket(math.pi * rho / 4)].append((rho, z, s))
    # ---------------------------------------------------------------- our data per range, per Paper 1 table
    ours = []
    a1 = rows("260919_A1v2_final_cs_vs_eta.csv")
    c = [0, 0, 0]; cr = [0, 0, 0]
    for r in a1: c[bucket(float(r["eta"]))] += 1; cr[bucket(float(r["eta_rec"]))] += 1
    ours.append(("A1 v2 c_s(eta), N = 100 (260919_A1v2_final_cs_vs_eta.csv), by eta_true", c))
    ours.append(("the same, by the recorded (nominal) eta", cr))
    for name, lab, key in (("260917_A2_cs_vs_N_extrapolation.csv", "A2 N -> infinity extrapolation (eta)", "eta"),
                           ("260916_A2_finite_size_forms.csv", "A2 finite-size forms (distinct eta)", "eta"),
                           ("261004_p1_confinement_cells.csv", "confinement cells (eta_true)", "eta_true"),
                           ("261005_p1_identity_afix_cells.csv", "identity cells, A-fixed (eta_lab; no eta_true column)", "eta_lab")):
        R = rows(name); es = sorted({float(r[key]) for r in R}); c = [0, 0, 0]
        for e in es: c[bucket(e)] += 1
        ours.append((lab, c))
    for name, lab in (("260919_A2_cs_per_mass.csv", "A2 per mass, (eta, N) state points"), ("260919_A2_cs_per_mass_famB.csv", "A2 famB per mass, (eta, N) state points")):
        R = rows(name); pts = sorted({(float(r["eta"]), int(float(r["N"]))) for r in R}); c = [0, 0, 0]
        for e, N in pts: c[bucket(e)] += 1
        ours.append((lab, c))
    # ---------------------------------------------------------------- the table
    print("## 1. The ranges\n")
    print("| range | rho = 4 eta / pi | reference that applies | fit accuracy as stated in the source (quoted) | the source's MD data in the range "
          "(Table 1, p. 3436: rho, Z, sigma(Z)) | our data points (A1 v2 c_s(eta), by eta_true / by recorded eta) | what Paper 1 does (draft) |\n"
          "|---|---|---|---|---|---|---|")
    ref = ["Kolafa-Rottner 2006, rho_max = 0.90 fit (Eq. 7, sec. 3.2), fitted to eta <= 0.7069: compared with data",
           "the same fit, inside its fit range but NOT compared (decision of 2026-10-02); the comparison is Engel et al. 2013 (qualitative)",
           "none of the KR fits (extrapolation); Engel et al. 2013 landmarks; a global EOS with hexatic and solid branches (Liu 2021) is OPEN"]
    acc = [f"{QUOTE['s'][1]} ({QUOTE['s'][0]}); {QUOTE['fit90'][1].split(';')[0]}\" ({QUOTE['fit90'][0]}). The source states no accuracy for Z', which c_s uses.",
           f"{QUOTE['fit90'][1]} ({QUOTE['fit90'][0]}); {QUOTE['corr'][1]} ({QUOTE['corr'][0]}); {QUOTE['mc'][1]} ({QUOTE['mc'][0]})",
           f"{QUOTE['extra'][1]} ({QUOTE['extra'][0]}); {QUOTE['loop'][1]} ({QUOTE['loop'][0]})"]
    rho_r = [f"<= {4 * ETA_CMP / math.pi:.4f}", f"({4 * ETA_CMP / math.pi:.4f}, 0.90]", "> 0.90"]
    a1c, a1r = ours[0][1], ours[1][1]
    paper = [f"compares (dev_KR_pct column of the table; draft:126-129: \"{draft_line(128)} {draft_line(129)}\"); the mass-independence test averages "
             f"\"the 24 densities with $\\eta \\leq 0.69$\" (draft:248), counted by recorded eta",
             f"shows the points, \"not counted as a deviation from it\" (draft:263-264); KR \"dashed from there to $\\eta = 0.705$, inside its fit range\" "
             f"(draft:476-477); the comparison is with Engel et al., qualitative (draft:433-466)",
             f"no KR curve: \"{draft_line(478)} {draft_line(479).split('.')[0]}.\" (draft:478-479); Engel et al.'s landmarks; a quantitative comparison "
             f"\"needs a global equation of state with hexatic and solid branches\" (draft:503)"]
    for k in range(3):
        sd = src[k]
        st = "; ".join(f"{rho:.2f}: {z} +- {s} (rel. {s / z:.1e})" for rho, z, s in sd) if sd else "none"
        print(f"| {names[k]} | {rho_r[k]} | {ref[k]} | {acc[k]} | {len(sd)} points{': ' + st if len(sd) <= 3 else f' (rho {sd[0][0]:.2f}-{sd[-1][0]:.2f}; relative sigma(Z) <= {max(s / z for _, z, s in sd):.1e})'} | "
              f"{a1c[k]} / {a1r[k]} | {paper[k]} |")
    print(f"\nThe boundary: one A1 v2 row, eta_true = 0.690460, has recorded eta 0.689999. It carries a KR deviation in the table (+43.191 %) and is "
          f"one of the draft's '24 densities with eta <= 0.69': the comparison range is applied by the recorded eta. By eta_true it lies in the "
          f"second range. Above: by eta_true / by recorded eta.")
    print("\n## 2. Our data points per Paper 1 table and range\n")
    print("| table | eta <= 0.69 | 0.69 < eta <= 0.7069 | eta > 0.7069 |\n|---|---|---|---|")
    for lab, c in ours:
        print(f"| {lab} | {c[0]} | {c[1]} | {c[2]} |")
    # ---------------------------------------------------------------- the fit-to-fit spread
    print("\n## 3. [DERIVATION] c_s from the paper's three published fits (kT = m = 1), and the 0.89 coefficient in the sanity script\n")
    print("| eta | rho | c_s, rho_max 0.88 | c_s, 0.89 | c_s, 0.90 (Paper 1) | largest relative difference to the 0.90 fit | c_s, 0.89 with the sanity script's x^57 coefficient | its relative change |\n"
          "|---|---|---|---|---|---|---|---|")
    for e in (0.40, 0.55, 0.60, 0.65, 0.67, 0.68, 0.69, 0.695, 0.70, 0.7069):
        v = {k: cs(e, KR_PDF[k]) for k in KR_PDF}
        dmax = max(abs(v[k] / v[0.90] - 1) for k in (0.88, 0.89)) if v[0.90] == v[0.90] else float("nan")
        vs = cs(e, SAN.KR[0.89])
        print(f"| {e:.4f} | {4 * e / math.pi:.4f} | {v[0.88]:.6f} | {v[0.89]:.6f} | {v[0.90]:.6f} | {dmax:.1e} | {vs:.6f} | {abs(vs / v[0.89] - 1):.1e} |")
    print(f"\n(the sanity script's rho_max = 0.89 x^57 coefficient: {SAN.KR[0.89][57]:.8e}; the paper, p. 3439: {KR_PDF[0.89][57]:.8e}. "
          f"Paper 1 uses the 0.90 fit: the sanity script's 0.90 coefficients equal the page's: "
          f"{'yes' if SAN.KR[0.90] == KR_PDF[0.90] else 'NO'}; its 0.88 coefficients: {'yes' if SAN.KR[0.88] == KR_PDF[0.88] else 'NO'}.)")


if __name__ == "__main__":
    main()
